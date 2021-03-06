/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "elliptic.h"

void MGLevel::Ax(occa::memory o_x, occa::memory o_Ax) {
  ellipticOperator(elliptic,lambda,
                    o_x,o_Ax, dfloatString); // "float" ); // hard coded for testing (should make an option)
}

void MGLevel::residual(occa::memory o_rhs, occa::memory o_x, occa::memory o_res) {
  ellipticOperator(elliptic,lambda,
                    o_x,o_res, dfloatString); // "float" ); // hard coded for testing (should make an option)

  // subtract r = b - A*x
  ellipticScaledAdd(elliptic, 1.f, o_rhs, -1.f, o_res);
}

void MGLevel::coarsen(occa::memory o_x, occa::memory o_Rx) {
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    elliptic->dotMultiplyKernel(mesh->Nelements*NpF, o_invDegree, o_x, o_x);

  elliptic->precon->coarsenKernel(mesh->Nelements, o_R, o_x, o_Rx);

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ogsGatherScatter(o_Rx, ogsDfloat, ogsAdd, elliptic->ogs);
    if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Rx);
  }
}

void MGLevel::prolongate(occa::memory o_x, occa::memory o_Px) {
  elliptic->precon->prolongateKernel(mesh->Nelements, o_R, o_x, o_Px);
}

void MGLevel::smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero) {
  if (stype==RICHARDSON) {
    this->smoothRichardson(o_rhs, o_x, x_is_zero);
  } else if (stype==CHEBYSHEV) {
    this->smoothChebyshev(o_rhs, o_x, x_is_zero);
  }
}

void MGLevel::smoother(occa::memory o_x, occa::memory o_Sx) {
  if (smtype==JACOBI) {
    this->smootherJacobi(o_x, o_Sx);
  } else if (smtype==LOCALPATCH) {
    this->smootherLocalPatch(o_x, o_Sx);
  }
}

void MGLevel::smoothRichardson(occa::memory &o_r, occa::memory &o_x, bool xIsZero) {

  occa::memory o_res = o_smootherResidual;

  if (xIsZero) {
    this->smoother(o_r, o_x);
    return;
  }

  dfloat one = 1.; dfloat mone = -1.;

  //res = r-Ax
  this->Ax(o_x,o_res);
  elliptic->scaledAddKernel(Nrows, one, o_r, mone, o_res);

  //smooth the fine problem x = x + S(r-Ax)
  this->smoother(o_res, o_res);
  elliptic->scaledAddKernel(Nrows, one, o_res, one, o_x);
}

void MGLevel::smoothChebyshev (occa::memory &o_r, occa::memory &o_x, bool xIsZero) {

  const dfloat theta = 0.5*(lambda1+lambda0);
  const dfloat delta = 0.5*(lambda1-lambda0);
  const dfloat invTheta = 1.0/theta;
  const dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  dfloat one = 1., mone = -1., zero = 0.0;

  occa::memory o_res = o_smootherResidual;
  occa::memory o_Ad  = o_smootherResidual2;
  occa::memory o_d   = o_smootherUpdate;

  if(xIsZero){ //skip the Ax if x is zero
    //res = Sr
    this->smoother(o_r, o_res);

    //d = invTheta*res
    elliptic->scaledAddKernel(Nrows, invTheta, o_res, zero, o_d);
  } else {
    //res = S(r-Ax)
    this->Ax(o_x,o_res);
    elliptic->scaledAddKernel(Nrows, one, o_r, mone, o_res);
    this->smoother(o_res, o_res);

    //d = invTheta*res
    elliptic->scaledAddKernel(Nrows, invTheta, o_res, zero, o_d);
  }

  for (int k=0;k<ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    if (xIsZero&&(k==0))
      elliptic->scaledAddKernel(Nrows, one, o_d, zero, o_x);
    else
      elliptic->scaledAddKernel(Nrows, one, o_d, one, o_x);

    //r_k+1 = r_k - SAd_k
    this->Ax(o_d,o_Ad);
    this->smoother(o_Ad, o_Ad);
    elliptic->scaledAddKernel(Nrows, mone, o_Ad, one, o_res);

    rho_np1 = 1.0/(2.*sigma-rho_n);
    dfloat rhoDivDelta = 2.0*rho_np1/delta;

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    elliptic->scaledAddKernel(Nrows, rhoDivDelta, o_res, rho_np1*rho_n, o_d);

    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  elliptic->scaledAddKernel(Nrows, one, o_d, one, o_x);
}

void MGLevel::smootherLocalPatch(occa::memory &o_r, occa::memory &o_Sr) {

  // occaTimerTic(mesh->device,"approxBlockJacobiSolveKernel");
  elliptic->precon->approxBlockJacobiSolverKernel(mesh->Nelements,
                            elliptic->precon->o_patchesIndex,
                            elliptic->precon->o_invAP,
                            elliptic->precon->o_invDegreeAP,
                            o_r,
                            o_Sr);
  // occaTimerToc(mesh->device,"approxBlockJacobiSolveKernel");
}

void MGLevel::smootherJacobi(occa::memory &o_r, occa::memory &o_Sr) {
  elliptic->dotMultiplyKernel(mesh->Np*mesh->Nelements,o_invDiagA,o_r,o_Sr);
}

