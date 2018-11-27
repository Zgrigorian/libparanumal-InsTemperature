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

#include "ins.h"
#include "mesh.h"

void insHeatGradient(ins_t *ins, occa::memory o_U, occa::memory o_T){
   
   mesh_t *mesh = ins->mesh;
   
   if(mesh->totalHaloPairs>0){
   ins->heatHaloExtractKernel(mesh->Nelements,
                              mesh->totalHaloPairs,
                              mesh->o_haloElementList,
                              ins->fieldOffset,
                              o_U,
                              o_T,
                              ins->o_hHaloBuffer);

   //copy exracted halo to HOST
   ins->o_hHaloBuffer.copyTo(ins->hSendBuffer);

   // start halo exchange
   meshHaloExchangeStart(mesh,
                         mesh->Np*(ins->NVfields+1)*sizeof(dfloat),
                         ins->hSendBuffer,
                         ins->hRecvBuffer);
   }      
   // Compute Volume Contribution
   ins->heatGradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_Dmatrices,
                            ins->o_T,
                            ins->o_Tx,
                            ins->o_Ty);
   if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
    
      ins->o_hHaloBuffer.copyFrom(ins->hRecvBuffer);
      
      ins->heatHaloScatterKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 ins->fieldOffset,
                                 o_U,
                                 o_T,
                                 ins->o_hHaloBuffer);
   }
   
   ins->heatGradientSurfaceKernel(mesh->Nelements,
                                  mesh->o_sgeo,
                                  mesh->o_LIFTT,
                                  mesh->o_vmapM,
                                  mesh->o_vmapP,
                                  mesh->o_EToB,
                                  ins->o_T,
                                  ins->o_Tx,
                                  ins->o_Ty);
}
