#ifndef PARALMOND_H
#define PARALMOND_H 1

typedef struct csr_t {

  int Nrows;
  int Ncols;

  int NlocalCols;

  //local
  int diagNNZ;
  int   *diagRowStarts;
  int   *diagCols;
  dfloat *diagCoefs;

  //non-local
  int offdNNZ;
  int   *offdRowStarts;
  int   *offdCols;
  dfloat *offdCoefs;

  dfloat *diagInv;

  dfloat *null;

  //storage for smoothing
  dfloat *scratch;

  int   *colMap;

  // MPI halo exchange info
  int  NHalo;
  int  NrecvTotal;  // number of elements to be sent in halo exchange
  int  NsendTotal;
  int  totalHaloPairs;
  int *haloElementList; // sorted list of elements to be sent in halo exchange
  int *NsendPairs;      // number of elements worth of data to send
  int *NrecvPairs;      // number of elements worth of data to recv
  int  NsendMessages;   // number of messages to send
  int  NrecvMessages;   // number of messages to recv
  dfloat *sendBuffer;

  void *haloSendRequests;
  void *haloRecvRequests;

} csr;


typedef struct ell_t {

  int Nrows;
  int Ncols;
  int nnzPerRow;
  int strideLength;
  int actualNNZ;

  occa::memory o_cols;
  occa::memory o_coefs;

} ell;

typedef struct coo_t {

  int Nrows;
  int Ncols;
  int nnz;

  // device memory
  occa::memory o_offsets;
  occa::memory o_cols;
  occa::memory o_coefs;

} coo;

typedef struct hyb_t {

  int Nrows;
  int Ncols;

  int NlocalCols;

  coo *C;
  ell *E;

  occa::memory o_diagInv;

  occa::memory o_null;

  // MPI halo exchange info
  int  NHalo;
  int *colMap;
  int  NrecvTotal;  // number of elements to be sent in halo exchange
  int  NsendTotal;
  int *haloElementList; // sorted list of elements to be sent in halo exchange
  occa::memory o_haloElementList;
  int *NsendPairs;      // number of elements worth of data to send
  int *NrecvPairs;      // number of elements worth of data to recv
  int  NsendMessages;   // number of messages to send
  int  NrecvMessages;   // number of messages to recv
  dfloat   *sendBuffer;
  dfloat   *recvBuffer;
  occa::memory o_haloBuffer;

  void *haloSendRequests;
  void *haloRecvRequests;

} hyb;


typedef struct dcsr_t {

  int Nrows;
  int Ncols;

  int NlocalCols;

  //local
  int diagNNZ;
  occa::memory o_diagRows;
  occa::memory o_diagCols;
  occa::memory o_diagCoefs;

  //non-local
  int offdNNZ;
  occa::memory o_offdRows;
  occa::memory o_offdCols;
  occa::memory o_offdCoefs;

  // MPI halo exchange info
  int  NHalo;
  int  NrecvTotal;  // number of elements to be sent in halo exchange
  int  NsendTotal;
  int  totalHaloPairs;
  int *haloElementList; // sorted list of elements to be sent in halo exchange
  int *NsendPairs;      // number of elements worth of data to send
  int *NrecvPairs;      // number of elements worth of data to recv
  int  NsendMessages;   // number of messages to send
  int  NrecvMessages;   // number of messages to recv
  dfloat   *sendBuffer;
  dfloat   *recvBuffer;

  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;

  void *haloSendRequests;
  void *haloRecvRequests;

} dcoo;

typedef enum {PCG=0,GMRES=1}KrylovType;
typedef enum {JACOBI=0,DAMPED_JACOBI=1,CHEBYSHEV=2}SmoothType;

typedef struct agmgLevel_t {
  int Nrows;
  int Ncols;

  int *globalRowStarts; //global partitioning of fine level
  int *globalAggStarts; //global partitioning of coarse level

  void **AxArgs;
  void **smoothArgs;
  void **smootherArgs;
  void **coarsenArgs;
  void **prolongateArgs;

  //operator call-backs
  void (*device_Ax)        (void **args, occa::memory &o_x, occa::memory &o_Ax);
  void (*device_smooth)    (void **args, occa::memory &o_r, occa::memory &o_x, bool x_is_zero);
  void (*device_smoother)  (void **args, occa::memory &o_r, occa::memory &o_Sr);
  void (*device_coarsen)   (void **args, occa::memory &o_x, occa::memory &o_Rx);
  void (*device_prolongate)(void **args, occa::memory &o_x, occa::memory &o_Px);

  //host versions
  void (*Ax)        (void **args, dfloat *x, dfloat *Ax);
  void (*smooth)    (void **args, dfloat *r, dfloat *x, bool x_is_zero);
  void (*smoother)  (void **args, dfloat *r, dfloat *Sr);
  void (*coarsen)   (void **args, dfloat *x, dfloat *Rx);
  void (*prolongate)(void **args, dfloat *x, dfloat *Px);

  //agmg operators
  csr *A;
  csr *P;
  csr *R;

  hyb  *deviceA;
  dcoo  *dcsrP;
  hyb  *deviceR;

  dfloat *rhs, *res, *x;

  dfloat *ckp1, *vkp1, *wkp1;

  occa::memory o_rhs, o_res, o_x;
  occa::memory o_ckp1, o_vkp1, o_wkp1;

  dfloat *smoother_params;
  dfloat *smootherResidual;
  dfloat *smootherResidual2;
  dfloat *smootherUpdate;
  occa::memory o_smootherResidual;
  occa::memory o_smootherResidual2;
  occa::memory o_smootherUpdate;
  int ChebyshevIterations;

  dfloat threshold;
  int numAggregates;
  SmoothType stype;

} agmgLevel;

typedef struct {
  agmgLevel **levels;
  int numLevels;

  KrylovType ktype;

  mesh_t *mesh;
  const char* options;

  //Matrix Free args
  void (*MatFreeAx)(void **args, occa::memory o_q, occa::memory o_Aq,const char* options);
  void **MatFreeArgs;

  //Coarse xxt solver
  void *ExactSolve;
  int coarseTotal;
  int coarseOffset;
  int *coarseOffsets;
  int *coarseCounts;
  dfloat *invCoarseA;
  dfloat *xCoarse, *rhsCoarse;

  bool nullSpace;
  dfloat nullSpacePenalty;

  occa::device device;

  occa::memory o_x;
  occa::memory o_Ax;
  occa::memory o_rho;

  occa::kernel ellAXPYKernel;
  occa::kernel ellZeqAXPYKernel;
  occa::kernel ellJacobiKernel;
  occa::kernel cooAXKernel;
  occa::kernel scaleVectorKernel;
  occa::kernel vectorAddKernel;
  occa::kernel vectorAddKernel2;
  occa::kernel setVectorKernel;
  occa::kernel sumVectorKernel;
  occa::kernel addScalarKernel;
  occa::kernel dotStarKernel;
  occa::kernel simpleDotStarKernel;
  occa::kernel haloExtract;
  occa::kernel agg_interpolateKernel;
  occa::kernel innerProdKernel;
  occa::kernel vectorAddInnerProdKernel;
  occa::kernel kcycleCombinedOp1Kernel;
  occa::kernel kcycleCombinedOp2Kernel;

} parAlmond_t;

parAlmond_t *parAlmondInit(mesh_t *mesh, const char* parAlmondOptions);

void parAlmondAgmgSetup(parAlmond_t* parAlmond,
                       int* rowStarts,
                       int  nnz,
                       int* Ai,
                       int* Aj,
                       dfloat* Avals,
                       bool nullSpace,
                       dfloat nullSpacePenalty);

void parAlmondPrecon(parAlmond_t* parAlmond, occa::memory o_x, occa::memory o_rhs);

int parAlmondFree(void* A);

#endif