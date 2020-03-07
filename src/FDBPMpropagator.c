/********************************************
 * FDBPMpropagator.c, in the C programming language, written for MATLAB MEX function generation
 * Can be compiled with GCC using
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' .\src\FDBPMpropagator.c ".\src\libut.lib" -R2018a"
 * ... or the Microsoft Visual C++ compiler (MSVC) with
 * "copyfile ./src/FDBPMpropagator.c ./src/FDBPMpropagator.cpp; mex COMPFLAGS='/Zp8 /GR /EHs /nologo /MD /openmp /W4 /WX /wd4204 /wd4100' .\src\FDBPMpropagator.cpp ".\src\libut.lib" -R2018a"
 * 
 * The source code in this file is written is such a way that it is
 * compilable by either C or C++ compilers, either with GCC, MSVC or
 * the Nvidia CUDA compiler called NVCC, which is based on MSVC. To
 * compile with CUDA GPU acceleration support, you must have MSVC
 * installed. As of January 2020, mexcuda does not work with MSVC 2019,
 * so I'd recommend MSVC 2017. You also need the Parallel Computing
 * Toolbox, which you will find in the MATLAB addon manager. To compile, run:
 * "copyfile ./src/FDBPMpropagator.c ./src/FDBPMpropagator_CUDA.cu; mexcuda -llibut COMPFLAGS='-use_fast_math -res-usage $COMPFLAGS' .\src\FDBPMpropagator_CUDA.cu -R2018a"
 ********************************************/
// printf("Reached line %d...\n",__LINE__);mexEvalString("drawnow;");mexEvalString("drawnow;");mexEvalString("drawnow;"); // For inserting into code for debugging purposes

#include <math.h>
#include <stdint.h>
#include "mex.h"
#ifdef _OPENMP
  #include "omp.h"
#endif
#ifdef __GNUC__ // This is defined for GCC and CLANG but not for Microsoft Visual C++ compiler
  #define min(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _b: _a;})
  #define max(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _a: _b;})
#endif
#ifndef __clang__
  #ifdef __cplusplus
  extern "C"
  #endif
  extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif
#ifdef __NVCC__
  #include <thrust/complex.h>
  typedef thrust::complex<float> floatcomplex;
  #include <nvml.h>
  #define TILE_DIM 32
#else
  #ifdef __GNUC__
    #include <complex.h>
    typedef float complex floatcomplex;
  #else
    #include <algorithm>
    using namespace std;
    #include <complex>
    typedef std::complex<float> floatcomplex;
  #endif
#endif

struct debug {
  double             dbls[3];
  unsigned long long ulls[3];
};

struct parameters {
  long nx;
  long ny;
  long nz;
  floatcomplex *Efinal;
  floatcomplex *E1;
  floatcomplex *E2;
  floatcomplex *Eyx;
  floatcomplex *b;
  floatcomplex *multiplier;
  floatcomplex ax;
  floatcomplex ay;
};

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep1a(struct parameters *P_global) {
  // Explicit part of substep 1 out of 2
  #ifdef __NVCC__

  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished

  __shared__ double tiledummy[TILE_DIM][TILE_DIM+1]; // We declare with double because a double is the same size as a float complex. +1 is to avoid memory bank conflicts
  floatcomplex *tile = (floatcomplex *)tiledummy;
  
  long xTiles = (P->nx + TILE_DIM - 1)/TILE_DIM;
  long yTiles = (P->ny + TILE_DIM - 1)/TILE_DIM;
  for(long tileNum=blockIdx.x; tileNum<xTiles*yTiles; tileNum += gridDim.x) {
    long tilexoffset = TILE_DIM*(tileNum%xTiles);
    long tileyoffset = TILE_DIM*(tileNum/xTiles);
    long ix = tilexoffset + threadIdx.x;
    long iy = tileyoffset + threadIdx.y;

    if(ix<P->nx && iy<P->ny) {
      long i = ix + iy*P->nx;
      tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] = P->E1[i];
      if(ix != 0      ) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i-1]     - P->E1[i])*P->ax;
      if(ix != P->nx-1) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i+1]     - P->E1[i])*P->ax;
      if(iy != 0      ) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i-P->nx] - P->E1[i])*P->ay*2;
      if(iy != P->ny-1) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i+P->nx] - P->E1[i])*P->ay*2;
    }
    __syncthreads();
    // Save transposed xy -> yx
    ix = tilexoffset + threadIdx.y;
    iy = tileyoffset + threadIdx.x;
    if(ix<P->nx && iy<P->ny) P->Eyx[iy + ix*P->ny] = tile[threadIdx.y + threadIdx.x*(TILE_DIM+1)];
    __syncthreads();
  }
  #else
  long ix,iy;
  struct parameters *P = P_global;
  #ifdef _OPENMP
  #pragma omp for schedule(dynamic)
  #endif
  for(iy=0; iy<P->ny; iy++) {
    for(ix=0; ix<P->nx; ix++) {
      long i = ix + iy*P->nx;

      P->E2[i] = P->E1[i];
      if(ix != 0      ) P->E2[i] += (P->E1[i-1]     - P->E1[i])*P->ax;
      if(ix != P->nx-1) P->E2[i] += (P->E1[i+1]     - P->E1[i])*P->ax;
      if(iy != 0      ) P->E2[i] += (P->E1[i-P->nx] - P->E1[i])*P->ay*2.0f;
      if(iy != P->ny-1) P->E2[i] += (P->E1[i+P->nx] - P->E1[i])*P->ay*2.0f;
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep1b(struct parameters *P_global) {
  // Implicit part of substep 1 out of 2
  #ifdef __NVCC__
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  for(long iy=threadNum;iy<P->ny;iy+=gridDim.x*blockDim.x*blockDim.y){
    for(long ix=0; ix<P->nx; ix++) {
      long i = iy + ix*P->ny;
      P->b[i] = 1;
      if(ix < P->nx-1) P->b[i] += P->ax;
      if(ix > 0) {
        P->b[i]        += P->ax;
        floatcomplex w  = -P->ax/P->b[i-P->ny];
        P->b[i]        += w*P->ax;
        P->Eyx[i]      -= w*P->Eyx[i-P->ny];
      }
    }

    for(long ix=P->nx-1; ix>=0; ix--) {
      long i = iy + ix*P->ny;
      P->Eyx[i] = (P->Eyx[i] + (ix == P->nx-1? 0: P->ax*P->Eyx[i+P->ny]))/P->b[i];
    }
  }
  #else
  struct parameters *P = P_global;
  long i,ix,iy;
  #ifdef _OPENMP
  long threadNum = omp_get_thread_num();
  #pragma omp for schedule(dynamic)
  #else
  long threadNum = 0;
  #endif
  for(iy=0; iy<P->ny; iy++) {
    // Thomson algorithm, sweeps up from 0 to nx-1 and then down from nx-1 to 0:
    // Algorithm is taken from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    for(ix=0; ix<P->nx; ix++) {
      long ib = ix + threadNum*P->nx;
      P->b[ib] = 1;
      if(ix < P->nx-1) P->b[ib] += P->ax;
      if(ix > 0) {
        P->b[ib]        += P->ax;
        floatcomplex w   = -P->ax/P->b[ib-1];
        P->b[ib]        += w*P->ax;
        i                = ix + iy*P->nx;
        P->E2[i]        -= w*P->E2[i-1];
      }
    }

    for(ix=P->nx-1; ix>=0; ix--) {
      long ib = ix + threadNum*P->nx;
      i = ix + iy*P->nx;
      P->E2[i] = (P->E2[i] + (ix == P->nx-1? 0: P->ax*P->E2[i+1]))/P->b[ib];
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep2a(struct parameters *P_global) {
  // Explicit part of substep 2 out of 2
  #ifdef __NVCC__
  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  __shared__ double tiledummy[TILE_DIM][TILE_DIM+1]; // We declare with double because a double is the same size as a float complex. +1 is to avoid memory bank conflicts
  floatcomplex *tile = (floatcomplex *)tiledummy;

  long xTiles = (P->nx + TILE_DIM - 1)/TILE_DIM;
  long yTiles = (P->ny + TILE_DIM - 1)/TILE_DIM;
  for(long tileNum=blockIdx.x; tileNum<xTiles*yTiles; tileNum += gridDim.x) {
    long tilexoffset = TILE_DIM*(tileNum%xTiles);
    long tileyoffset = TILE_DIM*(tileNum/xTiles);
    long ix = tilexoffset + threadIdx.y;
    long iy = tileyoffset + threadIdx.x;

    __syncthreads(); // All threads in the block wait for any previous tile usage to have completed
    if(ix<P->nx && iy<P->ny) tile[threadIdx.y + threadIdx.x*(TILE_DIM+1)] = P->Eyx[ix*P->ny + iy]; // load yx data and store in shared memory in tile, which is xy
    __syncthreads(); // All threads in the block wait for the copy to have finished

    ix = tilexoffset + threadIdx.x;
    iy = tileyoffset + threadIdx.y;
    if(ix<P->nx && iy<P->ny) {
      long i = ix + iy*P->nx;

      floatcomplex deltaE = 0;
      if(iy != 0      ) deltaE -= (P->E1[i-P->nx] - P->E1[i])*P->ay;
      if(iy != P->ny-1) deltaE -= (P->E1[i+P->nx] - P->E1[i])*P->ay;
      P->E2[i] = tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] + deltaE;
    }
  }
  #else
  struct parameters *P = P_global;
  long i,ix,iy;
  #ifdef _OPENMP
  #pragma omp for schedule(dynamic)
  #endif
  for(iy=0; iy<P->ny; iy++) {
    for(ix=0; ix<P->nx; ix++) {
      i = ix + iy*P->nx;

      if(iy != 0      ) P->E2[i] -= (P->E1[i-P->nx] - P->E1[i])*P->ay;
      if(iy != P->ny-1) P->E2[i] -= (P->E1[i+P->nx] - P->E1[i])*P->ay;
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep2b(struct parameters *P_global) {
  // Implicit part of substep 2 out of 2
  #ifdef __NVCC__
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  for(long ix=threadNum;ix<P->nx;ix+=gridDim.x*blockDim.x*blockDim.y) {
    for(long iy=0; iy<P->ny; iy++) {
      long i = ix + iy*P->nx;
      P->b[i] = 1;
      if(iy < P->ny-1) P->b[i] += P->ay;
      if(iy > 0) {
        P->b[i]         += P->ay;
        floatcomplex w  = -P->ay/P->b[i-P->nx];
        P->b[i]         += w*P->ay;
        P->E2[i]        -= w*P->E2[i-P->nx];
      }
    }

    for(long iy=P->ny-1; iy>=0; iy--) {
      long i = ix + iy*P->nx;
      P->E2[i] = (P->E2[i] + (iy == P->ny-1? 0: P->ay*P->E2[i+P->nx]))/P->b[i];
      if(iy < P->ny-1) P->E2[i+P->nx] *= P->multiplier[i];
    }
    P->E2[ix] *= P->multiplier[ix];
  }

  #else
  struct parameters *P = P_global;
  long i,ix,iy;
  #ifdef _OPENMP
  long threadNum = omp_get_thread_num();
  #pragma omp for schedule(dynamic)
  #else
  long threadNum = 0;
  #endif
  for(ix=0; ix<P->nx; ix++) {
    for(iy=0; iy<P->ny; iy++) {
      long ib = iy + threadNum*P->ny;
      P->b[ib] = 1;
      if(iy < P->ny-1) P->b[ib] += P->ay;
      if(iy > 0) {
        P->b[ib]        += P->ay;
        floatcomplex w   = -P->ay/P->b[ib-1];
        P->b[ib]        += w*P->ay;
        i                = ix + iy*P->nx;
        P->E2[i]        -= w*P->E2[i-P->nx];
      }
    }

    for(iy=P->ny-1; iy>=0; iy--) {
      long ib = iy + threadNum*P->ny;
      i = ix + iy*P->nx;
      P->E2[i] = (P->E2[i] + (iy == P->ny-1? 0: P->ay*P->E2[i+P->nx]))/P->b[ib];
    }
  }

  #ifdef _OPENMP
  #pragma omp for schedule(dynamic)
  #endif
  for(iy=0; iy<P->ny; iy++) {
    for(ix=0; ix<P->nx; ix++) {
      i = ix + iy*P->nx;
      P->E2[i] *= P->multiplier[i];
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void swapEPointers(struct parameters *P, long n) {
  #ifdef __NVCC__
  floatcomplex *temp = P->E1;
  P->E1 = P->E2;
  P->E2 = temp;
  #else
  if(n>0) { // Swap E1 and E2
    floatcomplex *temp = P->E1;
    P->E1 = P->E2;
    P->E2 = temp;
  } else if(P->nz%2) {
    P->E1 = P->E2;
    P->E2 = (floatcomplex *)malloc(P->nx*P->ny*sizeof(floatcomplex));
  } else {
    P->E1 = P->E2;
    P->E2 = P->Efinal;
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line) {
  if (code != cudaSuccess) {
    printf("GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    mexEvalString("drawnow;");
    while(true) {;}
  }
}

void createDeviceStructs(struct parameters *P, struct parameters **P_devptr,
                         struct debug *D, struct debug **D_devptr) {
  long N = P->nx*P->ny;
  struct parameters P_tempvar = *P;

  gpuErrchk(cudaMalloc(&P_tempvar.E1,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMemcpy(P_tempvar.E1,P->E1,N*sizeof(floatcomplex),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.E2,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMalloc(&P_tempvar.Eyx,N*sizeof(floatcomplex)));

  gpuErrchk(cudaMalloc(&P_tempvar.b,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMalloc(&P_tempvar.multiplier,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMemcpy(P_tempvar.multiplier,P->multiplier,N*sizeof(floatcomplex),cudaMemcpyHostToDevice));

  gpuErrchk(cudaMalloc(P_devptr, sizeof(struct parameters)));
  gpuErrchk(cudaMemcpy(*P_devptr,&P_tempvar,sizeof(struct parameters),cudaMemcpyHostToDevice));

  // Allocate and copy debug struct
  struct debug D_tempvar = *D;
  gpuErrchk(cudaMalloc(D_devptr, sizeof(struct debug)));
  gpuErrchk(cudaMemcpy(*D_devptr,&D_tempvar,sizeof(struct debug),cudaMemcpyHostToDevice));
}

void retrieveAndFreeDeviceStructs(struct parameters *P, struct parameters *P_dev,
                                  struct debug *D, struct debug *D_dev) {
  long N = P->nx*P->ny;
  struct parameters P_temp; gpuErrchk(cudaMemcpy(&P_temp,P_dev,sizeof(struct parameters),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(P->Efinal,P_temp.E2,N*sizeof(floatcomplex),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(P_temp.E1));
  gpuErrchk(cudaFree(P_temp.E2));
  gpuErrchk(cudaFree(P_temp.Eyx));
  gpuErrchk(cudaFree(P_temp.b));
  gpuErrchk(cudaFree(P_temp.multiplier));
  gpuErrchk(cudaFree(P_dev));

  gpuErrchk(cudaMemcpy(D, D_dev, sizeof(struct debug),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(D_dev));
}
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
  struct parameters P_var;
  struct parameters *P = &P_var;
  P->nx = (long)mxGetM(prhs[0]);
  P->ny = (long)mxGetN(prhs[0]);
  P->nz = (long)*(double *)mxGetData(mxGetField(prhs[1],0,"nz")); // Number of steps to perform
  P->E1 = (floatcomplex *)mxGetData(prhs[0]); // Input E field
  mwSize const *dimPtr = mxGetDimensions(prhs[0]);
  P->Efinal = (floatcomplex *)mxGetData(plhs[0] = mxCreateNumericArray(2,dimPtr,mxSINGLE_CLASS,mxCOMPLEX)); // Output E field
  #ifndef __NVCC__
  P->E2 = (floatcomplex *)(P->nz%2? P->Efinal: malloc(P->nx*P->ny*sizeof(floatcomplex)));
  #endif
  P->multiplier = (floatcomplex *)mxGetData(mxGetField(prhs[1],0,"multiplier")); // Array of multiplier values to apply to the E field after each step, which includes (1) the phase imparted by the refractive index differences, (2) the absorber outside the fiber and (3) the effects of fibre bending, if present
  P->ax = *(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ax"));
  P->ay = *(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ay"));
  
  bool ctrlc_caught = false;      // Has a ctrl+c been passed from MATLAB?
  #ifdef __NVCC__
  int temp, nBlocks; gpuErrchk(cudaOccupancyMaxPotentialBlockSize(&nBlocks,&temp,&substep1a,0,0));
  dim3 blockDims(TILE_DIM,TILE_DIM,1);

  struct parameters *P_dev;
  struct debug D_var = {{0.0,0.0,0.0},{0,0,0}};
  struct debug *D = &D_var;
  struct debug *D_dev;
  createDeviceStructs(P,&P_dev,D,&D_dev);
  #else
  #ifdef _OPENMP
  bool useAllCPUs = mxIsLogicalScalarTrue(mxGetField(prhs[1],0,"useAllCPUs"));
  long numThreads = useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1;
  #else
  long numThreads = 1;
  #endif
  P->b = (floatcomplex *)malloc(numThreads*max(P->nx,P->ny)*sizeof(floatcomplex));
  #ifdef _OPENMP
  #pragma omp parallel num_threads(useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1)
  #endif
  #endif
  {
    for(long iz=0; iz<P->nz; iz++) {
      if(ctrlc_caught) break;
      
      #ifdef __NVCC__
      substep1a<<<nBlocks, blockDims>>>(P_dev); // xy -> yx
      substep1b<<<nBlocks, blockDims>>>(P_dev); // yx -> yx
      substep2a<<<nBlocks, blockDims>>>(P_dev); // yx -> xy
      substep2b<<<nBlocks, blockDims>>>(P_dev); // xy -> xy
      gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
      #else
      substep1a(P);
      substep1b(P);
      substep2a(P);
      substep2b(P);
      #endif

      #ifdef _OPENMP
      #pragma omp master
      #endif
      {
        #ifndef __clang__
        if(utIsInterruptPending()) {
          ctrlc_caught = true;
          printf("\nCtrl+C detected, stopping.\n");
        }
        #endif

        if(iz+1 < P->nz) {
          #ifdef __NVCC__
          swapEPointers<<<1,1>>>(P_dev,iz);
          #else
          swapEPointers(P,iz);
          #endif
        }
      }
      #ifdef _OPENMP
      #pragma omp barrier
      #endif
    }
  }
  #ifdef __NVCC__
  gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
  retrieveAndFreeDeviceStructs(P,P_dev,D,D_dev);
//   printf("\nDebug: %.18e %.18e %.18e %llu %llu %llu\n          ",D->dbls[0],D->dbls[1],D->dbls[2],D->ulls[0],D->ulls[1],D->ulls[2]);
  #else
  if(P->nz > 1) free(P->E1);
  free(P->b);
  #endif

  return;
}
