/********************************************
 * FDBPMpropagatorSingles.c, in the C programming language, written for MATLAB MEX function generation
 * Can be compiled with GCC using
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' .\src\FDBPMpropagatorSingles.c ".\src\libut.lib" -R2018a"
 * ... or the Microsoft Visual C++ compiler (MSVC) with
 * "copyfile ./src/FDBPMpropagatorSingles.c ./src/FDBPMpropagatorSingles.cpp; mex COMPFLAGS='/Zp8 /GR /EHs /nologo /MD /openmp /W4 /WX /wd4204 /wd4100' .\src\FDBPMpropagatorSingles.cpp ".\src\libut.lib" -R2018a"
 * 
 * The source code in this file is written is such a way that it is
 * compilable by either C or C++ compilers, either with GCC, MSVC or
 * the Nvidia CUDA compiler called NVCC, which is based on MSVC. To
 * compile with CUDA GPU acceleration support, you must have MSVC
 * installed. As of January 2020, mexcuda does not work with MSVC 2019,
 * so I'd recommend MSVC 2017. You also need the Parallel Computing
 * Toolbox, which you will find in the MATLAB addon manager. To compile, run:
 * "copyfile ./src/FDBPMpropagatorSingles.c ./src/FDBPMpropagatorSingles_CUDA.cu; mexcuda -llibut COMPFLAGS='-use_fast_math -res-usage $COMPFLAGS' .\src\FDBPMpropagatorSingles_CUDA.cu -R2018a"
 ********************************************/
// printf("Reached line %d...\n",__LINE__);mexEvalString("drawnow;");mexEvalString("drawnow;");mexEvalString("drawnow;"); // For inserting into code for debugging purposes

#include <math.h>
#include <stdint.h>
#include "mex.h"
#define PI acosf(-1.0f)
#ifdef _OPENMP
  #include "omp.h"
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
  #define I thrust::complex<float>{0,1}
  #define CEXPF(x) (thrust::exp(x))
  #define MAX(x,y) (max(x,y))
  #include <nvml.h>
  #define TILE_DIM 32
#else
  #ifdef __GNUC__ // This is defined for GCC and CLANG but not for Microsoft Visual C++ compiler
    #define MAX(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _a: _b;})
    #include <complex.h>
    typedef float complex floatcomplex;
    #define CEXPF(x) (cexpf(x))
  #else
    #include <algorithm>
    #include <complex>
    typedef std::complex<float> floatcomplex;
    #define I std::complex<float>{0,1}
    #define CEXPF(x) (std::exp(x))
    #define MAX(x,y) (std::max(x,y))
  #endif
#endif

struct debug {
  double             dbls[3];
  unsigned long long ulls[3];
};

struct parameters {
  long Nx;
  long Ny;
  float dx;
  float dy;
  long iz_start;
  long iz_end;
  float taperPerStep;
  float twistPerStep;
  float d;
  float n_cladding;
  float n_0;
  long Nshapes;
  float *shapexs;
  float *shapeys;
  float *shapeRs;
  float *shapeTypes;
  float *shapeRIs;
  float *shapegs;
  float *shapexs_transformed;
  float *shapeys_transformed;
  float *shapeRs_transformed;
  floatcomplex *Efinal;
  floatcomplex *E1;
  floatcomplex *E2;
  floatcomplex *Eyx;
  float *n_out;
  floatcomplex *b;
  floatcomplex *multiplier;
  floatcomplex ax;
  floatcomplex ay;
  float rho_e;
  float RoC;
  float sinBendDirection; 
  float cosBendDirection; 
};

#ifdef __NVCC__ // If compiling for CUDA
__host__ __device__
#endif
float sqrf(float x) {return x*x;}

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
  
  long xTiles = (P->Nx + TILE_DIM - 1)/TILE_DIM;
  long yTiles = (P->Ny + TILE_DIM - 1)/TILE_DIM;
  for(long tileNum=blockIdx.x; tileNum<xTiles*yTiles; tileNum += gridDim.x) {
    long tilexoffset = TILE_DIM*(tileNum%xTiles);
    long tileyoffset = TILE_DIM*(tileNum/xTiles);
    long ix = tilexoffset + threadIdx.x;
    long iy = tileyoffset + threadIdx.y;

    if(ix<P->Nx && iy<P->Ny) {
      long i = ix + iy*P->Nx;
      tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] = P->E1[i];
      if(ix != 0      ) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i-1]     - P->E1[i])*P->ax;
      if(ix != P->Nx-1) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i+1]     - P->E1[i])*P->ax;
      if(iy != 0      ) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i-P->Nx] - P->E1[i])*P->ay*2;
      if(iy != P->Ny-1) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i+P->Nx] - P->E1[i])*P->ay*2;
    }
    __syncthreads();
    // Save transposed xy -> yx
    ix = tilexoffset + threadIdx.y;
    iy = tileyoffset + threadIdx.x;
    if(ix<P->Nx && iy<P->Ny) P->Eyx[iy + ix*P->Ny] = tile[threadIdx.y + threadIdx.x*(TILE_DIM+1)];
    __syncthreads();
  }
  #else
  long ix,iy;
  struct parameters *P = P_global;
  #ifdef _OPENMP
  #pragma omp for schedule(dynamic)
  #endif
  for(iy=0; iy<P->Ny; iy++) {
    for(ix=0; ix<P->Nx; ix++) {
      long i = ix + iy*P->Nx;

      P->E2[i] = P->E1[i];
      if(ix != 0      ) P->E2[i] += (P->E1[i-1]     - P->E1[i])*P->ax;
      if(ix != P->Nx-1) P->E2[i] += (P->E1[i+1]     - P->E1[i])*P->ax;
      if(iy != 0      ) P->E2[i] += (P->E1[i-P->Nx] - P->E1[i])*P->ay*2.0f;
      if(iy != P->Ny-1) P->E2[i] += (P->E1[i+P->Nx] - P->E1[i])*P->ay*2.0f;
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
  for(long iy=threadNum;iy<P->Ny;iy+=gridDim.x*blockDim.x*blockDim.y){
    for(long ix=0; ix<P->Nx; ix++) {
      long i = iy + ix*P->Ny;
      P->b[i] = 1;
      if(ix < P->Nx-1) P->b[i] += P->ax;
      if(ix > 0) {
        P->b[i]        += P->ax;
        floatcomplex w  = -P->ax/P->b[i-P->Ny];
        P->b[i]        += w*P->ax;
        P->Eyx[i]      -= w*P->Eyx[i-P->Ny];
      }
    }

    for(long ix=P->Nx-1; ix>=0; ix--) {
      long i = iy + ix*P->Ny;
      P->Eyx[i] = (P->Eyx[i] + (ix == P->Nx-1? 0: P->ax*P->Eyx[i+P->Ny]))/P->b[i];
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
  for(iy=0; iy<P->Ny; iy++) {
    // Thomson algorithm, sweeps up from 0 to Nx-1 and then down from Nx-1 to 0:
    // Algorithm is taken from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    for(ix=0; ix<P->Nx; ix++) {
      long ib = ix + threadNum*P->Nx;
      P->b[ib] = 1;
      if(ix < P->Nx-1) P->b[ib] += P->ax;
      if(ix > 0) {
        P->b[ib]        += P->ax;
        floatcomplex w   = -P->ax/P->b[ib-1];
        P->b[ib]        += w*P->ax;
        i                = ix + iy*P->Nx;
        P->E2[i]        -= w*P->E2[i-1];
      }
    }

    for(ix=P->Nx-1; ix>=0; ix--) {
      long ib = ix + threadNum*P->Nx;
      i = ix + iy*P->Nx;
      P->E2[i] = (P->E2[i] + (ix == P->Nx-1? 0: P->ax*P->E2[i+1]))/P->b[ib];
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

  long xTiles = (P->Nx + TILE_DIM - 1)/TILE_DIM;
  long yTiles = (P->Ny + TILE_DIM - 1)/TILE_DIM;
  for(long tileNum=blockIdx.x; tileNum<xTiles*yTiles; tileNum += gridDim.x) {
    long tilexoffset = TILE_DIM*(tileNum%xTiles);
    long tileyoffset = TILE_DIM*(tileNum/xTiles);
    long ix = tilexoffset + threadIdx.y;
    long iy = tileyoffset + threadIdx.x;

    __syncthreads(); // All threads in the block wait for any previous tile usage to have completed
    if(ix<P->Nx && iy<P->Ny) tile[threadIdx.y + threadIdx.x*(TILE_DIM+1)] = P->Eyx[ix*P->Ny + iy]; // load yx data and store in shared memory in tile, which is xy
    __syncthreads(); // All threads in the block wait for the copy to have finished

    ix = tilexoffset + threadIdx.x;
    iy = tileyoffset + threadIdx.y;
    if(ix<P->Nx && iy<P->Ny) {
      long i = ix + iy*P->Nx;

      floatcomplex deltaE = 0;
      if(iy != 0      ) deltaE -= (P->E1[i-P->Nx] - P->E1[i])*P->ay;
      if(iy != P->Ny-1) deltaE -= (P->E1[i+P->Nx] - P->E1[i])*P->ay;
      P->E2[i] = tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] + deltaE;
    }
  }
  #else
  struct parameters *P = P_global;
  long i,ix,iy;
  #ifdef _OPENMP
  #pragma omp for schedule(dynamic)
  #endif
  for(iy=0; iy<P->Ny; iy++) {
    for(ix=0; ix<P->Nx; ix++) {
      i = ix + iy*P->Nx;

      if(iy != 0      ) P->E2[i] -= (P->E1[i-P->Nx] - P->E1[i])*P->ay;
      if(iy != P->Ny-1) P->E2[i] -= (P->E1[i+P->Nx] - P->E1[i])*P->ay;
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
  for(long ix=threadNum;ix<P->Nx;ix+=gridDim.x*blockDim.x*blockDim.y) {
    for(long iy=0; iy<P->Ny; iy++) {
      long i = ix + iy*P->Nx;
      P->b[i] = 1;
      if(iy < P->Ny-1) P->b[i] += P->ay;
      if(iy > 0) {
        P->b[i]         += P->ay;
        floatcomplex w  = -P->ay/P->b[i-P->Nx];
        P->b[i]         += w*P->ay;
        P->E2[i]        -= w*P->E2[i-P->Nx];
      }
    }

    for(long iy=P->Ny-1; iy>=0; iy--) {
      long i = ix + iy*P->Nx;
      P->E2[i] = (P->E2[i] + (iy == P->Ny-1? 0: P->ay*P->E2[i+P->Nx]))/P->b[i];
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
  for(ix=0; ix<P->Nx; ix++) {
    for(iy=0; iy<P->Ny; iy++) {
      long ib = iy + threadNum*P->Ny;
      P->b[ib] = 1;
      if(iy < P->Ny-1) P->b[ib] += P->ay;
      if(iy > 0) {
        P->b[ib]        += P->ay;
        floatcomplex w   = -P->ay/P->b[ib-1];
        P->b[ib]        += w*P->ay;
        i                = ix + iy*P->Nx;
        P->E2[i]        -= w*P->E2[i-P->Nx];
      }
    }

    for(iy=P->Ny-1; iy>=0; iy--) {
      long ib = iy + threadNum*P->Ny;
      i = ix + iy*P->Nx;
      P->E2[i] = (P->E2[i] + (iy == P->Ny-1? 0: P->ay*P->E2[i+P->Nx]))/P->b[ib];
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void calcShapexyr(struct parameters *P, long iz) { // Called with only one thread
  // Calculate positions and sizes of geometric shapes
  float cosvalue = cosf(P->twistPerStep*iz);
  float sinvalue = sinf(P->twistPerStep*iz);
  float scaling = 1 - P->taperPerStep*iz;
  for(long iShape=0;iShape<P->Nshapes;iShape++) {
    P->shapexs_transformed[iShape] = scaling*(cosvalue*P->shapexs[iShape] - sinvalue*P->shapeys[iShape]);
    P->shapeys_transformed[iShape] = scaling*(sinvalue*P->shapexs[iShape] + cosvalue*P->shapeys[iShape]);
    P->shapeRs_transformed[iShape] = scaling*P->shapeRs[iShape];
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void applyMultiplier(struct parameters *P_global, long iz) {
  #ifdef __NVCC__
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying
  __syncthreads(); // All threads in the block wait for the copy to have finished

  if(P->taperPerStep || P->twistPerStep) {
    for(long i=threadNum;i<P->Nx*P->Ny;i+=gridDim.x*blockDim.x*blockDim.y) {
  #else
  struct parameters *P = P_global;
  if(P->taperPerStep || P->twistPerStep) {
    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for(long i=0;i<P->Nx*P->Ny;i++) {
  #endif
      // For each pixel, calculate refractive index and multiply
      long ix = i%P->Nx;
      float x = P->dx*(ix - (P->Nx-1)/2.0f);
      long iy = i/P->Nx;
      float y = P->dy*(iy - (P->Ny-1)/2.0f);

      float n = P->n_cladding;
      for(long iShape=0;iShape<P->Nshapes;iShape++) {
        switch((int)P->shapeTypes[iShape]) {
          case 1: // Step-index disk
            if(sqrf(x - P->shapexs_transformed[iShape]) + sqrf(y - P->shapeys_transformed[iShape]) < sqrf(P->shapeRs_transformed[iShape]))
              n = P->shapeRIs[iShape];
            break;
          case 2: { // Antialiased step-index disk
            float delta = MAX(P->dx,P->dy); // Width of antialiasing slope
            float r_diff = sqrtf(sqrf(x - P->shapexs_transformed[iShape]) + sqrf(y - P->shapeys_transformed[iShape])) - P->shapeRs_transformed[iShape] + delta/2.0f;
            if(r_diff < 0) {
              n = P->shapeRIs[iShape];
            } else if(r_diff < delta) {
              n = r_diff/delta*(P->n_cladding - P->shapeRIs[iShape]) + P->shapeRIs[iShape];
            }
            break;
          }
          case 3: { // Parabolic graded index disk
            float r_ratio_sqr = (sqrf(x - P->shapexs_transformed[iShape]) + sqrf(y - P->shapeys_transformed[iShape]))/sqrf(P->shapeRs_transformed[iShape]);
            if(r_ratio_sqr < 1)
              n = r_ratio_sqr*(P->n_cladding - P->shapeRIs[iShape]) + P->shapeRIs[iShape];
            break;
          }
          case 4: { // 2D Hyperbolic GRIN lens
              float r_ratio_sqr = (sqrf(x - P->shapexs_transformed[iShape])+sqrf(y - P->shapeys_transformed[iShape]))/sqrf(P->shapeRs_transformed[iShape]);
              float r_abs = sqrtf(sqrf(x - P->shapexs_transformed[iShape])+sqrf(y - P->shapeys_transformed[iShape]));
              if(r_ratio_sqr < 1)
                n =  2*P->shapeRIs[iShape] * exp(P->shapegs[iShape]*r_abs) / (exp(2*P->shapegs[iShape]*r_abs)+1); // GRINTECH: n = n_0 * sech(gr)
            break;
          }
          case 5: { // 1D (y) Hyperbolic GRIN lens
              float r_ratio_sqr = sqrf(y - P->shapeys_transformed[iShape])/sqrf(P->shapeRs_transformed[iShape]);
              float r_abs = y - P->shapeys_transformed[iShape];
              if(r_ratio_sqr < 1)
                n =  2*P->shapeRIs[iShape] * exp(P->shapegs[iShape]*r_abs) / (exp(2*P->shapegs[iShape]*r_abs)+1); 
            break;
          }
        }
      }
      float n_eff = n*(1-(sqrf(n)*(x*P->cosBendDirection+y*P->sinBendDirection)/2/P->RoC*P->rho_e))*exp((x*P->cosBendDirection+y*P->sinBendDirection)/P->RoC);
      if(iz == P->iz_end-1) P->n_out[i] = n_eff;
      P->E2[i] *= P->multiplier[i]*CEXPF(I*P->d*(sqrf(n_eff) - sqrf(P->n_0)));
    }
  } else {
  #ifdef __NVCC__
    for(long i=threadNum;i<P->Nx*P->Ny;i+=gridDim.x*blockDim.x*blockDim.y) P->E2[i] *= P->multiplier[i];
  #else
    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for(long i=0;i<P->Nx*P->Ny;i++) {
      P->E2[i] *= P->multiplier[i];
    }
  #endif
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void swapEPointers(struct parameters *P, long iz) {
  #ifdef __NVCC__
  floatcomplex *temp = P->E1;
  P->E1 = P->E2;
  P->E2 = temp;
  #else
  if(iz>P->iz_start) { // Swap E1 and E2
    floatcomplex *temp = P->E1;
    P->E1 = P->E2;
    P->E2 = temp;
  } else if((P->iz_end - P->iz_start)%2) {
    P->E1 = P->E2;
    P->E2 = (floatcomplex *)malloc(P->Nx*P->Ny*sizeof(floatcomplex));
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
  long N = P->Nx*P->Ny;
  struct parameters P_tempvar = *P;

  gpuErrchk(cudaMalloc(&P_tempvar.shapexs,P->Nshapes*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.shapexs,P->shapexs,P->Nshapes*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.shapeys,P->Nshapes*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.shapeys,P->shapeys,P->Nshapes*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.shapeRs,P->Nshapes*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.shapeRs,P->shapeRs,P->Nshapes*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.shapeTypes,P->Nshapes*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.shapeTypes,P->shapeTypes,P->Nshapes*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.shapeRIs,P->Nshapes*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.shapeRIs,P->shapeRIs,P->Nshapes*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.shapegs,P->Nshapes*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.shapegs,P->shapegs,P->Nshapes*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.shapexs_transformed,P->Nshapes*sizeof(float)));
  gpuErrchk(cudaMalloc(&P_tempvar.shapeys_transformed,P->Nshapes*sizeof(float)));
  gpuErrchk(cudaMalloc(&P_tempvar.shapeRs_transformed,P->Nshapes*sizeof(float)));

  gpuErrchk(cudaMalloc(&P_tempvar.E1,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMemcpy(P_tempvar.E1,P->E1,N*sizeof(floatcomplex),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.E2,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMalloc(&P_tempvar.Eyx,N*sizeof(floatcomplex)));

  if(P->taperPerStep || P->twistPerStep) gpuErrchk(cudaMalloc(&P_tempvar.n_out,N*sizeof(float)));

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
  long N = P->Nx*P->Ny;
  struct parameters P_temp; gpuErrchk(cudaMemcpy(&P_temp,P_dev,sizeof(struct parameters),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(P->Efinal,P_temp.E2,N*sizeof(floatcomplex),cudaMemcpyDeviceToHost));
  if(P->taperPerStep || P->twistPerStep) {
    gpuErrchk(cudaMemcpy(P->n_out,P_temp.n_out,N*sizeof(float),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(P_temp.n_out));
  }

  gpuErrchk(cudaFree(P_temp.shapexs));
  gpuErrchk(cudaFree(P_temp.shapeys));
  gpuErrchk(cudaFree(P_temp.shapeRs));
  gpuErrchk(cudaFree(P_temp.shapeTypes));
  gpuErrchk(cudaFree(P_temp.shapeRIs));
  gpuErrchk(cudaFree(P_temp.shapegs));
  gpuErrchk(cudaFree(P_temp.shapexs_transformed));
  gpuErrchk(cudaFree(P_temp.shapeys_transformed));
  gpuErrchk(cudaFree(P_temp.shapeRs_transformed));

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
  P->Nx = (long)mxGetM(prhs[0]);
  P->Ny = (long)mxGetN(prhs[0]);
  P->dx = *(float *)mxGetData(mxGetField(prhs[1],0,"dx"));
  P->dy = *(float *)mxGetData(mxGetField(prhs[1],0,"dy"));
  P->iz_start = *(long *)mxGetData(mxGetField(prhs[1],0,"iz_start"));
  P->iz_end = *(long *)mxGetData(mxGetField(prhs[1],0,"iz_end"));
  P->taperPerStep = *(float *)mxGetData(mxGetField(prhs[1],0,"taperPerStep"));
  P->twistPerStep = *(float *)mxGetData(mxGetField(prhs[1],0,"twistPerStep"));
  P->d = *(float *)mxGetData(mxGetField(prhs[1],0,"d"));
  P->n_cladding = *(float *)mxGetData(mxGetField(prhs[1],0,"n_cladding"));
  P->n_0 = *(float *)mxGetData(mxGetField(prhs[1],0,"n_0"));
  P->Nshapes = (long)mxGetM(mxGetField(prhs[1],0,"shapes"));
  P->shapexs = (float *)mxGetData(mxGetField(prhs[1],0,"shapes"));
  P->shapeys = P->shapexs + P->Nshapes;
  P->shapeRs = P->shapeys + P->Nshapes;
  P->shapeTypes = P->shapeRs + P->Nshapes;
  P->shapeRIs = P->shapeTypes + P->Nshapes;
  P->shapegs = P->shapeRIs + P->Nshapes;
  P->rho_e = *(float *)mxGetData(mxGetField(prhs[1],0,"rho_e"));
  P->RoC = *(float *)mxGetData(mxGetField(prhs[1],0,"RoC"));
  P->sinBendDirection = sin(*(float *)mxGetData(mxGetField(prhs[1],0,"bendDirection"))/180*PI);
  P->cosBendDirection = cos(*(float *)mxGetData(mxGetField(prhs[1],0,"bendDirection"))/180*PI);
  P->shapexs_transformed = (float *)malloc(P->Nshapes*sizeof(float));
  P->shapeys_transformed = (float *)malloc(P->Nshapes*sizeof(float));
  P->shapeRs_transformed = (float *)malloc(P->Nshapes*sizeof(float));
  P->E1 = (floatcomplex *)mxGetData(prhs[0]); // Input E field
  mwSize const *dimPtr = mxGetDimensions(prhs[0]);
  P->Efinal = (floatcomplex *)mxGetData(plhs[0] = mxCreateNumericArray(2,dimPtr,mxSINGLE_CLASS,mxCOMPLEX)); // Output E field
  P->n_out= (float *)mxGetData(plhs[1] = mxCreateNumericArray(2,dimPtr,mxSINGLE_CLASS,mxREAL)); // Output refractive index, calculated based on the geometric shapes definitions
  #ifndef __NVCC__
  P->E2 = (floatcomplex *)((P->iz_end - P->iz_start)%2? P->Efinal: malloc(P->Nx*P->Ny*sizeof(floatcomplex)));
  #endif
  floatcomplex *MatlabMultiplier = (floatcomplex *)mxGetData(mxGetField(prhs[1],0,"multiplier")); // Array of multiplier values to apply to the E field after each step, which includes (1) the absorber outside the fiber and (2) the effects of fibre bending, if present
  P->multiplier = (floatcomplex *)malloc(P->Nx*P->Ny*sizeof(floatcomplex)); // If there's tapering or twisting, this is just a copy of MatlabMultiplier. Otherwise, it will be MatlabMultiplier multiplied by the factor that takes the refractive index profile into account, cexpf(I*P->d*(sqrf(n) - sqrf(P->n_0))).
  P->ax = *(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ax"));
  P->ay = *(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ay"));
  
  if(P->taperPerStep || P->twistPerStep) {
    for(long i=0;i<P->Nx*P->Ny;i++) P->multiplier[i] = MatlabMultiplier[i];
  } else {
    for(long ix=0;ix<P->Nx;ix++) {
      float x = P->dx*(ix - (P->Nx-1)/2.0f);
      for(long iy=0;iy<P->Ny;iy++) {
        float y = P->dy*(iy - (P->Ny-1)/2.0f);
        long i = ix + iy*P->Nx;
        float n = P->n_cladding;
        for(long iShape=0;iShape<P->Nshapes;iShape++) {
          switch((int)P->shapeTypes[iShape]) {
            case 1: // Step-index disk
              if(sqrf(x - P->shapexs[iShape]) + sqrf(y - P->shapeys[iShape]) < sqrf(P->shapeRs[iShape]))
                n = P->shapeRIs[iShape];
              break;
            case 2: { // Antialiased step-index disk
              float delta = MAX(P->dx,P->dy); // Width of antialiasing slope
              float r_diff = sqrtf(sqrf(x - P->shapexs[iShape]) + sqrf(y - P->shapeys[iShape])) - P->shapeRs[iShape] + delta/2.0f;
              if(r_diff < 0) {
                n = P->shapeRIs[iShape];
              } else if(r_diff < delta) {
                n = r_diff/delta*(P->n_cladding - P->shapeRIs[iShape]) + P->shapeRIs[iShape];
              }
              break;
            }
            case 3: { // Parabolic graded index disk
              float r_ratio_sqr = (sqrf(x - P->shapexs[iShape]) + sqrf(y - P->shapeys[iShape]))/sqrf(P->shapeRs[iShape]);
              if(r_ratio_sqr < 1)
                n = r_ratio_sqr*(P->n_cladding - P->shapeRIs[iShape]) + P->shapeRIs[iShape];
              break;
            }
            case 4: { // 2D Hyperbolic GRIN lens
              float r_ratio_sqr = (sqrf(x - P->shapexs[iShape])+sqrf(y - P->shapeys[iShape]))/sqrf(P->shapeRs[iShape]);
              float r_abs = sqrtf(sqrf(x - P->shapexs[iShape])+sqrf(y - P->shapeys[iShape]));
              if(r_ratio_sqr < 1)
                n =  2*P->shapeRIs[iShape] * exp(P->shapegs[iShape]*r_abs) / (exp(2*P->shapegs[iShape]*r_abs)+1);  // GRINTECH: n = n_0 * sech(gr) = n_0*2*exp(gr)/(exp(2gr)+1)
              break;
            }
            case 5: { // 1D (y) Hyperbolic GRIN lens
              float r_ratio_sqr = sqrf(y - P->shapeys[iShape])/sqrf(P->shapeRs[iShape]);
              float r_abs = y - P->shapeys[iShape];
              if(r_ratio_sqr < 1)
                n =  2*P->shapeRIs[iShape] * exp(P->shapegs[iShape]*r_abs) / (exp(2*P->shapegs[iShape]*r_abs)+1); 
              break;
            }
          }
        }
        float n_eff = n*(1-(sqrf(n)*(x*P->cosBendDirection+y*P->sinBendDirection)/2/P->RoC*P->rho_e))*exp((x*P->cosBendDirection+y*P->sinBendDirection)/P->RoC);
        P->n_out[i] = n_eff;
        P->multiplier[i] = MatlabMultiplier[i]*CEXPF(I*P->d*(sqrf(n_eff) - sqrf(P->n_0)));
      }
    }
  }
  
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
  P->b = (floatcomplex *)malloc(numThreads*MAX(P->Nx,P->Ny)*sizeof(floatcomplex));
  #ifdef _OPENMP
  #pragma omp parallel num_threads(useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1)
  #endif
  #endif
  {
    for(long iz=P->iz_start; iz<P->iz_end; iz++) {
      if(ctrlc_caught) break;
      
      #ifdef __NVCC__
      substep1a<<<nBlocks, blockDims>>>(P_dev); // xy -> yx
      substep1b<<<nBlocks, blockDims>>>(P_dev); // yx -> yx
      substep2a<<<nBlocks, blockDims>>>(P_dev); // yx -> xy
      substep2b<<<nBlocks, blockDims>>>(P_dev); // xy -> xy
      if(P->taperPerStep || P->twistPerStep) calcShapexyr<<<1,1>>>(P_dev,iz);
      applyMultiplier<<<nBlocks, blockDims>>>(P_dev,iz); // xy -> xy
      gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
      #else
      substep1a(P);
      substep1b(P);
      substep2a(P);
      substep2b(P);
      if(P->taperPerStep || P->twistPerStep) calcShapexyr(P,iz);
      applyMultiplier(P,iz);
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

        if(iz+1 < P->iz_end) {
          #ifdef __NVCC__
          swapEPointers<<<1,1>>>(P_dev,iz);
          #else
          swapEPointers(P,iz);
          #endif
        }
        #endif
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
  if(P->E1 != mxGetData(prhs[0]) && P->E1 != P->Efinal) free(P->E1); // Part of the reason for checking this is to properly handle ctrl-c cases
  free(P->b);
  free(P->shapexs_transformed);
  free(P->shapeys_transformed);
  free(P->shapeRs_transformed);
  free(P->multiplier);
  #endif
  return;
}
