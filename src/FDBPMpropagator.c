/********************************************
 * FDBPMpropagator.c, in the C programming language, written for MATLAB MEX function generation
 * 
 ** Compiling on Windows
 * Can be compiled with GCC using
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' -outdir private .\src\FDBPMpropagator.c ".\src\libut.lib" -R2018a"
 * ... or the Microsoft Visual C++ compiler (MSVC) with
 * "copyfile ./src/FDBPMpropagator.c ./src/FDBPMpropagator.cpp; mex COMPFLAGS='/Zp8 /GR /EHs /nologo /MD /openmp /W4 /WX /wd4204 /wd4100' -outdir private .\src\FDBPMpropagator.cpp ".\src\libut.lib" -R2018a"
 * 
 * The source code in this file is written is such a way that it is
 * compilable by either C or C++ compilers, either with GCC, MSVC or
 * the Nvidia CUDA compiler called NVCC, which is based on MSVC. To
 * compile with CUDA GPU acceleration support, you must have MSVC
 * installed. As of January 2020, mexcuda does not work with MSVC 2019,
 * so I'd recommend MSVC 2017. You also need the Parallel Computing
 * Toolbox, which you will find in the MATLAB addon manager. To compile, run:
 * "copyfile ./src/FDBPMpropagator.c ./src/FDBPMpropagator_CUDA.cu; mexcuda -llibut COMPFLAGS='-use_fast_math -res-usage $COMPFLAGS' -outdir private .\src\FDBPMpropagator_CUDA.cu -R2018a"
 *
 ** Compiling on macOS
 * As of March 2021, the macOS compiler doesn't support libut (for ctrl+c 
 * breaking) or openmp (for multithreading).
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -std=c11 -Wall' -outdir private ./src/FDBPMpropagator.c -R2018a"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 *
 ** Compiling on Linux
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' -outdir private ./src/FDBPMpropagator.c ./src/libut.so -R2018a"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Use a package manager like apt to install GCC (on Ubuntu, part of the build-essential package)
 * 2. Type "mex -setup" in the MATLAB command window
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
  #define CREALF(x) (x.real())
  #define CIMAGF(x) (x.imag())
  #define MAX(x,y) (max(x,y))
  #define MIN(x,y) (min(x,y))
  #define FLOORF(x) (floor(x))
  #include <nvml.h>
  #define TILE_DIM 32
#else
  #ifdef __GNUC__ // This is defined for GCC and CLANG but not for Microsoft Visual C++ compiler
    #define MAX(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _a: _b;})
    #define MIN(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _b: _a;})
    #include <complex.h>
    typedef float complex floatcomplex;
    #define CEXPF(x) (cexpf(x))
    #define CREALF(x) (crealf(x))
    #define CIMAGF(x) (cimagf(x))
    #define FLOORF(x) (floorf(x))
  #else
    #include <algorithm>
    #include <complex>
    typedef std::complex<float> floatcomplex;
    #define I std::complex<float>{0,1}
    #define CEXPF(x) (std::exp(x))
    #define CREALF(x) (std::real(x))
    #define CIMAGF(x) (std::imag(x))
    #define MAX(x,y) (std::max(x,y))
    #define MIN(x,y) (std::min(x,y))
    #define FLOORF(x) (std::floor(x))
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
  float dz;
  long iz_start;
  long iz_end;
  float taperPerStep;
  float twistPerStep;
  float d;
  float n_0;
  floatcomplex *n_in;
  long  Nx_n;
  long  Ny_n;
  long  Nz_n;
  float dz_n;
  floatcomplex *Efinal;
  floatcomplex *E1;
  floatcomplex *E2;
  floatcomplex *Eyx;
  floatcomplex *n_out;
  floatcomplex *b;
  float *multiplier;
  floatcomplex ax;
  floatcomplex ay;
  float rho_e;
  float RoC;
  float sinBendDirection;
  float cosBendDirection;
  double precisePower;
  double EfieldPower;
};

#ifdef __NVCC__ // If compiling for CUDA
__host__ __device__
#endif
float sqrf(float x) {return x*x;}

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
__device__ double atomicAdd(double* address, double val) {
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

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
  double EfieldPowerThread = 0;
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
      EfieldPowerThread += sqrf(CREALF(P->E2[i])) + sqrf(CIMAGF(P->E2[i]));
    }
  }

  atomicAdd(&P_global->EfieldPower,EfieldPowerThread);

  #else
  double EfieldPowerThread = 0;
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
      EfieldPowerThread += sqrf(CREALF(P->E2[i])) + sqrf(CIMAGF(P->E2[i]));
    }
  }
  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  P->EfieldPower += EfieldPowerThread;
  #ifdef _OPENMP
  #pragma omp barrier
  #endif
  #endif
}

#ifdef __NVCC__ // 1 If compiling for CUDA
__global__
#endif // 1
void applyMultiplier(struct parameters *P_global, long iz, struct debug *D) {
  double precisePowerDiffThread = 0;
  #ifdef __NVCC__ // 1
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying
  __syncthreads(); // All threads in the block wait for the copy to have finished
  float fieldCorrection = sqrtf((float)(P->precisePower/P->EfieldPower));
  float cosvalue = cosf(-P->twistPerStep*iz); // Minus is because we go from the rotated frame to the source frame
  float sinvalue = sinf(-P->twistPerStep*iz);
  float scaling = 1/(1 - P->taperPerStep*iz); // Take reciprocal because we go from scaled frame to unscaled frame
  for(long i=threadNum;i<P->Nx*P->Ny;i+=gridDim.x*blockDim.x*blockDim.y) {
  #else // 1
  struct parameters *P = P_global;
  float fieldCorrection = sqrtf((float)(P->precisePower/P->EfieldPower));
  float cosvalue = cosf(-P->twistPerStep*iz); // Minus is because we go from the rotated frame to the source frame
  float sinvalue = sinf(-P->twistPerStep*iz);
  float scaling = 1/(1 - P->taperPerStep*iz); // Take reciprocal because we go from scaled frame to unscaled frame
  #ifdef _OPENMP // 2
  #pragma omp for schedule(dynamic)
  #endif // 2
  for(long i=0;i<P->Nx*P->Ny;i++) {
  #endif // 1
    long ix = i%P->Nx;
    float x = P->dx*(ix - (P->Nx-1)/2.0f);
    long iy = i/P->Nx;
    float y = P->dy*(iy - (P->Ny-1)/2.0f);
    floatcomplex n = 0;
    if(P->taperPerStep || P->twistPerStep) { // Rotate, scale, interpolate. If we are tapering or twisting, we know that the RIP is 2D
      float x_src = scaling*(cosvalue*x - sinvalue*y);
      float y_src = scaling*(sinvalue*x + cosvalue*y);
      float ix_src = MIN(MAX(0.0f,x_src/P->dx + (P->Nx - 1)/2.0f),(P->Nx - 1)*(1-FLT_EPSILON)); // Fractional index, coerced to be within the source window
      float iy_src = MIN(MAX(0.0f,y_src/P->dy + (P->Ny - 1)/2.0f),(P->Ny - 1)*(1-FLT_EPSILON));
      long ix_low = (long)FLOORF(ix_src);
      long iy_low = (long)FLOORF(iy_src);
      float ix_frac = ix_src - FLOORF(ix_src);
      float iy_frac = iy_src - FLOORF(iy_src);
      n = P->n_in[ix_low     + P->Nx*(iy_low    )]*(1 - ix_frac)*(1 - iy_frac) + 
          P->n_in[ix_low + 1 + P->Nx*(iy_low    )]*(    ix_frac)*(1 - iy_frac) +
          P->n_in[ix_low     + P->Nx*(iy_low + 1)]*(1 - ix_frac)*(    iy_frac) + 
          P->n_in[ix_low + 1 + P->Nx*(iy_low + 1)]*(    ix_frac)*(    iy_frac); // Bilinear interpolation
    } else if(P->Nz_n == 1) { // 2D RIP
      n = P->n_in[i];
    } else { // 3D RIP
      float z = iz*P->dz;
      long ix_n = MIN(MAX(0L,ix - (P->Nx - P->Nx_n)/2),P->Nx_n-1);
      long iy_n = MIN(MAX(0L,iy - (P->Ny - P->Ny_n)/2),P->Ny_n-1);
      float iz_n = MIN(MAX(0.0f,z/P->dz_n),(P->Nz_n - 1)*(1-FLT_EPSILON)); // Fractional index, coerced to be within the n window
      long iz_n_low = (long)FLOORF(iz_n);
      float iz_n_frac = iz_n - FLOORF(iz_n);
      n = P->n_in[ix_n + P->Nx_n*iy_n + P->Ny_n*P->Nx_n*(iz_n_low    )]*(1 - iz_n_frac) + 
          P->n_in[ix_n + P->Nx_n*iy_n + P->Ny_n*P->Nx_n*(iz_n_low + 1)]*(    iz_n_frac); // Linear interpolation in z
    }
    if(iz == P->iz_end-1) P->n_out[i] = n;
    float n_bend = CREALF(n)*(1-(sqrf(CREALF(n))*(x*P->cosBendDirection+y*P->sinBendDirection)/2/P->RoC*P->rho_e))*exp((x*P->cosBendDirection+y*P->sinBendDirection)/P->RoC);
    floatcomplex a = P->multiplier[i]*CEXPF(P->d*(CIMAGF(n) + (sqrf(n_bend) - sqrf(P->n_0))*I/(2*P->n_0))); // Multiplier includes only the edge absorber
    P->E2[i] *= fieldCorrection*a;
    float anormsqr = sqrf(CREALF(a)) + sqrf(CIMAGF(a));
    if(anormsqr > 1 - 10*FLT_EPSILON && anormsqr < 1 + 10*FLT_EPSILON) anormsqr = 1; // To avoid accumulating power discrepancies due to rounding errors
    precisePowerDiffThread += (sqrf(CREALF(P->E2[i])) + sqrf(CIMAGF(P->E2[i])))*(1 - 1/anormsqr);
  }

  #ifdef __NVCC__ // 1
  atomicAdd(&P_global->precisePower,precisePowerDiffThread);
  #else // 1
  #ifdef _OPENMP // 2
  #pragma omp atomic
  #endif // 2
  P->precisePower += precisePowerDiffThread;
  #ifdef _OPENMP // 2
  #pragma omp barrier
  #endif // 2
  #endif // 1
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void swapEPointers(struct parameters *P, long iz) {
  #ifdef __NVCC__
  P->EfieldPower = 0;
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
  long N_n = P->Nx_n*P->Ny_n*P->Nz_n;
  struct parameters P_tempvar = *P;

  gpuErrchk(cudaMalloc(&P_tempvar.E1,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMemcpy(P_tempvar.E1,P->E1,N*sizeof(floatcomplex),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.E2,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMalloc(&P_tempvar.Eyx,N*sizeof(floatcomplex)));

  gpuErrchk(cudaMalloc(&P_tempvar.n_out,N*sizeof(floatcomplex)));

  gpuErrchk(cudaMalloc(&P_tempvar.b,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMalloc(&P_tempvar.multiplier,N*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.multiplier,P->multiplier,N*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.n_in,N_n*sizeof(floatcomplex)));
  gpuErrchk(cudaMemcpy(P_tempvar.n_in,P->n_in,N_n*sizeof(floatcomplex),cudaMemcpyHostToDevice));

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
  gpuErrchk(cudaMemcpy(P->n_out,P_temp.n_out,N*sizeof(floatcomplex),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(P_temp.n_out));
  P->precisePower = P_temp.precisePower;

  gpuErrchk(cudaFree(P_temp.E1));
  gpuErrchk(cudaFree(P_temp.E2));
  gpuErrchk(cudaFree(P_temp.Eyx));
  gpuErrchk(cudaFree(P_temp.b));
  gpuErrchk(cudaFree(P_temp.multiplier));
  gpuErrchk(cudaFree(P_temp.n_in));
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
  P->dz = *(float *)mxGetData(mxGetField(prhs[1],0,"dz"));
  P->iz_start = *(long *)mxGetData(mxGetField(prhs[1],0,"iz_start"));
  P->iz_end = *(long *)mxGetData(mxGetField(prhs[1],0,"iz_end"));
  P->taperPerStep = *(float *)mxGetData(mxGetField(prhs[1],0,"taperPerStep"));
  P->twistPerStep = *(float *)mxGetData(mxGetField(prhs[1],0,"twistPerStep"));
  P->d = *(float *)mxGetData(mxGetField(prhs[1],0,"d"));
  P->n_0 = *(float *)mxGetData(mxGetField(prhs[1],0,"n_0"));
  P->n_in = (floatcomplex *)mxGetData(mxGetField(prhs[1],0,"n_mat"));
  mwSize nDims = mxGetNumberOfDimensions(mxGetField(prhs[1],0,"n_mat"));
  mwSize const *dimPtr = mxGetDimensions(mxGetField(prhs[1],0,"n_mat"));
  P->Nx_n = (long)dimPtr[0];
  P->Ny_n = (long)dimPtr[1];
  P->Nz_n = nDims > 2? (long)dimPtr[2]: 1;
  P->dz_n = *(float *)mxGetData(mxGetField(prhs[1],0,"dz_n"));
  P->rho_e = *(float *)mxGetData(mxGetField(prhs[1],0,"rho_e"));
  P->RoC = *(float *)mxGetData(mxGetField(prhs[1],0,"RoC"));
  P->sinBendDirection = sin(*(float *)mxGetData(mxGetField(prhs[1],0,"bendDirection"))/180*PI);
  P->cosBendDirection = cos(*(float *)mxGetData(mxGetField(prhs[1],0,"bendDirection"))/180*PI);
  P->E1 = (floatcomplex *)mxGetData(prhs[0]); // Input E field
  dimPtr = mxGetDimensions(prhs[0]);
  P->Efinal = (floatcomplex *)mxGetData(plhs[0] = mxCreateNumericArray(2,dimPtr,mxSINGLE_CLASS,mxCOMPLEX)); // Output E field
  P->n_out = (floatcomplex *)mxGetData(plhs[1] = mxCreateNumericArray(2,dimPtr,mxSINGLE_CLASS,mxCOMPLEX)); // Output refractive index
  P->precisePower = *(double *)mxGetData(mxGetField(prhs[1],0,"inputPrecisePower"));
  #ifndef __NVCC__
  P->E2 = (floatcomplex *)((P->iz_end - P->iz_start)%2? P->Efinal: malloc(P->Nx*P->Ny*sizeof(floatcomplex)));
  #endif
  P->multiplier = (float *)mxGetData(mxGetField(prhs[1],0,"multiplier")); // Array of multiplier values to apply to the E field after each step, due to the edge absorber outside the main simulation window
  P->ax = *(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ax"));
  P->ay = *(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ay"));
  
  bool ctrlc_caught = false;      // Has a ctrl+c been passed from MATLAB?
  P->EfieldPower = 0;
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
      applyMultiplier<<<nBlocks, blockDims>>>(P_dev,iz,D_dev); // xy -> xy
      gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
      #else
      substep1a(P);
      substep1b(P);
      substep2a(P);
      substep2b(P);
      applyMultiplier(P,iz,NULL);
      #endif

      #ifdef _OPENMP
      #pragma omp master
      #endif
      {
        P->EfieldPower = 0;
        
        #ifndef __clang__
        if(utIsInterruptPending()) {
          ctrlc_caught = true;
          printf("\nCtrl+C detected, stopping.\n");
        }
        #endif

        if(iz+1 < P->iz_end) {
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
  if(P->E1 != mxGetData(prhs[0]) && P->E1 != P->Efinal) free(P->E1); // Part of the reason for checking this is to properly handle ctrl-c cases
  free(P->b);
  #endif
  double *outputPrecisePowerPtr = (double *)mxGetData(plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL));
  *outputPrecisePowerPtr = P->precisePower;
  return;
}
