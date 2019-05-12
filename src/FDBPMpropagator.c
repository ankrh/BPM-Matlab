/********************************************
 * FDBPMpropagator.c, in the C programming language, written for MATLAB MEX function generation
 * Can be compiled using "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' .\src\FDBPMpropagator_floats.c ".\src\libut.lib" -R2018a"
 ********************************************/

#include <math.h>
#include <stdint.h>
#include <complex.h>
#include "mex.h"
#include "omp.h"
extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#define E1 (*E_src)
#define E2 (*E_tgt)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
	double complex *E_in = mxGetData(prhs[0]); // Input E field
    long nx = mxGetM(prhs[0]);
    long ny = mxGetN(prhs[0]);
    long nz = *mxGetDoubles(mxGetField(prhs[1],0,"nz")); // Number of steps to perform
	double complex *multiplier = mxGetData(mxGetField(prhs[1],0,"multiplier")); // Array of multiplier values to apply to the E field after each step, which includes (1) the phase imparted by the refractive index differences, (2) the absorber outside the fiber and (3) the effects of fibre bending, if present
    double complex ax = *(double complex *)mxGetData(mxGetField(prhs[1],0,"ax"));
    double complex ay = *(double complex *)mxGetData(mxGetField(prhs[1],0,"ay"));
	
    double complex *E_out = mxGetData(plhs[0] = mxCreateDoubleMatrix(nx,ny,mxCOMPLEX)); // Output E field

	double complex *E_temp = malloc(nx*ny*sizeof(double complex)); // Temporary E matrix
    
    double complex **E_src = &E_in; // Pointer to the pointer to whichever E matrix is to be read from
    double complex **E_tgt; // Pointer to the pointer to whichever E matrix is to be written to
    E_tgt = (nz%2)? &E_out: &E_temp; // If nz is odd then we can start by writing to the output E matrix (pointed to by plhs[0]), otherwise we start by writing to the temporary E matrix
    
    bool useAllCPUs = mxIsLogicalScalarTrue(mxGetField(prhs[1],0,"useAllCPUs"));
	bool ctrlc_caught = false;           // Has a ctrl+c been passed from MATLAB?
    #pragma omp parallel num_threads(useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1)
    {
        long i,ix,iy,iz;
        double complex w, b[nx>ny? nx: ny];
		
        for(iz=0; iz<nz; iz++) {
            if(ctrlc_caught) break;
			
            // Explicit part of substep 1 out of 3
            #pragma omp for schedule(auto)
			for(iy=0; iy<ny; iy++) {
				for(ix=0; ix<nx; ix++) {
					i = ix + iy*nx;

					E2[i] = E1[i];
					if(ix != 0   ) E2[i] += (E1[i-1]  - E1[i])*ax;
					if(ix != nx-1) E2[i] += (E1[i+1]  - E1[i])*ax;
					if(iy != 0   ) E2[i] += (E1[i-nx] - E1[i])*ay*2;
					if(iy != ny-1) E2[i] += (E1[i+nx] - E1[i])*ay*2;
				}
			}
			
			// Implicit part of substep 1 out of 3
            #pragma omp for schedule(auto)
			for(iy=0; iy<ny; iy++) {
				// Thomson algorithm, sweeps up from 0 to nx-1 and then down from nx-1 to 0:
				// Algorithm is taken from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
				for(ix=0; ix<nx; ix++) {
					b[ix] = 1;
					if(ix < nx-1) b[ix] += ax;
					if(ix > 0) {
						b[ix] +=  ax;
						w      = -ax/b[ix-1];
						b[ix] += w*ax;
						i = ix + iy*nx;
						E2[i] -= w*E2[i-1];
					}
				}

				for(ix=nx-1; ix>=0; ix--) {
					i = ix + iy*nx;
					E2[i] = (E2[i] + (ix == nx-1? 0: ax*E2[i+1]))/b[ix];
				}
			}
			
            // Explicit part of substep 2 out of 3
            #pragma omp for schedule(auto)
			for(iy=0; iy<ny; iy++) {
				for(ix=0; ix<nx; ix++) {
					i = ix + iy*nx;

					if(iy != 0   ) E2[i] -= (E1[i-nx] - E1[i])*ay;
					if(iy != ny-1) E2[i] -= (E1[i+nx] - E1[i])*ay;
				}
			}
			
			// Implicit part of substep 2 out of 3
            #pragma omp for schedule(auto)
			for(ix=0; ix<nx; ix++) {
				for(iy=0; iy<ny; iy++) {
					b[iy] = 1;
					if(iy < ny-1) b[iy] += ay;
					if(iy > 0) {
						b[iy] +=  ay;
						w      = -ay/b[iy-1];
						b[iy] += w*ay;
						i = ix + iy*nx;
						E2[i] -= w*E2[i-nx];
					}
				}

				for(iy=ny-1; iy>=0; iy--) {
					i = ix + iy*nx;
					E2[i] = (E2[i] + (iy == ny-1? 0: ay*E2[i+nx]))/b[iy];
				}
			}
			
            #pragma omp for schedule(auto)
			for(iy=0; iy<ny; iy++) {
				for(ix=0; ix<nx; ix++) {
					i = ix + iy*nx;
					E2[i] *= multiplier[i];
				}
			}

			#pragma omp master
			{
				if(utIsInterruptPending()) {
					ctrlc_caught = true;
					printf("\nCtrl+C detected, stopping.\n");
				}

				if(E_tgt == &E_out) {
					E_tgt = &E_temp;
					E_src = &E_out;
				} else {
					E_tgt = &E_out;
					E_src = &E_temp;
				}
			}
			#pragma omp barrier
		}
    }
    free(E_temp);
    return;
}
