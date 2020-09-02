//! @file
//!  Performance test 
//!
//! @author
//!     Elias El Yandouzi - 2020

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "perf.h"
#include "magic.h"

#define NLOOP 10

#define CONV magicfilter1d_naive_o3_
#define CONV_REF magicfilter1d_naive_ 


#define SIZE 65536
#define X 31
#define Y 32
#define Z 33
#define TOTAL (2*X)*(2*Y)*(2*Z)


double* init_vector(unsigned int n){
    double* res = malloc(n*sizeof(double));

    double random;

    for(int i = 0; i < n; i++){
        random = (double)(rand() % 100) / (double)100;
        res[i] = random;
        #ifdef DEBUG
            printf("Random value in init_vector (main.c:21) : %lu", random);
        #endif
    }

    return res;
}


void check(double *arr_in, double *arr_out, size_t dim) {
	for (size_t i = 0; i < dim; i++) {
		assert(fabs(arr_in[i]-arr_out[i]) < 1e-14);
	}
}

int main(int argc, char** argv){
    
    time_t t;
    srand((unsigned) time(&t));
    
    unsigned int n, ndat;
	double * data_in = calloc(sizeof(double), TOTAL);
	double * data_out = calloc(sizeof(double), TOTAL);

	n = 2*X;
	ndat = 2*Y*2*Z;

	printf("- Running perf test -\n");
	struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
   	for(int i = 0; i < NLOOP; i++){
		magicfilter1d_naive_(&n, &ndat, data_in, data_out);
	}
	clock_gettime(CLOCK_REALTIME, &end);
	flop_compute("MagicFiler_naive_ : ", NLOOP*(8*16*2+2*(n-16))*ndat, (end.tv_sec- start.tv_sec)*1e9 + (end.tv_nsec-start.tv_nsec));


	clock_gettime(CLOCK_REALTIME, &start);
   	for(int i = 0; i < NLOOP; i++){
		magicfilter1d_naive_o1_(&n, &ndat, data_in, data_out);
	}
	clock_gettime(CLOCK_REALTIME, &end);
	flop_compute("MagicFiler_naive_o1_ : ", NLOOP*(8*16*2+2*(n-16))*ndat, (end.tv_sec- start.tv_sec)*1e9 + (end.tv_nsec-start.tv_nsec));

	clock_gettime(CLOCK_REALTIME, &start);
   	for(int i = 0; i < NLOOP; i++){
		magicfilter1d_naive_o2_(&n, &ndat, data_in, data_out);
	}
	clock_gettime(CLOCK_REALTIME, &end);
	flop_compute("MagicFiler_naive_o2_ : ", NLOOP*(8*16*2+2*(n-16))*ndat, (end.tv_sec- start.tv_sec)*1e9 + (end.tv_nsec-start.tv_nsec));

	clock_gettime(CLOCK_REALTIME, &start);
   	for(int i = 0; i < NLOOP; i++){
		magicfilter1d_naive_o3_(&n, &ndat, data_in, data_out);
	}
	clock_gettime(CLOCK_REALTIME, &end);
	flop_compute("MagicFiler_naive_o3_ : ", NLOOP*(8*16*2+2*(n-16))*ndat, (end.tv_sec- start.tv_sec)*1e9 + (end.tv_nsec-start.tv_nsec));

    return 0;
}
