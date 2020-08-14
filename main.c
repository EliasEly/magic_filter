#include <stdlib.h>
#include <assert.h>
#include <math.h>

extern void magicfilter1d_sse_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);
extern void magicfilter1d_t_sse_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);
extern void magicfilter1d_naive_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);
extern void magicfilter1d_naive_o1_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);
extern void magicfilter1d_t_naive_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);

#define CONV_REF magicfilter1d_naive_o1_
#define CONV magicfilter1d_sse_ 
#define X 31
#define Y 32
#define Z 33
#define TOTAL (2*X)*(2*Y)*(2*Z)

void init_rand(double *arr, size_t dim) {
	for (size_t i = 0; i < dim; i++) {
		arr[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
	}
}

void check(double *arr_in, double *arr_out, size_t dim) {
	for (size_t i = 0; i < dim; i++) {
		assert(fabs(arr_in[i]-arr_out[i]) < 1e-14);
	}
}

int main() {
        
	double * data_in = malloc(sizeof(double) * TOTAL);
	double * data_tmp = calloc(sizeof(double), TOTAL);
	double * data_out = calloc(sizeof(double), TOTAL);
	double * data_out2 = calloc(sizeof(double), TOTAL);

	init_rand(data_in, TOTAL);

	int n, ndat;

	n = 2*X;
	ndat = 2*Y*2*Z;
	CONV(&n, &ndat, data_in, data_out);
	n = 2*Y;
	ndat = 2*Z*2*X;
	CONV(&n, &ndat, data_out, data_tmp);
	n = 2*Z;
	ndat = 2*X*2*Y;
	CONV(&n, &ndat, data_tmp, data_out);

	n = 2*X;
	ndat = 2*Y*2*Z;
	CONV_REF(&n, &ndat, data_in, data_out2);
	n = 2*Y;
	ndat = 2*Z*2*X;
	CONV_REF(&n, &ndat, data_out2, data_tmp);
	n = 2*Z;
	ndat = 2*X*2*Y;
	CONV_REF(&n, &ndat, data_tmp, data_out2);

	check(data_out, data_out2, TOTAL);
	return 0;
}
