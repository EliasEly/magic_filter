/*
* This file provides the prototype for magicfilter.
* Useful for the calls in main.c
*/

void magicfilter1d_naive_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);

void magicfilter1d_naive_bis_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);

void magicfilter1d_t_naive_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);

void magicfilter1d_naive_o1_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);

void magicfilter1d_naive_o2_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);

void magicfilter1d_naive_o3_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);

void magicfilter1d_sse_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);

void magicfilter1d_t_sse_(unsigned int *n, unsigned int *ndat, double const *source, double *dest);
