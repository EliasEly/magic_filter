//! @file
//!  Magic filter not optimised

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "magic.h"

#define __START_TRACE() { asm volatile (".inst 0x2520e020"); }
#define __STOP_TRACE() { asm volatile (".inst 0x2520e040"); }

const double filter[] __attribute__ ((aligned (16))) = { 8.4334247333529341094733325815816e-7,
                       -0.1290557201342060969516786758559028e-4,
                        0.8762984476210559564689161894116397e-4,
                       -0.30158038132690463167163703826169879e-3,
                        0.174723713672993903449447812749852942e-2,
                       -0.942047030201080385922711540948195075e-2,
                        0.2373821463724942397566389712597274535e-1,
                        0.612625895831207982195380597e-1,
                        0.9940415697834003993178616713,
                       -0.604895289196983516002834636e-1,
                       -0.2103025160930381434955489412839065067e-1,
                        0.1337263414854794752733423467013220997e-1,
                       -0.344128144493493857280881509686821861e-2,
                        0.49443227688689919192282259476750972e-3,
                       -0.5185986881173432922848639136911487e-4,
                        2.72734492911979659657715313017228e-6};
const double filter_u[] __attribute__ ((aligned (16))) = { 2.72734492911979659657715313017228e-6,
                        8.4334247333529341094733325815816e-7,
                       -0.1290557201342060969516786758559028e-4,
                        0.8762984476210559564689161894116397e-4,
                       -0.30158038132690463167163703826169879e-3,
                        0.174723713672993903449447812749852942e-2,
                       -0.942047030201080385922711540948195075e-2,
                        0.2373821463724942397566389712597274535e-1,
                        0.612625895831207982195380597e-1,
                        0.9940415697834003993178616713,
                       -0.604895289196983516002834636e-1,
                       -0.2103025160930381434955489412839065067e-1,
                        0.1337263414854794752733423467013220997e-1,
                       -0.344128144493493857280881509686821861e-2,
                        0.49443227688689919192282259476750972e-3,
                       -0.5185986881173432922848639136911487e-4,
                        2.72734492911979659657715313017228e-6};
const double filter_reverse[] __attribute__ ((aligned (16))) = {
                        2.72734492911979659657715313017228e-6,
                       -0.5185986881173432922848639136911487e-4,
                        0.49443227688689919192282259476750972e-3,
                       -0.344128144493493857280881509686821861e-2,
                        0.1337263414854794752733423467013220997e-1,
                       -0.2103025160930381434955489412839065067e-1,
                       -0.604895289196983516002834636e-1,
                        0.9940415697834003993178616713,
                        0.612625895831207982195380597e-1,
                        0.2373821463724942397566389712597274535e-1,
                       -0.942047030201080385922711540948195075e-2,
                        0.174723713672993903449447812749852942e-2,
                       -0.30158038132690463167163703826169879e-3,
                        0.8762984476210559564689161894116397e-4,
                       -0.1290557201342060969516786758559028e-4,
                        8.4334247333529341094733325815816e-7
};
const double filter_reverse_u[] __attribute__ ((aligned (16))) = {
                        8.4334247333529341094733325815816e-7,
                        2.72734492911979659657715313017228e-6,
                       -0.5185986881173432922848639136911487e-4,
                        0.49443227688689919192282259476750972e-3,
                       -0.344128144493493857280881509686821861e-2,
                        0.1337263414854794752733423467013220997e-1,
                       -0.2103025160930381434955489412839065067e-1,
                       -0.604895289196983516002834636e-1,
                        0.9940415697834003993178616713,
                        0.612625895831207982195380597e-1,
                        0.2373821463724942397566389712597274535e-1,
                       -0.942047030201080385922711540948195075e-2,
                        0.174723713672993903449447812749852942e-2,
                       -0.30158038132690463167163703826169879e-3,
                        0.8762984476210559564689161894116397e-4,
                       -0.1290557201342060969516786758559028e-4,
                        8.4334247333529341094733325815816e-7
};

void magicfilter1d_naive_(int *n, int *ndat, double const *source, double *dest) {
  double tmp;
  unsigned int i,j,k;
  for(i=0;i<(*ndat);i++){
    for(j=0;j<(*n);j++) {
      tmp=0;
      for(k=0;k<16;k++){
        tmp+=source[(j-8+k+(*n))%(*n)]*filter[k];
      }
      dest[j*(*ndat)]=tmp;
    }
    dest += 1;
    source += (*n);
  } 
}

void magicfilter1d_naive_bis_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
  double tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  unsigned int i,j,k;
  for(i=0;i<(*ndat);i++){
    for(j=0;j<(*n);j++) {
      tmp=0, tmp1=0, tmp2=0, tmp3=0, tmp4=0, tmp5=0, tmp6=0, tmp7=0;
      for(k=0;k<16;k+=8){
        tmp+=source[(j-8+k+(*n))%(*n)]*filter[k];
        tmp1+=source[(j-8+k+1+(*n))%(*n)]*filter[k+1];
        tmp2+=source[(j-8+k+2+(*n))%(*n)]*filter[k+2];
        tmp3+=source[(j-8+k+3+(*n))%(*n)]*filter[k+3];
        tmp4+=source[(j-8+k+4+(*n))%(*n)]*filter[k+4];
        tmp5+=source[(j-8+k+5+(*n))%(*n)]*filter[k+5];
        tmp6+=source[(j-8+k+6+(*n))%(*n)]*filter[k+6];
        tmp7+=source[(j-8+k+7+(*n))%(*n)]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7;
    }
    dest += 1;
    source += (*n);
  } 
}

void magicfilter1d_naive_o1_(const int* restrict n, const int* restrict ndat,  const double* restrict source, double* restrict dest) {
  double tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  unsigned int i,j,k,index;
  for(i=0;i<(*ndat);i++){
    for(j=0;j<(*n);j++) {
      index = j-8+(*n);
      tmp=0, tmp1=0, tmp2=0, tmp3=0, tmp4=0, tmp5=0, tmp6=0, tmp7=0;
      for(k=0;k<16;k+=8){
        tmp+=source[(index+k)%(*n)]*filter[k];
        tmp1+=source[(index+k+1)%(*n)]*filter[k+1];
        tmp2+=source[(index+k+2)%(*n)]*filter[k+2];
        tmp3+=source[(index+k+3)%(*n)]*filter[k+3];
        tmp4+=source[(index+k+4)%(*n)]*filter[k+4];
        tmp5+=source[(index+k+5)%(*n)]*filter[k+5];
        tmp6+=source[(index+k+6)%(*n)]*filter[k+6];
        tmp7+=source[(index+k+7)%(*n)]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7;
    }
    dest += 1;
    source += (*n);
  } 
}

void magicfilter1d_naive_o2_(const int* restrict n, const int* restrict ndat,  const double* restrict source, double* restrict dest) {
  double tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  unsigned int i,j,k,index;
  for(i=0;i<(*ndat);i++){
    for(j=0; j<8; j++){
      index = j-8+(*n);
      tmp=0, tmp1=0, tmp2=0, tmp3=0, tmp4=0, tmp5=0, tmp6=0, tmp7=0;
      for(k=0;k<16;k+=8){
        tmp+=source[(index+k)%(*n)]*filter[k];
        tmp1+=source[(index+k+1)%(*n)]*filter[k+1];
        tmp2+=source[(index+k+2)%(*n)]*filter[k+2];
        tmp3+=source[(index+k+3)%(*n)]*filter[k+3];
        tmp4+=source[(index+k+4)%(*n)]*filter[k+4];
        tmp5+=source[(index+k+5)%(*n)]*filter[k+5];
        tmp6+=source[(index+k+6)%(*n)]*filter[k+6];
        tmp7+=source[(index+k+7)%(*n)]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7;
    }
    for(j=8;j<(*n)-8;j++) {
      index = j-8;
      tmp=0, tmp1=0, tmp2=0, tmp3=0, tmp4=0, tmp5=0, tmp6=0, tmp7=0;
      for(k=0;k<16;k+=8){
        tmp+=source[(index+k)]*filter[k];
        tmp1+=source[(index+k+1)]*filter[k+1];
        tmp2+=source[(index+k+2)]*filter[k+2];
        tmp3+=source[(index+k+3)]*filter[k+3];
        tmp4+=source[(index+k+4)]*filter[k+4];
        tmp5+=source[(index+k+5)]*filter[k+5];
        tmp6+=source[(index+k+6)]*filter[k+6];
        tmp7+=source[(index+k+7)]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7;
    }
    for(j=(*n)-8; j<(*n); j++){
      index = j-8+(*n);
      tmp=0, tmp1=0, tmp2=0, tmp3=0, tmp4=0, tmp5=0, tmp6=0, tmp7=0;
      for(k=0;k<16;k+=8){
        tmp+=source[(index+k)%(*n)]*filter[k];
        tmp1+=source[(index+k+1)%(*n)]*filter[k+1];
        tmp2+=source[(index+k+2)%(*n)]*filter[k+2];
        tmp3+=source[(index+k+3)%(*n)]*filter[k+3];
        tmp4+=source[(index+k+4)%(*n)]*filter[k+4];
        tmp5+=source[(index+k+5)%(*n)]*filter[k+5];
        tmp6+=source[(index+k+6)%(*n)]*filter[k+6];
        tmp7+=source[(index+k+7)%(*n)]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7;
    }

    dest += 1;
    source += (*n);
  } 
}

void magicfilter1d_naive_o3_(const int* restrict n, const int* restrict ndat,  const double* restrict source, double* restrict dest) {
	unsigned int* idx = malloc(((*n)+16)*sizeof(int));
	unsigned int x = 0;
	while(x < (*n)+16){
		*(idx+x)= (x-8+(*n))%(*n);
		x++;
	}

  double tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  unsigned int i,j,k,index;
  for(i=0;i<(*ndat);i++){
    for(j=0; j<8; j++){
      index = j;
      tmp=0, tmp1=0, tmp2=0, tmp3=0, tmp4=0, tmp5=0, tmp6=0, tmp7=0;
      for(k=0;k<16;k+=8){
        tmp+=source[idx[index+k]]*filter[k];
        tmp1+=source[idx[(index+k+1)]]*filter[k+1];
        tmp2+=source[idx[(index+k+2)]]*filter[k+2];
        tmp3+=source[idx[(index+k+3)]]*filter[k+3];
        tmp4+=source[idx[(index+k+4)]]*filter[k+4];
        tmp5+=source[idx[(index+k+5)]]*filter[k+5];
        tmp6+=source[idx[(index+k+6)]]*filter[k+6];
        tmp7+=source[idx[(index+k+7)]]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7;
    }
    for(j=8;j<(*n)-8;j++) {
      index = j-8;
      tmp=0, tmp1=0, tmp2=0, tmp3=0, tmp4=0, tmp5=0, tmp6=0, tmp7=0;
      for(k=0;k<16;k+=8){
        tmp+=source[(index+k)]*filter[k];
        tmp1+=source[(index+k+1)]*filter[k+1];
        tmp2+=source[(index+k+2)]*filter[k+2];
        tmp3+=source[(index+k+3)]*filter[k+3];
        tmp4+=source[(index+k+4)]*filter[k+4];
        tmp5+=source[(index+k+5)]*filter[k+5];
        tmp6+=source[(index+k+6)]*filter[k+6];
        tmp7+=source[(index+k+7)]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7;
    }
    for(j=(*n)-8; j<(*n); j++){
      index = j;
      tmp=0, tmp1=0, tmp2=0, tmp3=0, tmp4=0, tmp5=0, tmp6=0, tmp7=0;
      for(k=0;k<16;k+=8){
        tmp+=source[idx[index+k]]*filter[k];
        tmp1+=source[idx[(index+k+1)]]*filter[k+1];
        tmp2+=source[idx[(index+k+2)]]*filter[k+2];
        tmp3+=source[idx[(index+k+3)]]*filter[k+3];
        tmp4+=source[idx[(index+k+4)]]*filter[k+4];
        tmp5+=source[idx[(index+k+5)]]*filter[k+5];
        tmp6+=source[idx[(index+k+6)]]*filter[k+6];
        tmp7+=source[idx[(index+k+7)]]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7;
    }

    dest += 1;
    source += (*n);
  } 
}

void magicfilter1d_naive_o4_(const int* restrict n, const int* restrict ndat,  const double* restrict source, double* restrict dest) {
	int* idx = malloc(((*n)+16)*sizeof(int));
	int x = 0;
	while(x < (*n)+16){
		*(idx+x)= (x-8+(*n))%(*n);
		x++;
	}

  double tmp[8];
  int i,j,k,index;
  for(i=0;i<(*ndat);i++){ 
    for(j=0; j<8; j++){
      index = j;
      for(int i = 0; i<8; i++)
        tmp[i] = 0;
      for(k=0;k<16;k+=8){
        tmp[0]+=source[idx[(index+k)]]*filter[k];
        tmp[1]+=source[idx[(index+k+1)]]*filter[k+1];
        tmp[2]+=source[idx[(index+k+2)]]*filter[k+2];
        tmp[3]+=source[idx[(index+k+3)]]*filter[k+3];
        tmp[4]+=source[idx[(index+k+4)]]*filter[k+4];
        tmp[5]+=source[idx[(index+k+5)]]*filter[k+5];
        tmp[6]+=source[idx[(index+k+6)]]*filter[k+6];
        tmp[7]+=source[idx[(index+k+7)]]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
    }
    for(j=8;j<(*n)-8;j++) {
      index = j-8;
      for(int i = 0; i<8; i++)
        tmp[i] = 0;
      for(k=0;k<16;k+=8){
        tmp[0]+=source[index+k]*filter[k];
        tmp[1]+=source[(index+k+1)]*filter[k+1];
        tmp[2]+=source[(index+k+2)]*filter[k+2];
        tmp[3]+=source[(index+k+3)]*filter[k+3];
        tmp[4]+=source[(index+k+4)]*filter[k+4];
        tmp[5]+=source[(index+k+5)]*filter[k+5];
        tmp[6]+=source[(index+k+6)]*filter[k+6];
        tmp[7]+=source[(index+k+7)]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp[0]+ tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
    }
 for(j=(*n)-8; j<(*n); j++){
     index = j;
     for(int i = 0; i<8; i++)
        tmp[i] = 0;
      for(k=0;k<16;k+=8){
        tmp[0]+=source[idx[(index+k)]]*filter[k];
        tmp[1]+=source[idx[(index+k+1)]]*filter[k+1];
        tmp[2]+=source[idx[(index+k+2)]]*filter[k+2];
        tmp[3]+=source[idx[(index+k+3)]]*filter[k+3];
        tmp[4]+=source[idx[(index+k+4)]]*filter[k+4];
        tmp[5]+=source[idx[(index+k+5)]]*filter[k+5];
        tmp[6]+=source[idx[(index+k+6)]]*filter[k+6];
        tmp[7]+=source[idx[(index+k+7)]]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp[0]+ tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
    }
    dest += 1;
    source += (*n);
  }
}

void magicfilter1d_t_naive_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
  double tmp;
  unsigned int i,j,k;
  for(i=0;i<(*ndat);i++){
    for(j=0;j<(*n);j++) {
      tmp=0;
      for(k=0;k<16;k++){
        tmp+=source[(j-7+k+(*n))%(*n)]*filter_reverse[k];
      }
      dest[j*(*ndat)]=tmp;
    }

    dest += 1;
    source += (*n);
  } 
}

