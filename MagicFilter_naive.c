//! @file
//!  Magic filter not optimised

#include <time.h>
#include <stdio.h>
#include "magic_exec.h"

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

void magicfilter1d_naive_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
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

//couper boucle j en 3 (0 à 8) (8 à n-8) ...
//(8 à n-8) pas besoin de modulo
//(0 à 8) (n-8 à n) possible de calculer les indices

void magicfilter1d_naive_bis_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
  //should create an array of size n which would contain value ranging from 0 to n-1 (-8 à n+8) -> rajouter 8 pour indice valide
  //goal->avoid modulo computation in k loop 
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
//module empeche vect -> s'en debarasser 
//256 vecteur donc pas de 4 (double-64 bits)

//utiliser clock_gettime pour mesurer perf
void magicfilter1d_naive_o1_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
  //should create an array of size n which would contain value ranging from 0 to n-1
  //goal->avoid modulo computation in k loop 

  //sortir les 8 premières et 8 dernières itérations dans j-loop 
  //derouler i_loop (une ou deux fois)
  //derouler j_loop (pas aussi important)
  double tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  unsigned int i,j,k,index;
  for(i=0;i<(*ndat);i++){
    for(j=0;j<(*n);j++) {
      index = j-8+(*n);
      tmp=0, tmp1=0, tmp2=0, tmp3=0, tmp4=0, tmp5=0, tmp6=0, tmp7=0;
      for(k=0;k<16;k+=8){
        tmp+=source[index+k%(*n)]*filter[k];
        tmp1+=source[index+k+1%(*n)]*filter[k+1];
        tmp2+=source[index+k+2%(*n)]*filter[k+2];
        tmp3+=source[index+k+3%(*n)]*filter[k+3];
        tmp4+=source[index+k+4%(*n)]*filter[k+4];
        tmp5+=source[index+k+5%(*n)]*filter[k+5];
        tmp6+=source[index+k+6%(*n)]*filter[k+6];
        tmp7+=source[index+k+7%(*n)]*filter[k+7];
      }
      dest[j*(*ndat)]=tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7;
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

