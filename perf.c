#include <stdio.h>
#include <x86intrin.h>
#include <time.h>

// The cycle time
static const float cycle_duration = (float) 1 / (float) 3.1;

// How long takes the chrono to start and stop
static unsigned long long int residue;

void perf_init(){
    unsigned long long int start, end;

    start = _rdtsc();

    end = _rdtsc();

    residue = end - start;

    return;
}


void flop_compute(char* message, unsigned long long int nb_operation, unsigned long long int cycles){
    printf("%s %llu operations\n", message, nb_operation);
    printf("    %5.3lf GFLOP/s\n", ((double) nb_operation) / (((double) (cycles - residue)) * cycle_duration));
}
