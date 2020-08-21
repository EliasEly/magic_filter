#include <stdio.h>
#include <time.h>

// The cycle time
static const float cycle_duration = (float) 1 / (float) 3.1;

// How long takes the chrono to start and stop
static double residue;

void perf_init(){
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
    clock_gettime(CLOCK_REALTIME, &end);

    residue = (end.tv_sec- start.tv_sec)*1e9 + (end.tv_nsec-start.tv_nsec);

    return;
}


void flop_compute(char* message, unsigned long long int nb_operation, double cycles){
    printf("%s %llu operations\n", message, nb_operation);
    printf("    %5.3lf GFLOP/s\n", ((double) nb_operation) / (((double) (cycles - residue)) * cycle_duration));
}
