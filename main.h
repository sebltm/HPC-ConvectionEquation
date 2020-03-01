//
// Created by 660046669 on 2020/02/10
//

#ifndef HIGH_PERFORMANCE_COMPUTING_MAIN_H
#define HIGH_PERFORMANCE_COMPUTING_MAIN_H

#define SIZE_SIDE 1000
#define DELTA_X 0.001
#define DELTA_Y 0.001
#define TIME_STEP 0.05
#define FLUID_VELOCITY 0.01
#define IO_WRITE 1
#define TIMES_WRITE 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main();
double** initialiseArray();
void destroyArray(double **);
double fieldRateOfChange(double **, double **);
double addArrays(double **, double **);

#endif //HIGH_PERFORMANCE_COMPUTING_MAIN_H
