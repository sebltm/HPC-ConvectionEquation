//
// Created by 660046669 on 2020/02/10
//

#include "main.h"
#include <omp.h>

int main() {
    // Initiate time variables and start the timer
    double tstart, tend, timeRoC, timeAdd;
    tstart = omp_get_wtime();

    // Initialise other variables and the arrays
    int currentStep, i, j = 0;
    double **scalarField = initialiseArray();
    double **rateOfChange = initialiseArray();
    int maxTimeSteps = 1500;


    // If either of our arrays couldn't be allocated, exit
    if(scalarField == NULL || rateOfChange == NULL) {
        fprintf(stderr, "Error: memory allocation for the scalar field failed\n");
        return -1;
    }

    // Calculate for timestep 0 to maxTimeSteps
    for(currentStep = 0; currentStep < maxTimeSteps; currentStep++) {

        // Calculate rate of change in "rateOfChange" and add to "scalarField" array
        timeRoC = fieldRateOfChange(scalarField, rateOfChange);
        timeAdd = addArrays(scalarField, rateOfChange);

#if TIMES_WRITE
        // Output time taken to calculate scalar field change and adding arrays, for each time step
        // Not compiled if TIMES_WRITE flag is set to 0 in main.h
        printf("Time step %d\n, time for rate of change calc = %f, time for adding arrays = %f\n", currentStep, timeRoC, timeAdd);
#endif
    }

    // Stop the timer and output total time
    tend = omp_get_wtime();
    printf("Total time %f seconds with %d threads\n", tend - tstart, omp_get_num_threads());


#if IO_WRITE
    // Output final array
    // Not compiled if IO_WRITE flag is 0 in main.h

    FILE *fp;
    fp = fopen("output.dat", "w");

    for(i = 0; i < SIZE_SIDE; i++) {
        for(j = 0; j < SIZE_SIDE; j++) {
            fprintf(fp, "%f %f %f\n", (double)i/1000.0, (double)j/1000.0, scalarField[i][j]);
        }
    }

    fclose(fp);
#endif

    destroyArray(scalarField);
    destroyArray(rateOfChange);

    return 0;
}

double fieldRateOfChange(double **currentArray, double **newArray) {
    int i, j = 0;
    double xDim, yDim;
    double tstart, tend;

    tstart = omp_get_wtime();
#pragma omp parallel for schedule(dynamic) default(none) shared(currentArray, newArray) private(i, j, xDim, yDim)
    for(i = 0; i < SIZE_SIDE; i++) {
        for(j = 0; j < SIZE_SIDE; j++) {

            if((i - 1) < 0) {
                // Boundary condition on x
                xDim = FLUID_VELOCITY * (currentArray[i][j]) / DELTA_X;
            } else {
                // General condition of y
                xDim = FLUID_VELOCITY * (currentArray[i][j] - currentArray[i - 1][j]) / DELTA_X;
            }

            if((j - 1) < 0) {
                // Boundary condition on y
                yDim = FLUID_VELOCITY * (currentArray[i][j]) / DELTA_Y;

            } else {
                // General case of y
                yDim = FLUID_VELOCITY * (currentArray[i][j] - currentArray[i][j - 1]) / DELTA_Y;

            }

            newArray[i][j] = -xDim-yDim;
        }
    }

    tend = omp_get_wtime();
    return tend-tstart;
}

double addArrays(double **currentArray, double **newArray) {
    int i, j;
    double tstart, tend;

    tstart = omp_get_wtime();

    // Add the array elements in parallel
#pragma omp parallel for schedule(dynamic) default(none) shared(currentArray, newArray) private(i, j)
    for(i = 0; i < SIZE_SIDE; i++) {
        for(j = 0; j < SIZE_SIDE; j++) {
            currentArray[i][j] += (newArray[i][j] * TIME_STEP);
        }
    }

    return tend-tstart;
}

void destroyArray(double** array) {
    // Function to destroy the double pointer array

    free(array[0]);
    free(array);
}

double** initialiseArray() {
    /*
     * Initialise a 2D array with values corresponding to the initial conditions of the problem
     * The 2D array is created so that the entire array is stored in a contiguous memory block
     */

    int i, j;
    double x, y, x0, y0, width;
    double **array = (double **)malloc(sizeof(double *) * SIZE_SIDE);
    double *arrayData = (double *)malloc(sizeof(double) * SIZE_SIDE * SIZE_SIDE);

    // In case of error while allocating memory, return null and let main function handle error
    if(array == NULL || arrayData == NULL) {
        free(array);
        free(arrayData);
        return NULL;
    }

    // Initialise the "rows" of the array
    for(i = 0; i < SIZE_SIDE; i++) {
        array[i] = &(arrayData[i * SIZE_SIDE]);
    }

    x0 = y0 = 0.1;
    width = 0.03;

    // Initiate array in parallel
#pragma omp parallel for schedule(dynamic) default(none) shared(array, x0, y0, width) private(i, j, x, y)
    for(i = 0; i < SIZE_SIDE; i++) {
        for(j = 0; j < SIZE_SIDE; j++) {

            // x = i/1000 and y = i/1000
            x = DELTA_X * (double)(i + 1);
            y = DELTA_Y * (double)(j + 1);

            // Initial condition function
            double top = pow(x - x0, 2.0) + pow(y - y0, 2.0);
            double bottom = 2.0 * pow(width, 2.0);
            array[i][j] = exp(-top/bottom);
        }
    }

#if IO_WRITE
    // Output initial conditions
    // Not compiled if IO_WRITE flag is 0 in main.h

    FILE *fp;
    fp = fopen("initial.dat", "w");

    for(i = 0; i < SIZE_SIDE; i++) {
        for(j = 0; j < SIZE_SIDE; j++) {
            fprintf(fp, "%f %f %f\n", (double)i/1000.0, (double)j/1000.0, array[i][j]);
        }
    }

    fclose(fp);
#endif
    return array;
}