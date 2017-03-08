#ifndef UTILITIES_H
#define UTILITIES_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "postProcessing.h"
#include "structs.h"

void printVecD(double* vec, int n);
void printVecI(int* vec, int n);
void printMatI(int* mat, int n, int m);
void printMatD(double* mat, int n, int m);
void csvWriteD(double* mat, int n, int m, char* path);
void csvWriteI(int* mat, int n, int m, char* path);
void csvWriteLayer(double* mat, int n, int m, int k, char* path);
void writeTimeSeries(double* v, int n, char* path);
void* launchWriteThread(void* pdata_void);
void writeResults(FlowData* flow, LatticeConsts* lc, int iter, pthread_t* printThread, PrintData* pdata);
void csvWriteOmega(double* ux, double* uy, int n, int m, char* path);
void csvWritePres(double* rho, int n, int m, char* path);
void printProgression(int iter, int nIter);

#endif