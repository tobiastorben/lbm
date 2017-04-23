#ifndef OUTPUT_H
#define OUTPUT_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "structs.h"

void writeTimeSeries(double* v, int n, char* path);
void* launchWriteThread(void* pdata_void);
void writeResults(FlowData* flow, LatticeConsts* lc, int iter, pthread_t* printThread, PrintData* pdata);
void writeVorticity(double* ux, double* uy, int n, int m, double dt, char* path);
void writePres(double* rho, int n, int m, double c, double rhoPhys, char* path);
void writeU(double* mat, int n, int m, double c, char* path);
void printProgression(int iter, int nIter);
double* calcF(int ny, SimParams* params, double* rho, double* ux, double* uy);

#endif