#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include <stdlib.h>

double* calcF(double* F, double* rho, int nx, int ny, int* bbCells, int nBBcells, int* bbCellMat);

#endif