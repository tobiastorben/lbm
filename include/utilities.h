#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

void printVecD(double* vec, int n);
void printVecI(int* vec, int n);
void printMatI(int* mat, int n, int m);
void printMatD(double* mat, int n, int m);
void csvWriteD(double* mat, int n, int m, char* path);
void csvWriteI(int* mat, int n, int m, char* path);
void csvWriteLayer(double* mat, int n, int m, int k, char* path);