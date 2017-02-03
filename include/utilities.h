int* readObstacleData(int nBBcells);
void printVecD(double* vec, int n);
void printVecI(int* vec, int n);
void printMatI(int** mat, int n, int m);
void printMatD(double** mat, int n, int m);
//double*** initF(int nf,int nx, int ny, double uIn);
double** initUx(int nx,int ny, double uIn);
double** initUy(int nx,int ny);