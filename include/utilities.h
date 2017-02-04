int* readObstacleData(int nBBcells);
void printVecD(double* vec, int n);
void printVecI(int* vec, int n);
void printMatI(int** mat, int n, int m);
void printMatD(double** mat, int n, int m);
double*** initFin(int nf,int nx, int ny, int* ex, int* ey,
				double** ux, double** uy, double* w, double** rho);
double*** initFout(int nx,int ny, int nf);
double** initUx(int nx,int ny, double uIn);
double** initUy(int nx,int ny);
double** initRho(int nx, int ny, double rho);
void updateU(double** ux, double** uy, double*** fIn, double** rho,
			int* ex, int* ey,int nx, int ny, int nf);
void updateRho(double** rho, double*** fIn, int nx, int ny, int nf);
void inlet(double** ux, double** uy, double** rho, double*** fIn, double uIn, int nx, int ny);
void outlet(double** ux, double** uy, double** rho, double*** fIn, int nx, int ny);
void collide(double** ux, double** uy, double***fIn, double*** fOut, double** rho, int* ex,
				  int* ey, int nx,int ny, int nf, double tau, double* w);
void bounce(double*** fIn, double*** fOut,int nx,int ny,int nf, int** bbCells, int nBBcells, int* opposite);