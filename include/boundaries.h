void inlet(double* ux, double* uy, double* rho, double* fIn, double uIn, int nx, int ny);
void outlet(double* ux, double* uy, double* rho, double* fIn, int nx, int ny);
void bounce(double* fIn, double* fOut,int nx,int ny,int nf, int* bbCells, int nBBcells, int* opposite);