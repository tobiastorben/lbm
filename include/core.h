void updateU(double* ux, double* uy, double* fIn, double* rho,
			double* ex, double* ey,int nx, int ny, int nf);
void updateRho(double* rho, double* fIn, int nx, int ny, int nf);
void collide(double* ux, double* uy, double*fIn, double* fOut, double* rho, double* ex,
				  double* ey, int nx,int ny, int nf, double tau, double* w);

void stream(double* fIn, double* fOut,int* ex,int* ey, int nx, int ny, int nf);