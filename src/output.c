#include "output.h"

//------------------------------------------------------------------------------
//LBM    Function: writeU  
//------------------------------------------------------------------------------
//PURPOSE:	Writes velocity component to file in csv format
//USAGE:	writeU(u,n,m,c,path)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			u		double*     	Matrix containing velocity compent at each
//									node
//			n     	int				Number of rows
//			m		int 			Number of colums
//			c		double			Scale factor between lattice units and
//									physical units
//			path	char*			Path to be written to
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void writeU(double* u, int n, int m, double c, char* path) {
	FILE* fp;
	fp = fopen(path,"w");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		//Convert velocity to physical units, and write to file
		fprintf(fp,"%.3e",c*u[m*i + j]);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}
//------------------------------------------------------------------------------
//LBM    Function: writeTimeSeries
//------------------------------------------------------------------------------
//PURPOSE:	Writes a vector of values to a row in a file
//USAGE:	writeTimeSeries(v,n,path)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			v		double*     	Vector to be written
//			n     	int				Number of rows
//			path	char*			Path to be written to
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void writeTimeSeries(double* v, int n, char* path) {
	FILE* fp = fopen(path,"a");
	//Print the vector to the row in a file, in CSV format
	for (int i = 0; i < (n-1); i++){
		fprintf(fp,"%.5e, ", v[i]);
		}
	fprintf(fp,"%.5e\n", v[n-1]);
	fclose(fp);
}

//------------------------------------------------------------------------------
//LBM    Function: writePres
//------------------------------------------------------------------------------
//PURPOSE:	Writes pressure field to file in csv format
//USAGE:	writePres(rho,n,m,c,rhoPhys,path)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			rho		double*     	Matrix containing the density at each
//									node, in lattice units
//			n     	int				Number of rows
//			m		int 			Number of colums
//			c		double			Scale factor between lattice units and
//									physical units
//			rhoPhys	double			Refrence density in physical units
//			path	char*			Path to be written to
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void writePres(double* rho, int n, int m, double c, double rhoPhys, char* path) {
	double cSq;
	FILE* fp;
	
	fp = fopen(path,"w");
	cSq = c*c;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		//Convert pressure to physical units, and write to file
		fprintf(fp,"%.3e", rhoPhys*cSq*(1.0/3.0)*(rho[m*i + j]-1));
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}

//------------------------------------------------------------------------------
//LBM    Function: writeVorticity
//------------------------------------------------------------------------------
//PURPOSE:	Calculates vorticity at domain interior by a central differance
//			approximation, and writes it to a file in csv format.
//USAGE:	writeVorticity(ux,uy,n,m,dt,path)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			ux		double*     	Matrix containing x component of velocity
//									at each node, in lattice units
//			uy		double*     	Matrix containing y component of velocity
//									at each node, in lattice units
//			n     	int				Number of rows
//			m		int 			Number of colums
//			dt		double			Time step in physical units physical units
//			path	char*			Path to be written to
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void writeVorticity(double* ux, double* uy, int n, int m, double dt, char* path) {
	FILE* fp;
	double vort,coeff;
	
	fp = fopen(path,"w");
	coeff = 1.0/(4.0*dt);
	
	for (int i = 0; i < (n-1); i++){
		for (int j = 0; j < (m-1); j++){
		//Calculate second order central difference approximation to vorticity
		vort = coeff*(uy[(i+1)*m +j] - uy[(i-1)*m+j] - (ux[i*m +j+1] - ux[i*m +j-1]));
		//Write to file
		fprintf(fp,"%.3e", vort);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}

//------------------------------------------------------------------------------
//LBM    Function: writeResults
//------------------------------------------------------------------------------
//PURPOSE:	Copy the values to be written to a new memory location, and launch a
//			thread to write them, while the execution of the solver continues in
//			parallell.
//USAGE:	writeResults(flow,lc,iter,printThread,pdata)
//ARGUMENTS:
//			Name 	 	Type     		Description
//.............................................................................
//			flow		FlowData*		The field variables of the flow
//			lc			LatticeConsts*	The constants of the D2Q9 lattice
//			iter    	int				Current iteration
//			printThread	pthread_t* 		Pointer to print thread
//			pdata		PrintData*		All the data needed for the print thread
//.............................................................................
//CALLS:				
//			launchWriteThread			Starts the thread that writes results
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void writeResults(FlowData* flow, LatticeConsts* lc, int iter, pthread_t* printThread, PrintData* pdata) {
	int nx,ny;
	//Wait if the last print thread has not finished yet. Do not wait if this is the first print
	if (iter != pdata->params->startWrite) pthread_join(*printThread,NULL);
	pdata->iter = iter;
	nx = lc->nx;
	ny = lc->ny;
	
	//Copy data to new memory locations, so they can be printed while the original matrices are updated
	memcpy(pdata->uxCpy,flow->ux,nx*ny*sizeof(double));
	memcpy(pdata->uyCpy,flow->uy,nx*ny*sizeof(double));
	memcpy(pdata->rhoCpy,flow->rho,nx*ny*sizeof(double));
	
	//Launch print thread
	pthread_create(printThread,NULL,launchWriteThread,pdata);
	
	return;
}

//------------------------------------------------------------------------------
//LBM    Function: launchWriteThread
//------------------------------------------------------------------------------
//PURPOSE:	Entry point for print thread. Checks which output should be written,
//			and calls the corresponding function to do the writing.
//USAGE:	launchWriteThread(flow,lc,iter,printThread,pdata)
//ARGUMENTS:
//			Name 	 	Type     		Description
//.............................................................................
//			pdata		PrintData*		All the data needed for the print thread
//.............................................................................
//CALLS:				
//			writeU						Writes velocity componen to file
//			writePres					Writes pressure field to file
//			writeVorticity				Write vorticity field to file
//			calcF						Calculte the lift and drag force on the
//										obstacle
//			writeTimeSeries				Write vector to a line in a file

//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void* launchWriteThread(void* pdata_void) {
	
	double *F,*ux,*uy,*rho,dx,dt,c,rhoPhys;
	int nx,ny,iter,*outputSelect;
	char *fName,*path,*outDir;
	PrintData* pdata;
	SimParams* params;
	
	pdata = (PrintData*) pdata_void;
	fName =  (char*) malloc(100);
	path = (char*) malloc(100);
	nx = pdata->nx;
	ny = pdata->ny;
	ux = pdata->uxCpy;
	uy = pdata->uyCpy;
	rho = pdata->rhoCpy;
	params = pdata->params;
	outputSelect = params->outputSelect;
	outDir = params->outDir;
	iter = pdata->iter;
	dt = params->dtPhys;
	dx = params->dxPhys;	
	rhoPhys = params->rhoPhys;
	c = dx/dt;
	
	if (outputSelect[0]) {//Write x-velocity
		strcpy(path,outDir);
		sprintf(fName,"/ux%d.csv", iter);
		strcat(path,fName);
		writeU(ux,nx,ny,c,path);
	}
	
	if (outputSelect[1]) {//Write y-velocity
		strcpy(path,outDir);
		sprintf(fName,"/uy%d.csv", iter);
		strcat(path,fName);
		writeU(uy,nx,ny,c,path);
	}
	
	if (outputSelect[2]) {//Write pressure
		strcpy(path,outDir);
		sprintf(fName,"/p%d.csv", iter);
		strcat(path,fName);
		writePres(rho,nx,ny,c,rhoPhys,path);
	}
	
	if (outputSelect[3]) {//Write vorticity
		strcpy(path,outDir);
		sprintf(fName,"/vort%d.csv", iter);
		strcat(path,fName);
		writeVorticity(ux,uy,nx,ny,dt,path);
	}
	
	if (outputSelect[4]) {//Write Force vector on obstacle
		strcpy(path,outDir);
		strcpy(fName, "/F.csv");
		strcat(path,fName);
		F = calcF(ny,params,rho,ux,uy);
		writeTimeSeries(F,2,path);
	}
	free(path);
	return NULL;
}

//------------------------------------------------------------------------------
//LBM    Function: printProgression
//------------------------------------------------------------------------------
//PURPOSE:	Prints percentage progression to console
//USAGE:	printProgression(iter,nIter)
//ARGUMENTS:
//			Name 	 	Type     	Description
//.............................................................................
//			iter     	int			Current iteration
//			nIter		int 		Total number of iterations
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void printProgression(int iter, int nIter) {
	double progression;
	progression = 100*iter/nIter;
	printf("%2.0f%%\b\b\b", progression);
	fflush(stdout);
}

//------------------------------------------------------------------------------
//LBM    Function: calcF
//------------------------------------------------------------------------------
//PURPOSE:	Calculate the total force of the obstacle. This is done by traversing
//			the boolean mask of the obstacle, and integrating the pressure of all
//			cells that are exposed to the fluid. The skin friction is calculated
//			by using a forward difference approximation of the the derivative of
//			the tangetial velocity, in normal directio. This value is mulitplied
//			with the dynamic viscosity, to obtain the shear stress, from Newtons
//			law of wet friction. This stress is integrated over the surface.
//USAGE:	F = calcF(ny,params,rho,ux,uy)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			ny		int  	   		Number of nodes in y direction
//			rho     double*			Density matrix
//			ux		double*			x-velocity field
//			uy		double*			y-velocity field
//.............................................................................
//RETURNS:
//			F		double*			2-vector containing forces [Fx, Fy]
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
double* calcF(int ny, SimParams* params, double* rho, double* ux, double* uy) {
	double p,*F,c,dx,dt,rhoPhys,visc,fric;
	int i,j,k,*bbCells,*bbCellMat,nBBcells;
	
	bbCells = params->bbCells;
	bbCellMat = params->bbCellMat;
	nBBcells = params->nBBcells;
	dx = params->dxPhys;
	dt = params->dtPhys;
	rhoPhys = params->rhoPhys;
	visc = params->nuPhys;
	c = dx/dt;
	F = calloc(2,sizeof(double));

	for (k = 0; k < nBBcells; k++) {
		i = bbCells[k];
		j = bbCells[1*nBBcells +k];
		
		if (bbCellMat[(i+1)*ny + j] == 0) {//Check if the node is exposed to the fluid to its right
			p = c*c*rhoPhys*(1.0/3.0)*(rho[(i+1)*ny + j]-1);//Calculate pressure
			fric = c*rhoPhys*visc*(uy[(i+2)*ny + j] - uy[(i+1)*ny + j]);//Calculate shear stress due to skin friction, by a first order approximation of Newtons law of Friction
			F[0] -= p*dx;//Add pressure contribution
			F[1] += fric;//Add skin friction contribution
		}
		
		//Analogously for left, top and bottom
		if (bbCellMat[(i-1)*ny + j] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[(i-1)*ny + j]-1);
			fric = c*rhoPhys*visc*(uy[(i-2)*ny + j] - uy[(i-1)*ny + j]);
			F[0] += p*dx;
			F[1] += fric;
		}
		
		if (bbCellMat[i*ny + j+1] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[i*ny + j+1]-1);
			fric = c*rhoPhys*visc*(ux[i*ny + j+2] - ux[i*ny + j+1]);
			F[1] -= p*dx;
			F[0] += fric;
		}
		
		if (bbCellMat[i*ny + j-1] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[i*ny + j-1]-1);
			fric = c*rhoPhys*visc*(ux[i*ny + j-2] - ux[i*ny + j-1]);
			F[1] += p*dx;
			F[0] += fric;
		}
	}
	return F;
}
