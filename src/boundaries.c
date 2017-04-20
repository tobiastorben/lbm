#include "boundaries.h"

//This file contains functions to apply the boundary conditions, as specified in the
//input file. It also contain the function to collide the flow with the obstacle.

//------------------------------------------------------------------------------
//LBM    	Function: southXY
//------------------------------------------------------------------------------
//PURPOSE:	Applies Zou/He velocity boundary conditions to south boundary. Assigns
//			the prescribed velocities, and calculates the density and the unknown
//			distributions using the BCs and the known distributions from the
//			streaming step.			
//USAGE:	southXY(flow,lc,startX,endX,uxBC,uyBC)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			startX	int				Index of first block column
//			endX	int				Index of last block column
//			uxBC	double			The value of x velocity at boundary
//			uyBC	double			The value of y velocity at boundary
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void southXY(FlowData* flow, LatticeConsts * lc,int startX, int endX, double uxBC, double uyBC) {
	double *fIn,*rho;
	int nx,ny,nxny;

	fIn = flow->fIn;
	rho = flow->rho;
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	
	for (int i = startX; i <= endX;i++) {
		rho[ny*i] = (1.0/(1.0-uyBC))*(fIn[i*ny] + fIn[nxny + i*ny] + fIn[3*nxny + i*ny]
					+ 2.0*(fIn[4*nxny + i*ny] + fIn[7*nxny + i*ny] + fIn[8*nxny + i*ny]));
		fIn[2*nxny + i*ny] = fIn[4*nxny + i*ny] + (2.0/3.0)*(rho[i*ny])*uyBC;
		fIn[5*nxny + i*ny] = fIn[7*nxny + i*ny] - 0.5*(fIn[nxny + i*ny]-fIn[3*nxny + i*ny]) + 0.5*(rho[i*ny])*uxBC + (1.0/6.0)*(rho[i*ny])*uyBC;
		fIn[6*nxny + i*ny] = fIn[8*nxny + i*ny] + 0.5*(fIn[nxny + i*ny]-fIn[3*nxny + i*ny]) - 0.5*(rho[i*ny])*uxBC + (1.0/6.0)*(rho[i*ny])*uyBC;
		flow->ux[ny*i] = uxBC;
		flow->uy[ny*i] = uyBC;		
	}
}

//------------------------------------------------------------------------------
//LBM    	Function: southPX
//------------------------------------------------------------------------------
//PURPOSE:	Applies Zou/He velocity/pressure BCs to south boundary. Assigns the
//			prescribed velocity/pressure, and calculates the unknown normal
//			velocity, and the unknown distributions using the BCs, and the known
//			distributions from the streaming step.			
//USAGE:	southPX(flow,lc,startX,endX,rhoBC,uxBC)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			startX	int				Index of first block column
//			endX	int				Index of last block column
//			rhoBC	double			The value the density at the boundary
//			uxBC	double			The value of x velocity at boundary
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void southPX(FlowData* flow, LatticeConsts * lc,int startX, int endX, double rhoBC, double uxBC) {
	double *fIn,*uy;
	int nx,ny,nxny;
	
	fIn = flow->fIn;
	uy = flow->uy;
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	
	for (int i = startX; i <= endX;i++) {
		flow->ux[ny*i] = uxBC;
		flow->rho[ny*i] = rhoBC;
		uy[i*ny] = 1 - (fIn[i*ny] + fIn[nxny + i*ny] + fIn[3*nxny + i*ny]
					+ 2.0*(fIn[4*nxny + i*ny] + fIn[7*nxny + i*ny] + fIn[8*nxny + i*ny]))/rhoBC;
		fIn[2*nxny + i*ny] = fIn[4*nxny + i*ny] + (2.0/3.0)*(rhoBC)*uy[ny*i];
		fIn[5*nxny + i*ny] = fIn[7*nxny + i*ny] - 0.5*(fIn[nxny + i*ny]-fIn[3*nxny + i*ny]) + 0.5*(rhoBC)*uxBC + (1.0/6.0)*(rhoBC)*uy[ny*i];
		fIn[6*nxny + i*ny] = fIn[8*nxny + i*ny] + 0.5*(fIn[nxny + i*ny]-fIn[3*nxny + i*ny]) - 0.5*(rhoBC)*uxBC + (1.0/6.0)*(rhoBC)*uy[ny*i];
		flow->ux[ny*i] = uxBC;
	}
}

//------------------------------------------------------------------------------
//LBM    	Function: northXY
//------------------------------------------------------------------------------
//PURPOSE:	Applies Zou/He velocity boundary conditions to north boundary. Assigns
//			the prescribed velocities, and calculates the density and the unknown
//			distributions using the BCs and the known distributions from the
//			streaming step.			
//USAGE:	northXY(flow,lc,startX,endX,uxBC,uyBC)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			startX	int				Index of first block column
//			endX	int				Index of last block column
//			uxBC	double			The value of x velocity at boundary
//			uyBC	double			The value of y velocity at boundary
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void northXY(FlowData* flow, LatticeConsts* lc,int startX, int endX, double uxBC, double uyBC) {
	double *fIn,*rho;
	int nx,ny,nxny;
	
	fIn = flow->fIn;
	rho = flow->rho;
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	
	for (int i = startX; i <= endX;i++) {
		rho[ny*i + ny -1] = (1.0/(1.0+uyBC))*(fIn[ny*i + ny -1] + fIn[nxny + ny*i + ny -1] + fIn[3*nxny + ny*i + ny -1]
					+ 2.0*(fIn[2*nxny + ny*i + ny -1] + fIn[5*nxny + ny*i + ny -1] + fIn[6*nxny + ny*i + ny -1]));
		fIn[4*nxny + ny*i + ny -1] = fIn[2*nxny + ny*i + ny -1] - (2.0/3.0)*(rho[ny*i + ny -1])*uyBC;
		fIn[8*nxny + ny*i + ny -1] = fIn[6*nxny + ny*i + ny -1] + 0.5*(fIn[3*nxny + ny*i + ny -1]-fIn[nxny + ny*i + ny -1]) + 0.5*(rho[ny*i + ny -1])*uxBC - (1.0/6.0)*(rho[ny*i + ny -1])*uyBC;
		fIn[7*nxny + ny*i + ny -1] = fIn[5*nxny + ny*i + ny -1] + 0.5*(fIn[nxny + ny*i + ny -1]-fIn[3*nxny + ny*i + ny -1]) - 0.5*(rho[ny*i + ny -1])*uxBC - (1.0/6.0)*(rho[ny*i + ny -1])*uyBC;
		flow->ux[ny*i + ny-1] = uxBC;
		flow->uy[ny*i + ny-1] = uyBC;		
	}
}

//------------------------------------------------------------------------------
//LBM    	Function: northPX
//------------------------------------------------------------------------------
//PURPOSE:	Applies Zou/He velocity/pressure BCs to north boundary. Assigns the
//			prescribed velocity/pressure, and calculates the unknown normal
//			velocity, and the unknown distributions using the BCs, and the known
//			distributions from the streaming step.			
//USAGE:	northPX(flow,lc,startX,endX,rhoBC,uxBC)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			startX	int				Index of first block column
//			endX	int				Index of last block column
//			rhoBC	double			The value the density at the boundary
//			uxBC	double			The value of x velocity at boundary
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void northPX(FlowData* flow, LatticeConsts * lc,int startX, int endX, double rhoBC, double uxBC) {
	double *fIn,*uy;
	int nx,ny,nxny;

	fIn = flow->fIn;
	uy = flow->uy;
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	
	for (int i = startX; i <= endX;i++) {
		flow->ux[ny*i + ny-1] = uxBC;
		flow->rho[ny*i + ny-1] = rhoBC;
		uy[ny*i + ny-1] = -1.0 + (fIn[ny*i + ny -1] + fIn[nxny + ny*i + ny -1] + fIn[3*nxny + ny*i + ny -1]
					+ 2.0*(fIn[2*nxny + ny*i + ny -1] + fIn[5*nxny + ny*i + ny -1] + fIn[6*nxny + ny*i + ny -1]))/rhoBC;					
		fIn[4*nxny + ny*i + ny -1] = fIn[2*nxny + ny*i + ny -1] - (2.0/3.0)*(rhoBC)*uy[ny*i + ny-1];
		fIn[8*nxny + ny*i + ny -1] = fIn[6*nxny + ny*i + ny -1] + 0.5*(fIn[3*nxny + ny*i + ny -1]-fIn[nxny + ny*i + ny -1]) + 0.5*(rhoBC)*uxBC - (1.0/6.0)*(rhoBC)*uy[ny*i + ny-1];
		fIn[7*nxny + ny*i + ny -1] = fIn[5*nxny + ny*i + ny -1] + 0.5*(fIn[nxny + ny*i + ny -1]-fIn[3*nxny + ny*i + ny -1]) - 0.5*(rhoBC)*uxBC - (1.0/6.0)*(rhoBC)*uy[ny*i + ny-1];		
	}
}
//------------------------------------------------------------------------------
//LBM    	Function: westXY
//------------------------------------------------------------------------------
//PURPOSE:	Applies Zou/He velocity boundary conditions to west boundary. Assigns
//			the prescribed velocities, and calculates the density and the unknown
//			distributions using the BCs and the known distributions from the
//			streaming step. The corner nodes are treaded by extrapolation of the
//			neighbour in x-direction.	
//USAGE:	westXY(flow,lc,uxBC,uyBC)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			uxBC	double			The value of x velocity at boundary
//			uyBC	double			The value of y velocity at boundary
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void westXY(FlowData* flow, LatticeConsts* lc, double uxBC, double uyBC) {
	double sum1, sum2, rhoJ,*fIn,*rho;
	int nx,ny,k;

	nx = lc->nx;
	ny = lc->ny;	
	fIn = flow->fIn;
	rho = flow->rho;

	for (int j = 1; j<ny-1;j++) {
		flow->ux[j] = uxBC;
		flow->uy[j] = uyBC;
		sum1 = fIn[j] + fIn[2*nx*ny +j] + fIn[4*nx*ny+j];
		sum2 = fIn[3*nx*ny +j] + fIn[6*nx*ny + j] + fIn[7*nx*ny + j];
		rhoJ = (1.0/(1.0-uxBC))*(sum1 + 2.0*sum2);
		fIn[1*nx*ny+j] = fIn[3*nx*ny+j] + (2.0/3.0)*rhoJ*uxBC;
		fIn[5*nx*ny+j] = fIn[7*nx*ny+j] + 0.5*(fIn[4*nx*ny+j]-fIn[2*nx*ny+j]) + 0.5*(rhoJ*uyBC) + (1.0/6)*(rhoJ*uxBC);
		fIn[8*nx*ny+j] = fIn[6*nx*ny+j] + 0.5*(fIn[2*nx*ny+j]-fIn[4*nx*ny+j]) - 0.5*(rhoJ*uyBC) + (1.0/6)*(rhoJ*uxBC);
		flow->rho[j] = rhoJ;
	}
	//Bottom corner
	rho[0*ny + 0] = rho[1*ny + 0];
	flow->ux[0*ny + 0] = flow->ux[1*ny + 0];
	flow->uy[0*ny + 0] = flow->uy[1*ny + 0];
	for (k = 0; k < 9; k++){
	fIn[k*nx*ny + 0*ny + 0] = fIn[k*nx*ny + 1*ny + 0];
	}
	
	//Top corner
	rho[0*ny + ny-1] = rho[1*ny + ny-1];
	flow->ux[0*ny + ny-1] = flow->ux[1*ny + ny-1];
	flow->uy[0*ny + ny-1] = flow->uy[1*ny + ny-1];
	for (k = 0; k < 9; k++){
	fIn[k*nx*ny + 0*ny + ny-1] = fIn[k*nx*ny + 1*ny + ny-1];
	}	
}

//------------------------------------------------------------------------------
//LBM    	Function: westPY
//------------------------------------------------------------------------------
//PURPOSE:	Applies Zou/He velocity/pressure BCs to west boundary. Assigns the
//			prescribed velocity/pressure, and calculates the unknown normal
//			velocity, and the unknown distributions using the BCs, and the known
//			distributions from the streaming step. The corner nodes are treaded
//			 by extrapolation of the neighbour in x-direction.	
//USAGE:	westPY(flow,lc,rhoBC,uyBC)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			rhoBC	double			The value the density at the boundary
//			uyBC	double			The value of x velocity at boundary
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void westPY(FlowData* flow, LatticeConsts* lc, double rhoBC, double uyBC) {
	double sum1,sum2,*fIn,*ux,*rho;
	int nx,ny,k;

	nx = lc->nx;
	ny = lc->ny;	
	fIn = flow->fIn;
	ux = flow->ux;
	rho = flow->rho;

	for (int j = 1; j<ny-1; j++) {
		flow->rho[j] = rhoBC;
		flow->uy[j] = uyBC;
		sum1 = fIn[j] + fIn[2*nx*ny +j] + fIn[4*nx*ny+j];
		sum2 = fIn[3*nx*ny +j] + fIn[6*nx*ny + j] + fIn[7*nx*ny + j];
		ux[j] = 1.0-(sum1 + 2.0*sum2)/rhoBC;
		fIn[1*nx*ny+j] = fIn[3*nx*ny+j] + (2.0/3.0)*rhoBC*ux[j];
		fIn[5*nx*ny+j] = fIn[7*nx*ny+j] + 0.5*(fIn[4*nx*ny+j]-fIn[2*nx*ny+j]) + 0.5*(rhoBC*uyBC) + (1.0/6)*(rhoBC*ux[j]);
		fIn[8*nx*ny+j] = fIn[6*nx*ny+j] + 0.5*(fIn[2*nx*ny+j]-fIn[4*nx*ny+j]) - 0.5*(rhoBC*uyBC) + (1.0/6)*(rhoBC*ux[j]);
	}
	
	//Bottom corner
	rho[0*ny + 0] = rho[1*ny + 0];
	flow->ux[0*ny + 0] = flow->ux[1*ny + 0];
	flow->uy[0*ny + 0] = flow->uy[1*ny + 0];
	for (k = 0; k < 9; k++){
	fIn[k*nx*ny + 0*ny + 0] = fIn[k*nx*ny + 1*ny + 0];
	}
	
	//Top corner
	rho[0*ny + ny-1] = rho[1*ny + ny-1];
	flow->ux[0*ny + ny-1] = flow->ux[1*ny + ny-1];
	flow->uy[0*ny + ny-1] = flow->uy[1*ny + ny-1];
	for (k = 0; k < 9; k++){
	fIn[k*nx*ny + 0*ny + ny-1] = fIn[k*nx*ny + 1*ny + ny-1];
	}	
}

//------------------------------------------------------------------------------
//LBM    	Function: eastPY
//------------------------------------------------------------------------------
//PURPOSE:	Applies Zou/He velocity/pressure BCs to east boundary. Assigns the
//			prescribed velocity/pressure, and calculates the unknown normal
//			velocity, and the unknown distributions using the BCs, and the known
//			distributions from the streaming step. The corner nodes are treaded
//			 by extrapolation of the neighbour in x-direction.	
//USAGE:	eastPY(flow,lc,rhoBC,uyBC)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			rhoBC	double			The value the density at the boundary
//			uyBC	double			The value of x velocity at boundary
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void eastPY(FlowData* flow, LatticeConsts* lc, double rhoBC, double uyBC) {
	double sum1,sum2,*fIn,*ux,*rho;
	int nx,ny,k;
	
	nx = lc->nx;
	ny = lc->ny;	
	fIn = flow->fIn;
	ux = flow->ux;
	rho = flow->rho;
	
	for (int j = 1; j<ny-1; j++) {
		flow->rho[(nx-1)*ny +j] = rhoBC;
		flow->uy[(nx-1)*ny +j] = uyBC;
		sum1 = fIn[(nx-1)*ny+j] + fIn[2*nx*ny+(nx-1)*ny+j] + fIn[4*nx*ny+(nx-1)*ny+j];
		sum2 = fIn[1*nx*ny+(nx-1)*ny+j] + fIn[5*nx*ny+(nx-1)*ny+j] + fIn[8*nx*ny+(nx-1)*ny+j];
		ux[(nx-1)*ny +j] = -1.0 + (sum1 + 2.0*sum2)/rhoBC;
		fIn[3*nx*ny+(nx-1)*ny + j] = fIn[1*nx*ny+(nx-1)*ny + j] - (2.0/3)*ux[(nx-1)*ny +j];
		fIn[7*nx*ny+(nx-1)*ny + j] = fIn[5*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[2*nx*ny+(nx-1)*ny + j]-fIn[4*nx*ny+(nx-1)*ny + j]) - 0.5*rhoBC*uyBC - (1.0/6)*rhoBC*ux[(nx-1)*ny +j];
		fIn[6*nx*ny+(nx-1)*ny + j] = fIn[8*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[4*nx*ny+(nx-1)*ny + j]-fIn[2*nx*ny+(nx-1)*ny + j]) + 0.5*rhoBC*uyBC - (1.0/6)*rhoBC*ux[(nx-1)*ny +j];
	}
	//Bottom corner
	rho[(nx-1)*ny + 0] = rho[(nx-2)*ny + 0];
	flow->ux[(nx-1)*ny + 0] = flow->ux[(nx-2)*ny + 0];
	flow->uy[(nx-1)*ny + 0] = flow->uy[(nx-2)*ny + 0];
	for (k = 0; k < 9; k++){
	fIn[k*nx*ny + (nx-1)*ny + 0] = fIn[k*nx*ny + (nx-2)*ny + 0];
	}
	
	//Top corner
	rho[(nx-1)*ny + ny-1] = rho[(nx-2)*ny + ny-1];
	flow->ux[(nx-1)*ny + ny-1] = flow->ux[(nx-2)*ny + ny-1];
	flow->uy[(nx-1)*ny + ny-1] = flow->uy[(nx-2)*ny + ny-1];
	for (k = 0; k < 9; k++){
	fIn[k*nx*ny + (nx-1)*ny + ny-1] = fIn[k*nx*ny + (nx-2)*ny + ny-1];
	}
}
//------------------------------------------------------------------------------
//LBM    	Function: eastXY
//------------------------------------------------------------------------------
//PURPOSE:	Applies Zou/He velocity boundary conditions to west boundary. Assigns
//			the prescribed velocities, and calculates the density and the unknown
//			distributions using the BCs and the known distributions from the
//			streaming step. The corner nodes are treaded by extrapolation of the
//			neighbour in x-direction.	
//USAGE:	eastflow,lc,uxBC,uyBC)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			uxBC	double			The value of x velocity at boundary
//			uyBC	double			The value of y velocity at boundary
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void eastXY(FlowData* flow, LatticeConsts* lc, double uxBC, double uyBC) {
	double sum1,sum2,*fIn,*rho,rhoJ;
	int nx,ny,k;
	
	nx = lc->nx;
	ny = lc->ny;
	fIn = flow->fIn;
	rho = flow->rho;
	
	for (int j = 1; j<ny-1;j++) {
		flow->ux[(nx-1)*ny +j] = uxBC;
		flow->uy[(nx-1)*ny +j] = uyBC;
		sum1 = fIn[(nx-1)*ny+j] + fIn[2*nx*ny+(nx-1)*ny+j] + fIn[4*nx*ny+(nx-1)*ny+j];
		sum2 = fIn[1*nx*ny+(nx-1)*ny+j] + fIn[5*nx*ny+(nx-1)*ny+j] + fIn[8*nx*ny+(nx-1)*ny+j];
		rhoJ = (sum1 + 2.0*sum2)/(1.0+uxBC);
		fIn[3*nx*ny+(nx-1)*ny + j] = fIn[1*nx*ny+(nx-1)*ny + j] - (2.0/3)*rhoJ*uxBC;
		fIn[7*nx*ny+(nx-1)*ny + j] = fIn[5*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[2*nx*ny+(nx-1)*ny + j]-fIn[4*nx*ny+(nx-1)*ny + j]) - 0.5*rhoJ*uyBC - (1.0/6)*rhoJ*uxBC;
		fIn[6*nx*ny+(nx-1)*ny + j] = fIn[8*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[4*nx*ny+(nx-1)*ny + j]-fIn[2*nx*ny+(nx-1)*ny + j]) + 0.5*rhoJ*uyBC - (1.0/6)*rhoJ*uxBC;
		rho[(nx-1)*ny +j] = rhoJ;
	}
	
	//Bottom corner
	rho[(nx-1)*ny + 0] = rho[(nx-2)*ny + 0];
	flow->ux[(nx-1)*ny + 0] = flow->ux[(nx-2)*ny + 0];
	flow->uy[(nx-1)*ny + 0] = flow->uy[(nx-2)*ny + 0];
	for (k = 0; k < 9; k++){
	fIn[k*nx*ny + (nx-1)*ny + 0] = fIn[k*nx*ny + (nx-2)*ny + 0];
	}
	
	//Top corner
	rho[(nx-1)*ny + ny-1] = rho[(nx-2)*ny + ny-1];
	flow->ux[(nx-1)*ny + ny-1] = flow->ux[(nx-2)*ny + ny-1];
	flow->uy[(nx-1)*ny + ny-1] = flow->uy[(nx-2)*ny + ny-1];
	for (k = 0; k < 9; k++){
	fIn[k*nx*ny + (nx-1)*ny + ny-1] = fIn[k*nx*ny + (nx-2)*ny + ny-1];
	}
}

//------------------------------------------------------------------------------
//LBM    	Function: bounce
//------------------------------------------------------------------------------
//PURPOSE:	Apply the collision to the obstacle. This is done by on-wall bounce
//			back. This reduces to assigning the outgoing distribution in one
//			direction, to the incoming distribution in the opposite direction.
//			This process is done for all '1' nodes in the obstacle mask. This
//			function is not executed in parallell, but demands very little 
//			CPU time.
//USAGE:	bounce(flow,lc,params)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			params	SimParams*		The simulation parameters
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void bounce(FlowData* flow, LatticeConsts* lc, SimParams* params) {
	int i,j,k,l,nx,ny,*opp;
	
	nx = lc->nx;
	ny = lc->ny;
	opp = lc->opposite;
	
	for (l = 0; l < params->nBBcells; l++) {
		for (k = 0; k < 9; k++){
				i = params->bbCells[l];
				j = params->bbCells[params->nBBcells+l];
				flow->fOut[nx*ny*k + ny*i + j] = flow->fIn[nx*ny*opp[k] + ny*i + j];
		}
	}	
}