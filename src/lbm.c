//Main program
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utilities.h"

//Indexing: ux(i,j) = ux(x = i*h,y=j*h);

int main(int argc, char** argv) {
	//Define grid
	int nx = 400;
	int ny = 100;
	//Define obstacle //Remove??
	int cylX  = nx/5+1;
	int cylY = ny/2 + 3;
	int cylR = ny/10+1;
	
	//Define flow parameters
	double uIn = 0.1;
	double Re = 1000;
	double nu = uIn*cylR/Re;
	double omega = 1/(3*nu+0.5);//Relaxation parameter of BGK collision operator
	
	//Define simulations parameters
	int nIter = 400000;
	int cyclesPerWrite = 5;
	
	
	//D2Q9 lattice constants
	double w[] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};//Weights
	int ex[] = {0,1,0,-1,0,1,-1,-1,1};//x component of lattice vectors 
	int ey[] = {0,0,1,0,-1,1,1,1,-1};//y component of lattice vectors
	int opp[] = {0,3,4,1,2,7,8,5,6};//Index of opposite vectors
	
	//BCs
	int nBBcells = 1177;
	int* bbCells = readObstacleData(nBBcells);//Array of indices to bounce back nodes
	int indInlet = 1;
	int indOutlet = nx;
	
	//IC: Initialze flow as Poiseuille flow
	double H = ny-2;
	double** ux = initUx(nx,ny,uIn);
	double** uy = initUy(nx,ny);
	//double*** fIn = initFlow(10,10,10,1.0);
	printVecD(ux[0],ny);
  
  return 0;
}




