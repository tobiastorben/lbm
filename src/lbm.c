//Main program
#include "lbm.h"
#include "utilities.h"
#include "initialization.h"
#include "core.h"
#include <sys/time.h>


int main(int argc, char* argv[]) {	
	SimParams params;
	FlowData flow;
	LatticeConsts lc;
	ThreadData tdata[4];
	time_t t1,t2;
	int iter,nx,ny;
	char* inPath;
	char *fName,*path,*outDir;
	
	
	//Parse input file
	if (argc == 2) inPath = argv[1];
	else {
		inPath = (char*) malloc(100);
		strcpy(inPath,"..\\input\\input.in");
	}
	
	//Initialize simulation
	initialize(&flow,&params,&lc,tdata,inPath);
	
	fName =  (char*) malloc(100);
	path = (char*) malloc(100);
	outDir = params.outDir;
	nx = lc.nx;
	ny = lc.ny;

	//Time execution
	//t1 = clock();
	struct timeval begin, end;
	gettimeofday(&begin, NULL);

	//Main loop. Marches flow in time
	printf("Solving flow:  ");
	for (iter = 0; iter < params.nIter; iter++){		
		//Do one time step
		step(&flow,&lc,&params,tdata);
		
		/*///Write results to file
		if (iter >= (params.startWrite) && iter % (params.cyclesPerWrite) == 0) {//TODO : OPTIMIZE!!
		writeResults(&flow, &lc, &params, iter);
		}
		
		strcpy(path,outDir);
		sprintf(fName,"\\rho%d.csv", iter);
		strcat(path,fName);
		csvWriteD(flow.rho,nx,ny,path);
		
		strcpy(path,outDir);
		sprintf(fName,"\\fIn%d.csv", iter);
		strcat(path,fName);
		csvWriteLayer(flow.fIn,nx,ny,4,path);
		
		strcpy(path,outDir);
		sprintf(fName,"\\fOut%d.csv", iter);
		strcat(path,fName);
		csvWriteLayer(flow.fOut,nx,ny,4,path);
				
		//Print progression to console
		printProgression(iter,params.nIter);*/
	}
	//t2 = clock();
	gettimeofday(&end, NULL);
	//get the total number of ms that the code took:
	double elapsed = (end.tv_sec - begin.tv_sec) + 
              ((end.tv_usec - begin.tv_usec)/1000000.0);
	printf("100%%\n\nElapsed time: %.4f s\n\n", elapsed);
	writeResults(&flow, &lc, &params, iter);//TODO : Set outputSelect to ones
	printf("Results written to: %s\n", params.outDir);
	return 0;
}



