//Main program
#include "lbm.h"
#include "utilities.h"
#include "initialization.h"
#include "core.h"


int main(int argc, char* argv[]) {	
	SimParams params;
	FlowData flow;
	LatticeConsts lc;
	time_t t1,t2;
	double progression;
	int iter;
	char* inPath;
	
	//Parse input file
	if (argc == 2) inPath = argv[1];
	else {
		inPath = (char*) malloc(100);
		strcpy(inPath,"..\\input\\input.in");
	}
	
	//Initialize simulation
	initialize(&flow,&params,&lc,inPath);

	//Time execution
	t1 = clock();
	
	//Main loop. Marches flow in time
	printf("Solving flow:  ");
	for (iter = 0; iter < params.nIter; iter++){		
		//Do one time step
		step(&flow,&lc,&params);
		
		/*//Write results to file
		if (iter >= (params.startWrite) && iter % (params.cyclesPerWrite) == 0) {//TODO : OPTIMIZE!!
		writeResults(&flow, &lc, &params, iter);
		}
				
		//Print progression to console
		printProgression(iter,params.nIter);*/
	}
	t2 = clock();
	printf("100%%\n\nElapsed time: %f.2 s\n\n", ((float) t2-(float) t1)/1000.0);
	writeResults(&flow, &lc, &params, iter);//TODO : Set outputSelect to ones
	printf("Results written to: %s\n", params.outDir);
	return 0;
}



