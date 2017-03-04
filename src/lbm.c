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
	ThreadData* tdata;
	struct timeval begin, end;
	int iter;
	char* inPath;
		
	//Parse input file
	if (argc == 2) inPath = argv[1];
	else {
		inPath = (char*) malloc(100);
		strcpy(inPath,"..\\input\\input.in");
	}
	
	//Initialize simulation
	initialize(&flow,&params,&lc,&tdata,inPath);

	//Time execution	
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
				
		//Print progression to console
		printProgression(iter,params.nIter);*/
	}
	gettimeofday(&end, NULL);
	double elapsed = (end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0);
	printf("100%%\n\nElapsed time: %.2f s\n\n", elapsed);
	writeResults(&flow, &lc, &params, iter);//TODO : Set outputSelect to ones
	printf("Results written to: %s\n", params.outDir);
	return 0;
}



