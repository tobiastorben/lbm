//Main program
#include "lbm.h"
#include "utilities.h"
#include "initialization.h"
#include "boundaries.h"
#include "core.h"
#include "postProcessing.h"
#include "inputParser.h"

int main(int argc, char* argv[]) {	
	SimParams params;
	FlowData flow;
	LatticeConsts lc;
	time_t t1,t2;
	double progression;
	int iter;
	char* inPath;

	if (argc == 2) inPath = argv[1];
	else {
		inPath = (char*) malloc(100);
		strcpy(inPath,"..\\input\\input.in");
	}
	
	initialize(&flow,&params,&lc,inPath);

	//Time execution
	t1 = clock();
	
	//Main loop. Marches flow in time
	printf("Solving flow:  ");
	for (iter = 0; iter < params.nIter; iter++){
		//Update macroscopic variables	
		updateRho(&flow,&lc);		
		updateU(&flow,&lc);
		
		//Apply BCs
		inlet(&flow,&lc,&params);//Poiseuille (Zou/He)		
		outlet(&flow,&lc);//Constant pressure (Zou/He)
		
		//Collide (Bhatnagar-Gross-Kroot model)
		collide(&flow,&lc,&params);//Particle-Particle collisions		
		bounce(&flow,&lc,&params);//Bounce-back collision with boundary
		
		//Streaming
		stream(&flow,&lc);
		
		//Write results and progression
		if (iter >= (params.startWrite) && iter % (params.cyclesPerWrite) == 0) {
			writeResults(&flow, &lc, &params, iter);
		}
		progression = 100*iter/(params.nIter);
		printf("%2.0f%%\b\b\b", progression);
		fflush(stdout);
	}
	t2 = clock();
	printf("100%%\n\nElapsed time: %f.2 s\n\n", ((float) t2-(float) t1)/1000.0);
	writeResults(&flow, &lc, &params, iter);//TODO : Set outputSelect to ones
	printf("Results written to: %s\n", params.outDir);
	return 0;
}



