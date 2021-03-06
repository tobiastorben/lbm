#include "lbm.h"
//------------------------------------------------------------------------------
//LBM    Function: main   
//------------------------------------------------------------------------------
//PURPOSE:	Program entry point. Declares structs and variables, calls 
//			initialization, steps through all iterations and writes output.
//USAGE:	error = main(argc,argv)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			argc    int     		Command line argument count
//			argv    char**  		Arrary of cmd line arguments
//.............................................................................
//RETURNS:
//			Name 	 Type     		Description
//.............................................................................
//			error    int      		Error code. 0 indicates success
//.............................................................................
//CALLS:				
//			initialize			Initializes simulation
//			step				Progresses the simulation one time step
//			writeResults		Writes results to file
//			printProgression	Print progression to console
//
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
int main(int argc, char* argv[]) {	
	SimParams params;
	FlowData flow;
	LatticeConsts lc;
	ThreadData* tdata;
	PrintData pdata;
	BoundaryData bcdata;
	pthread_t printThread;
	struct timeval begin, end;
	int iter;
	double elapsed;
	char* inPath;
	
	//Parse input file
	if (argc == 2) inPath = argv[1];
	else {
		inPath = (char*) malloc(100);
		strcpy(inPath,"../input/input.in");
	}
	
	//Initialize simulation
	initialize(&flow,&params,&lc,&tdata,&pdata,&bcdata,inPath);
	//Time execution	
	gettimeofday(&begin, NULL);
	
	//Main loop. Marches flow in time
	printf("Solving flow:  ");
	for (iter = 1; iter <= params.nIter; iter++){	
		
		//Write results to file
		if (iter >= (params.startWrite) && iter % (params.cyclesPerWrite) == 0) {
		writeResults(&flow,&lc,iter,&printThread,&pdata);
		}
				
		//Print progression to console
		printProgression(iter,params.nIter);
		
		//Do one time step
		step(&flow,&lc,&params,tdata);		
	}
	pthread_join(printThread,NULL);
	gettimeofday(&end, NULL);
	elapsed = (end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0);
	printf("\n\nElapsed time: %.2f s\n\n", elapsed);
	return 0;
}