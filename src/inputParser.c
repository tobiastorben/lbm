#include "inputParser.h"

int parseInput(char* inPath, SimParams* params) {
	FILE* fp = fopen(inPath,"r");
	char* line =  (char*) malloc(1000);
	params->obstaclePath = (char*) malloc(1000);
	params->outDir = (char*) malloc(1000);
	int* outputSelect =  (int*) calloc(5,sizeof(int));
	char* token;
	int count = 0;
	int outDirIsSet = 0;
	
	while (!feof(fp)){
		fgets(line,1000,fp);
		token = strtok(line, "=");
		if (token[0] == '%' || token[0] == '\r' || token[0] == '\n'){
			continue;
		}
		else if(!strcmp(token,"obstaclePath")){
		strcpy(params->obstaclePath,strtok(strtok(NULL, "="),(NULL, "\r")));		
			count++;
		}
				
		else if(!strcmp(token,"dx")){
			params->dxPhys = atof(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"dt")){
			params->dtPhys = atof(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"viscosity")){
			params->nuPhys = atof(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"density")){
			params->rhoPhys = atof(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"inletVel")){
			params->u0Phys= atof(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"nIter")){
			params->nIter = atoi(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"cyclesPerWrite")){
			params->cyclesPerWrite = atoi(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"startWrite")){
			params->startWrite = atoi(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"outDir")){
			strcpy(params->outDir,strtok(strtok(NULL, "="),(NULL, "\r")));
			outDirIsSet = 1;
		}
		
		else if(!strcmp(token,"nThreads")){
			params->nThreads = atoi(strtok(NULL, "="));	
			count++;
		}
				
		else if(!strcmp(token,"outputSelect")){
			token = strtok(NULL, ",\r\n");
			while( token != NULL ){			
				if (!strcmp(token,"ux")) {
					outputSelect[0] =1;
				}
				else if (!strcmp(token,"uy")) {
					outputSelect[1] = 1;
				}
				else if (!strcmp(token,"p")) {
					outputSelect[2] = 1;
				}
				else if (!strcmp(token,"vorticity")) {
					outputSelect[3] = 1;
				}
				else if (!strcmp(token,"F")) {
					outputSelect[4] = 1;
				}
				token = strtok(NULL, ",\r\n");				
		}
		}
		
	}
	if (!outDirIsSet) strcpy(params->outDir,"../res");
	params->outputSelect=outputSelect;
	free(line);
	fclose(fp);
	return !(count == 10);
}