#include "inputParser.h"

int parseInput(char* inPath, SimParams* params, BoundaryData* bcdata) {
	FILE* fp;
	char *line,*token;
	int *outputSelect,count,outDirIsSet;
	
	fp = fopen(inPath,"r");
	line =  (char*) malloc(1000);
	params->obstaclePath = (char*) malloc(1000);
	params->outDir = (char*) malloc(1000);
	outputSelect =  (int*) calloc(5,sizeof(int));
	count = 0;
	outDirIsSet = 0;
	
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
		
		else if(!strcmp(token,"uRef")){
			params->uRef = atof(strtok(NULL, "="));
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
		
		else if(!strcmp(token,"startVelX")){
			params->startVelX = atof(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"startVelY")){
			params->startVelY = atof(strtok(NULL, "="));
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
		
		else if(!strcmp(token,"west")){
			token = strtok(NULL, ",\r\n");
			if (!strcmp(token,"pres")) {
				bcdata->westFun = &westPY;
				bcdata->westBCType = 1;
			}
			else {
				bcdata->westFun = &westXY;
				bcdata->westBCType = 0;
			}
			bcdata->westBC[0] = atof(strtok(NULL, ",\r\n"));
			bcdata->westBC[1] = atof(strtok(NULL, ",\r\n"));	
			count++;
			}

			else if(!strcmp(token,"north")){
			token = strtok(NULL, ",\r\n");
			if (!strcmp(token,"pres")) {
				bcdata->northFun = &northPX;
				bcdata->northBCType = 1;
			}
			else {
				bcdata->northFun = &northXY;
				bcdata->northBCType = 0;
			}
			bcdata->northBC[0] = atof(strtok(NULL, ",\r\n"));
			bcdata->northBC[1] = atof(strtok(NULL, ",\r\n"));	
			count++;
			}
			
			else if(!strcmp(token,"east")){
			token = strtok(NULL, ",\r\n");
			if (!strcmp(token,"pres")) {
				bcdata->eastFun = &eastPY;
				bcdata->eastBCType = 1;
			}
			else {
				bcdata->eastFun = &eastXY;
				bcdata->eastBCType = 0;
			}
			bcdata->eastBC[0] = atof(strtok(NULL, ",\r\n"));
			bcdata->eastBC[1] = atof(strtok(NULL, ",\r\n"));	
			count++;
			}
			
			else if(!strcmp(token,"south")){
			token = strtok(NULL, ",\r\n");
			if (!strcmp(token,"pres")) {
				bcdata->southFun = &southPX;
				bcdata->southBCType = 1;
			}
			else {
				bcdata->southFun = &southXY;
				bcdata->southBCType = 0;
			}
			bcdata->southBC[0] = atof(strtok(NULL, ",\r\n"));
			bcdata->southBC[1] = atof(strtok(NULL, ",\r\n"));	
			count++;
			}
		}		
	if (!outDirIsSet) strcpy(params->outDir,"../res");
	params->outputSelect=outputSelect;
	free(line);
	fclose(fp);
	return !(count == 16);
}