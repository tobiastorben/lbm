#include "inputParser.h"

int parseInput(char* path, char* obstaclePath, double* Re_p, double* L_p, int* nIter_p, int* cyclesPerWrite_p, int* startWrite_p, int* outputSelect, char* outDir) {
	FILE* fp = fopen(path,"r");
	char* line = malloc(1000);
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
		strcpy(obstaclePath,strtok(strtok(NULL, "="),(NULL, "\r")));
			count++;
		}
				
		else if(!strcmp(token,"Re")){
			*Re_p = atof(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"L")){
			*L_p = atof(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"nIter")){
			*nIter_p = atoi(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"cyclesPerWrite")){
			*cyclesPerWrite_p = atoi(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"startWrite")){
			*startWrite_p = atoi(strtok(NULL, "="));
			count++;
		}
		
		else if(!strcmp(token,"outDir")){
			strcpy(outDir,strtok(strtok(NULL, "="),(NULL, "\r")));
			outDirIsSet = 1;
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
	if (!outDirIsSet) strcpy(outDir,"..\\res");
	free(line);
	return !(count >= 6);
}