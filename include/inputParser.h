#ifndef INPUTPARSER_H
#define INPUTPARSER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "structs.h"
#include "boundaries.h"

int parseInput(char* inPath, SimParams* params, BoundaryData* bcdata);

#endif