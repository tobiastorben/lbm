#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lbm.h"
#include "inputParser.h"

void mapObstacleCells(LatticeConsts* lc, SimParams* params);
void initUx(LatticeConsts* lc, FlowData* flow, SimParams* params);
void initRho(LatticeConsts* lc, FlowData* flow);
void initFin(LatticeConsts* lc, FlowData* flow);
void readObstacle(LatticeConsts* lc, SimParams* params);
void nonDimensionalize(LatticeConsts* lc, SimParams* params);
void initialize(FlowData* flow, SimParams* params, LatticeConsts* lc, char* inPath);
void setLatticeConstants(LatticeConsts* lc);

#endif