#ifndef CORE_H
#define CORE_H

#include "structs.h"
#include "boundaries.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void step(FlowData* flow, LatticeConsts* lc, SimParams* params, ThreadData* tdata);
void updateU(FlowData* flow, LatticeConsts* lc, int startX, int endX);
void updateRho(FlowData* flow, LatticeConsts* lc, int startX, int endX);
void collide(FlowData* flow, LatticeConsts* lc,SimParams* params, int startX, int endX);
void* updateBlock(void* tdata_void);
void streamBlockInterior(FlowData* flow,LatticeConsts* lc,int startX, int endX);
void streamFirstBlockInterior(FlowData* flow,LatticeConsts* lc,int endX);
void streamLastBlockInterior(FlowData* flow,LatticeConsts* lc,int startX);
void streamBlockBoundaries(FlowData* flow,LatticeConsts* lc,SimParams* params);
void* updateFirstBlock(void* tdata_void);
void* updateLastBlock(void* tdata_void);
void stream(FlowData* flow, LatticeConsts* lc);

#endif