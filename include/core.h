#ifndef CORE_H
#define CORE_H

#include "structs.h"
#include "boundaries.h"

void step(FlowData* flow, LatticeConsts* lc, SimParams* params, ThreadData* tdata);
void updateU(FlowData* flow, LatticeConsts* lc, int startX, int endX);
void updateRho(FlowData* flow, LatticeConsts* lc, int startX, int endX);
void collide(FlowData* flow, LatticeConsts* lc,SimParams* params, int startX, int endX);
void* updateBlock(void* tdata_void);
void streamBlockInterior(FlowData* flow,LatticeConsts* lc,int startX, int endX);
void streamBlockBoundaries(FlowData* flow,LatticeConsts* lc,SimParams* params);
void* updateFirstBlock(void* tdata_void);
void* updateLastBlock(void* tdata_void);

#endif