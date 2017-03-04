#ifndef CORE_H
#define CORE_H

#include "lbm.h"
#include "boundaries.h"
#include <pthread.h>

void step(FlowData* flow, LatticeConsts* lc, SimParams* params, ThreadData* tdata);
void updateU(FlowData* flow, LatticeConsts* lc, int startX, int endX);
void updateRho(FlowData* flow, LatticeConsts* lc, int startX, int endX);
void collide(FlowData* flow, LatticeConsts* lc,SimParams* params, int startX, int endX);
void stream(FlowData* flow, LatticeConsts* lc);
void* streamAndUpdateBlock(void* tdata_void);
void* collideBlock(void* tdata_void);
void streamBlockInterior(FlowData* flow,LatticeConsts* lc,int startX, int endX);
void streamBlockBoundaries(FlowData* flow,LatticeConsts* lc,SimParams* params);
#endif