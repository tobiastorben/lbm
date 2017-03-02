#ifndef CORE_H
#define CORE_H

#include "lbm.h"
#include "boundaries.h"

void step(FlowData* flow, LatticeConsts* lc, SimParams* params);
void updateU(FlowData* flow, LatticeConsts* lc);
void updateRho(FlowData* flow, LatticeConsts* lc);
void collide(FlowData* flow, LatticeConsts* lc, SimParams* params);
void stream(FlowData* flow, LatticeConsts* lc);

#endif