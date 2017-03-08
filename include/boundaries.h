#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "structs.h"

void inlet(FlowData* flow, LatticeConsts* lc, SimParams* params);
void outlet(FlowData* flow, LatticeConsts* lc);
void bounce(FlowData* flow, LatticeConsts* lc, SimParams* params);

#endif