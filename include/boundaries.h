#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "structs.h"
void westXY(FlowData* flow, LatticeConsts* lc, SimParams* params);
void westPY(FlowData* flow, LatticeConsts* lc, SimParams* params);
void eastPY(FlowData* flow, LatticeConsts* lc);
void eastXY(FlowData* flow, LatticeConsts* lc);
void southPX(FlowData* flow, LatticeConsts * lc, SimParams* params,int startX, int endX);
void southXY(FlowData* flow, LatticeConsts * lc, SimParams* params,int startX, int endX);
void northXY(FlowData* flow, LatticeConsts * lc, SimParams* params,int startX, int endX);
void northPX(FlowData* flow, LatticeConsts * lc, SimParams* params,int startX, int endX);
void bounce(FlowData* flow, LatticeConsts* lc, SimParams* params);
#endif