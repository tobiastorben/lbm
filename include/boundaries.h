#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "structs.h"
void westXY(FlowData* flow, LatticeConsts* lc, double uxBC, double uyBC);
void westPY(FlowData* flow, LatticeConsts* lc, double rhoBC, double uyBC);
void eastPY(FlowData* flow, LatticeConsts* lc, double rhoBC, double uyBC);
void eastXY(FlowData* flow, LatticeConsts* lc, double uxBC, double uyBC);
void southPX(FlowData* flow, LatticeConsts * lc,int startX, int endX, double rhoBC, double uxBC);
void southXY(FlowData* flow, LatticeConsts * lc,int startX, int endX, double uxBC, double uyBC);
void northXY(FlowData* flow, LatticeConsts * lc,int startX, int endX, double uxBC, double uyBC);
void northPX(FlowData* flow, LatticeConsts * lc,int startX, int endX, double rhoBC, double uxBC);
void bounce(FlowData* flow, LatticeConsts* lc, SimParams* params);

#endif