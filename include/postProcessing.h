#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include <stdlib.h>
#include "lbm.h"

double* calcF(FlowData* flow, LatticeConsts* lc, SimParams* params);

#endif