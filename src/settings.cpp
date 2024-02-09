#include "include/settings.h"

namespace settings
{
    double sim::T_min = 0.25;
    double sim::T_step = 0.25;
    int sim::MCS_sample = 100;
    int sim::MCS_therm = 50;
    int sim::MCS_total =10;

    bool model::isPeriodic = true;
    double model::field = 0;
    int model::nClusters = 1;
    int model::nSpins = 24*model::nClusters;

    int counter::swap_trials = 0;
}