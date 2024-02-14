#include "include/settings.h"

namespace settings
{
    double sim::T_min = 1.15;
    double sim::T_step = 0.1;
    int sim::MCS_sample = 100;
    int sim::MCS_swap = 10;
    int sim::MCS_therm = 500;
    int sim::MCS_total = 10000;
    int sim::MCS_decorr = 10;

    bool model::isPeriodic = true;
    double model::field = 0;
    int model::nClusters = 4;
    int model::nSpins = 24*model::nClusters;

    int counter::swap_trials = 0;
}