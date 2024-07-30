#include "include/settings.h"

namespace settings
{
    double sim::T_min = 0.05;
    double sim::T_max = 1.0;
    double sim::ratio = sim::T_min/sim::T_max;
    int sim::MCS_sample = 100;
    int sim::MCS_swap = 1;
    int sim::MCS_therm = 500;
    int sim::MCS_total = 300;
    int sim::MCS_decorr = 10;

    bool model::isPeriodic = true;
    double model::field = 0.8;
    int model::nClusters = 15;
    int model::nSpins = 24*model::nClusters;

    int counter::swap_trials = 0;
}