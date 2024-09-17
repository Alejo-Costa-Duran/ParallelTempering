#include "include/settings.h"

namespace settings
{
    double sim::T_min = 0.01;
    double sim::T_max = 3.0;
    double sim::ratio = sim::T_min/sim::T_max;
    int sim::MCS_sample = 1;
    int sim::MCS_swap = 1;
    int sim::MCS_therm = 10000;
    int sim::MCS_total = 1000000;
    int sim::MCS_decorr = 10;
    bool sim::storeCorrelations = false;
    bool sim::ladderUpdate = true;

    bool model::isPeriodic = true;
    double model::field = 0.5;

    int model::nClusters = 8;
    int model::neighIdx = 1;
    int model::nSpins = 24*model::nClusters;
    int model::distances = 5;
    

    int counter::swap_trials = 0;
}