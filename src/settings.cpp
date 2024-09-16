#include "include/settings.h"

namespace settings
{
    double sim::T_min = 0.01;
    double sim::T_max = 3.0;
    double sim::ratio = sim::T_min/sim::T_max;
    int sim::MCS_sample = 1;
    int sim::MCS_swap = 5;
    int sim::MCS_therm = 10000;
    int sim::MCS_total = 1000000;
    int sim::MCS_decorr = 10;
    bool sim::storeCorrelations = false;
    bool sim::ladderUpdate = false;

    bool model::isPeriodic = true;
    double model::field = 0.0;

    int model::nClusters = 512;
    int model::neighIdx = 3;
    int model::nSpins = 24*model::nClusters;
    int model::distances = 13;
    

    int counter::swap_trials = 0;
}