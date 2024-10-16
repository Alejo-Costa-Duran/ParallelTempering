#include "include/settings.h"

namespace settings
{
    double sim::T_min = 0.01;
    double sim::T_max = 4.0;
    double sim::ratio = sim::T_min/sim::T_max;
    int sim::MCS_sample = 1;
    int sim::MCS_swap = 25;
    int sim::MCS_therm = 100;
    int sim::MCS_total = 100;
    int sim::MCS_decorr = 10;
    bool sim::storeCorrelations = false;
    bool sim::ladderUpdate = true;

    bool model::isPeriodic = true;
    double model::field = 0;

    int model::nClusters = 64;
    int model::neighIdx = 2;
    int model::nSpins = 24*model::nClusters;
    int model::distances = 5;
    

    int counter::swap_trials = 0;
}