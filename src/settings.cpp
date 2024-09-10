#include "include/settings.h"

namespace settings
{
    double sim::T_min = 0.01;
    double sim::T_max = 3.0;
    double sim::ratio = sim::T_min/sim::T_max;
    int sim::MCS_sample = 10;
    int sim::MCS_swap = 1;
    int sim::MCS_therm = 10000;
    int sim::MCS_total = 100000;
    int sim::MCS_decorr = 10;
    extern bool storeCorrelations = false;
    extern bool ladderUpdate = false;

    bool model::isPeriodic = true;
    double model::field = 0.0;

    int model::nClusters = 64;
    int model::neighIdx = 2;
    int model::nSpins = 24*model::nClusters;
    int model::distances = 9;
    

    int counter::swap_trials = 0;
}