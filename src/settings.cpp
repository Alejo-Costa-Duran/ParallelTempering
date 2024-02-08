#include "include/settings.h"

namespace settings
{
    double sim::T_min = 0.05;
    double sim::T_max = 7.0;
    int sim::MCS_sample = 100;
    int sim::MCS_therm = 50;

    bool model::isPeriodic = true;
    double model::field = 0;
    int model::nClusters = 1;
    int model::nSpins = 24*model::nClusters;
}