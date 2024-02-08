#ifndef SETTINGS_H
#define SETTINGS_H

namespace settings
{
    namespace sim 
    {
        extern double T_min;
        extern double T_max;
        extern int MCS_therm;
        extern int MCS_sample;
    }
    namespace model
    {
        extern bool isPeriodic;
        extern double field;
        extern int nClusters;
        extern int nSpins;
    }
}

#endif