#ifndef SETTINGS_H
#define SETTINGS_H

namespace settings
{
    namespace sim 
    {
        extern double T_min;
        extern double T_max;
        extern double ratio;
        extern int MCS_therm;
        extern int MCS_sample;
        extern int MCS_swap;
        extern int MCS_total;
        extern int MCS_decorr;
        extern inline int mod(int x, int y)
        {
        if( x%y < 0)
        {
            return x%y+y;
        }
        else
        {
            return x%y;
        }
        }
    }
    namespace model
    {
        extern bool isPeriodic;
        extern double field;
        extern int nClusters;
        extern int nSpins;
        extern int neighIdx;
        extern int distances;
    }
    namespace counter
    {
        extern int swap_trials;
    }
}

#endif