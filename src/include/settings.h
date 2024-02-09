#ifndef SETTINGS_H
#define SETTINGS_H

namespace settings
{
    namespace sim 
    {
        extern double T_min;
        extern double T_step;
        extern int MCS_therm;
        extern int MCS_sample;
        extern int MCS_total;
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
    }
    namespace counter
    {
        extern int swap_trials;
    }
}

#endif