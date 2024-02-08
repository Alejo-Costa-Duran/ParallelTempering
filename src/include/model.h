#ifndef MODEL
#define MODEL
#include "rngenerator.h"
#include <vector>
#include "settings.h"

class model
{
    public:
        bool isPeriodic; 
        model();
        int nSpins;
        std::vector<int> lattice;
        double E, E_trial;
        int M, M_trial;
        double delE;
        double field = 0;
        std::vector< std::vector<int>> neighboursMatrix;
        
        void randomize_lattice();
        void trialMove(int trialSite);
        void acceptMove(int trialSite);
        double compute_E();
        int compute_M();
        void set_E();
        void set_M();
};

#endif