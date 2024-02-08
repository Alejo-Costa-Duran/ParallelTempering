#ifndef MODEL
#define MODEL
#include "rngenerator.h"
#include <vector>

class model
{
    public:
        model();
        int nSpins;
        std::vector<int> lattice;
        double E, E_trial;
        double M, M_trial;
        double h = 0;
        std::vector<int> neighboursMatrix;
        
        void randomize_lattice();
        void flip();
        void accept_state();
        void set_E();
        void set_M();
};

#endif