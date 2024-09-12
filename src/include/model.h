#ifndef MODEL_H
#define MODEL_H
#include "rngenerator.h"
#include <vector>
#include "settings.h"
#include <iostream>
#include <queue>
#include <limits> // For infinity
#include <climits> // For infinity


class model
{
    public:
        bool isPeriodic;
        model(); 
        model(int rank,int* shared_neighbours);
        int nSpins;
        int rank;
        int* lattice;
        double E, E_trial;
        int M, M_trial;
        double delE;
        double field;

        void randomize_lattice();
        void trialMove(int trialSite,int *shared_neighbours);
        void acceptMove(int trialSite);
        double compute_E(int *shared_neighbours);
        int compute_M();
        void set_E(int *shared_neighbours);
        void set_M();
        std::vector<double> compute_Correlations(int *shared_distances);

};

#endif