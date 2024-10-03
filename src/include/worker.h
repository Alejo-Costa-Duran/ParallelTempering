#ifndef WORKER_H
#define WORKER_H
#include "model.h"
#include <mpi.h>
#include <math.h>
#include "counters.h"
#include "settings.h"



class worker
{
    private:
        std::vector<std::vector<double>> probabilities_list;
        
        inline void compute_probabilities();
        void cooldown(int *shared_neighbours);
        void start_counters();
        inline bool performTrialMove(double temp,int trialSite, int *shared_neighbours);
        void storeThermodynamicData(bool storeCorrelations, int *distanceMatrix);
        void storeCounters();

    public:
        worker();
        worker(int rank, int numWorkers,int *shared_neighbours);
        bool thermalized;
        
        int worker_dn;
        int worker_up;
        int world_rank;                   //Thread number
        int world_size;
        int T_id;
        model modelo;    
        std::vector<double> energies;
        std::vector<int> magnetizations;
        std::vector<int> accept_timeseries;
        std::vector<std::vector<double>> correlations;
        std::vector<int> mcs_timeseries;
        std::vector<double> counters_t;
        std::vector<double> t_timeseries;
        std::vector<double> temperatures;
        double T;

        void thermalization(double temp, int *shared_neighbours);
        void sampling(int *shared_neighbours,int *shared_distanceMatrix, bool storeCorrelations, bool ladderUpdate);
        void sweep(double temp,int *shared_neighbours, bool counterFlag); 
        void swap_workers();
};
#endif