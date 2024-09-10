#ifndef WORKER_H
#define WORKER_H
#include "model.h"
#include <mpi.h>
#include <math.h>
#include "counters.h"
#include "settings.h"



class worker
{
    public:
        worker();
        worker(int rank, int numWorkers,int *shared_neighbours);
        int world_rank;                   //Thread number
        int world_size;
        bool thermalized;
        std::vector<double> energies;
        std::vector<int> magnetizations;
        std::vector<int> accept_timeseries;
        std::vector<std::vector<double>> correlations;
        std::vector<int> mcs_timeseries;
        std::vector<double> counters_t;
        std::vector<double> t_timeseries;
        std::vector<double> t_ladder;
        std::vector<double> temperatures;
        std::vector<std::vector<double>> probabilities_list;
        int T_id;
        double T;
        int worker_dn;
        int worker_up;
        int acceptedSteps;

        bool ladderUpdate;
        std::vector<int> T_id_list;
        void cooldown(int *shared_neighbours);
        void start_counters();
        void sweep(double temp,int *shared_neighbours);
        void thermalization(double temp, int *shared_neighbours);
        void sampling(int *shared_neighbours,int *shared_distanceMatrix);
        void compute_probabilities();
        int numWorkers;
        model modelo;

};

void swap_workers(worker &work, std::vector<worker> &workers);

#endif