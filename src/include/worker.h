#ifndef WORKER_H
#define WORKER_H
#include "model.h"
#include <mpi.h>
#include <math.h>
#include "counters.h"



class worker
{
    public:
        worker();
        worker(int rank, int numWorkers);
        int world_rank;                   //Thread number
        int world_size;
        bool thermalized;
        std::vector<double> e_timeseries;
        std::vector<int> magn_timeseries;
        std::vector<double> t_timeseries;
        std::vector<double> t_ladder;
        std::vector<double> temperatures;
        int T_id;
        double T;
        int worker_dn;
        int worker_up;
        int acceptedSteps;

        std::vector<int> T_id_list;
        void cooldown();
        void start_counters();
        void sweep(double temp);
        void thermalization(double temp);
        void sampling();
        void broadcastDecision(bool shouldSwap);
        void swapConfigurations();
        int numWorkers;
        model modelo;
        int rank;

};

void swap_workers(worker &work, std::vector<worker> &workers);

#endif