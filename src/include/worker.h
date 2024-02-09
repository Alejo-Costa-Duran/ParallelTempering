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
        int T_id;
        double T;
        int worker_dn;
        int worker_up;
        int acceptedSteps;

        std::vector<double> temperatures;
        std::vector<int> T_id_list;

        void start_counters();
        void sweep();
        void thermalization();
        void sampling();
        void broadcastDecision(bool shouldSwap);
        void swapConfigurations();
        model modelo;
    private:
        int rank;
        int numWorkers;
};

void swap_workers(worker &work, std::vector<worker> &workers);

#endif