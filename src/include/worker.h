#ifndef WORKER
#define WORKER
#include "model.h"
#include <math.h>

class worker
{
    public:
        worker(int rank, int numWorkers);
        int world_rank;                   //Thread number
        int world_size;
        int T_id;
        double T;
        int acceptedSteps;

        std::vector<double> temperatures;
        std::vector<int> T_id_list;

        void step();
        void sweep();
        void gatherData();
        void broadcastDecision(bool shouldSwap);
        void swapConfigurations();
        model modelo;
    private:
        int rank;
        int numWorkers;
};



#endif