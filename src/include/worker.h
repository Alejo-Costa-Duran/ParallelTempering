#ifndef WORKER
#define WORKER
#include "model.h"

class worker
{
    public:
        worker(int rank, int numWorkers);
        int world_rank;                   //Thread number
        int world_size;                 //Total number of threads
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