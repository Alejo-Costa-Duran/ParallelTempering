#ifndef WORKER
#define WORKER

class worker
{
    public:
        worker(int rank, int numWorkers);
        void runSimulation();
        void gatherData();
        void broadcastDecision(bool shouldSwap);
        void swapConfigurations();
    private:
        int rank;
        int numWorkers;
};



#endif