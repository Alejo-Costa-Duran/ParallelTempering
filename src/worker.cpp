#include "include/worker.h"
#include <mpi.h>

worker::worker(int rank, int numWorkers)
{
    world_rank = rank;
    world_size = numWorkers;
    acceptedSteps = 0;
    T = settings::sim::T_max;
}

void worker::sweep()
{
    for(int idx = 0; idx<modelo.nSpins; idx++)
    {
        int trialSite = rn_gen::rand_site();
        modelo.trialMove(trialSite);
        if(modelo.delE<0)
        {   
            acceptedSteps+=1;
            modelo.acceptMove(trialSite);
        }
        else
        {
            double rn = rn_gen::rand_double();
            double prob = exp(-modelo.delE/T);
            if(rn<prob)
            {
                acceptedSteps+=1;
                modelo.acceptMove(trialSite);
            }

        }
        
    }
}