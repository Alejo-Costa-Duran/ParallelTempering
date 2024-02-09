#include "include/worker.h"

worker::worker(){}

worker::worker(int rank, int numWorkers)
{
    modelo = model(rank);
    world_rank = rank;
    world_size = numWorkers;
    acceptedSteps = 0;
    T = settings::sim::T_min+(rank-1)*settings::sim::T_step;
    T_id = rank;
    thermalized = false;
    worker_dn = settings::sim::mod(rank-1,world_size);
    worker_up = settings::sim::mod(rank+1,world_size);
    start_counters();
}

void worker::start_counters() {
    counter::MCS                = 0;
    counter::accepts            = 0;
    counter::trials             = 0;
    counter::samples            = 0;
    counter::swap_trials        = 0;
    counter::swap_accepts       = 0;
    counter::store              = 0;
    timer::cout                 = 0;
    timer::comp                 = 0;
    timer::swap 				= 0;
    timer::move 				= 0;
    timer::sync 				= 0;
}



void worker::sweep()
{
    for(int idx = 0; idx<modelo.nSpins; idx++)
    {
        int trialSite = rn_gen::rand_site();
        modelo.trialMove(trialSite);
        if(modelo.delE<0)
        {   
            counter::accepts+=1;
            modelo.acceptMove(trialSite);
        }
        else
        {
            double rn = rn_gen::rand_double();
            double prob = exp(-modelo.delE/T);
            if(rn<prob)
            {
                counter::accepts+=1;
                modelo.acceptMove(trialSite);
            }
        }
    }
    counter::MCS+=1;
}

void worker::thermalization()
{
    for(int therm_step = 0; therm_step < settings::sim::MCS_therm; therm_step ++)
    {
        sweep();
    }
    thermalized = true;
}

void worker::sampling()
{
    for(int samp_step = 0; samp_step <settings::sim::MCS_sample; samp_step ++)
    {
        sweep();
        e_timeseries.push_back(modelo.E);
        magn_timeseries.push_back(modelo.M);
        t_timeseries.push_back(T);
    }
}
