#include "include/worker.h"
#include <iostream>

worker::worker(){}

worker::worker(int rank, int numWorkers)
{
    for(int idx=0; idx<numWorkers; idx++)
    {
        t_ladder.push_back(settings::sim::T_min*(pow(settings::sim::ratio,1.0*idx/(1-numWorkers))));
    }
    for(int idx = 0; idx<numWorkers; idx++)
    {
        temperatures.push_back(t_ladder[idx]);
        T_id_list.push_back(idx);
    }
    modelo = model(rank);
    world_rank = rank;
    world_size = numWorkers;
    acceptedSteps = 0;
    T = temperatures[rank];
    T_id = rank;
    thermalized = false;
    worker_dn = settings::sim::mod(rank-1,world_size);
    worker_up = settings::sim::mod(rank+1,world_size);
    start_counters();
    cooldown();
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

void worker::cooldown()
{
    double highTemp = 5;
    while(highTemp > T)
    {
        thermalization(highTemp);
        highTemp-=0.1;
    }
    thermalized = false;
}

void worker::sweep(double temp)
{
    for(int idx = 0; idx<modelo.nSpins; idx++)
    {
        counter::MCS += 1;
        int trialSite = rn_gen::rand_site();
        modelo.trialMove(trialSite);
        if(modelo.delE<=0)
        {   
            counter::accepts+=1;
            modelo.acceptMove(trialSite);
        }
        else
        {
            double rn = rn_gen::rand_double();
            double prob = exp(-modelo.delE/temp);
            if(rn<prob)
            {
                counter::accepts+=1;
                modelo.acceptMove(trialSite);
            }
        }
    }
}

void worker::thermalization(double temp)
{
    for(int therm_step = 0; therm_step < settings::sim::MCS_therm; therm_step ++)
    {
        for(int idx = 0; idx<modelo.nSpins; idx++)
    {
        int trialSite = rn_gen::rand_site();
        modelo.trialMove(trialSite);
        if(modelo.delE<=0)
        {   
            modelo.acceptMove(trialSite);
        }
        else
        {
            double rn = rn_gen::rand_double();
            double prob = exp(-modelo.delE/temp);
            if(rn<prob)
            {
                modelo.acceptMove(trialSite);
            }
        }
    }
    }
    thermalized = true;
}

void worker::sampling()
{
    for(int samp_step = 0; samp_step <settings::sim::MCS_sample; samp_step ++)
    {
        for(int wait = 0; wait<settings::sim::MCS_decorr; wait++)
        {sweep(T);}
        e_timeseries.push_back(modelo.E);
        magn_timeseries.push_back(modelo.M);
        t_timeseries.push_back(T);
    }
}
