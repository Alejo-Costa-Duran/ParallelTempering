#include "include/worker.h"
#include <iostream>

worker::worker(){}

worker::worker(int rank, int numWorkers,int *shared_neighbours)
{
    world_rank = rank;
    world_size = numWorkers;
    worker_dn = settings::sim::mod(rank-1,world_size);
    worker_up = settings::sim::mod(rank+1,world_size);

    for(int idx=0; idx<numWorkers; idx++)
    {
        t_ladder.push_back(settings::sim::T_min*(pow(settings::sim::ratio,1.0*idx/(1-numWorkers))));//idx*1.0*(settings::sim::T_max-settings::sim::T_min)/(numWorkers-1)+settings::sim::T_min);
    }
    t_ladder = {0.01, 0.1131375925, 0.2192512725, 0.2897279738, 0.3408229442, 0.3796016988, 0.4089561089, 0.4348109125, 0.459629004, 0.4866419968, 0.5164514458, 0.5499549047, 0.5857143595, 0.6238921119, 0.664056418, 0.7063748626, 0.7508436232, 0.79682986, 0.843975676, 0.8925446972, 0.9432834456, 0.9968726634, 1.05371731, 1.11420416, 1.1788821600000001, 1.248377824, 1.323236132, 1.40397498, 1.49119837, 1.5855957760000001, 1.687754818, 1.797967044, 1.916501189, 2.044449512, 2.184010185, 2.33889142, 2.514784558, 2.7153148959999998, 2.935504338, 3.175375747, 3.451545256, 3.775331019, 4.17206792, 4.627818632, 5.1152590920000005, 5.697152635, 6.525809684, 7.80650979, 9.187644902999999, 10.0};
    for(int idx = 0; idx<numWorkers; idx++)
    {
        temperatures.push_back(t_ladder[idx]);
        T_id_list.push_back(idx);
    }

    compute_probabilities();
    start_counters();
    modelo = model(rank,shared_neighbours);

    acceptedSteps = 0;
    T = temperatures[rank];
    T_id = rank;
    ladderUpdate = false;
    thermalized = false;


    cooldown(shared_neighbours);
    thermalization(T,shared_neighbours);
}

void worker::compute_probabilities()
{
    for(int temperature_idx = 0; temperature_idx<numWorkers; temperature_idx++)
    {
        std::vector<double> prob_list(15,0);
        double currentTemperature = temperatures[temperature_idx];
        for(int neighbourSum = 0; neighbourSum<15; neighbourSum++)
        {
            prob_list[neighbourSum] = exp(-2.0*(-7+neighbourSum)/currentTemperature);
        } 
        probabilities_list.push_back(prob_list);
    }
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

void worker::cooldown(int *shared_neighbours)
{
    double highTemp = 15;
    while(highTemp > T)
    {
        thermalization(highTemp,shared_neighbours);
        highTemp-=0.5;
    }
    thermalized = false;
}

void worker::sweep(double temp,int *shared_neighbours)
{
    for(int idx = 0; idx<modelo.nSpins; idx++)
    {
        counter::MCS += 1;
        int trialSite = rn_gen::rand_site();
        modelo.trialMove(trialSite,shared_neighbours);
        if(modelo.delE<=0)
        {   
            counter::accepts+=1;
            modelo.acceptMove(trialSite);
        }
        else
        {
            double rn = rn_gen::rand_double();
            double prob = probabilities_list[T_id][int(modelo.delE/2)+7];
            
            if(rn<prob)
            {
                counter::accepts+=1;
                modelo.acceptMove(trialSite);
            }
        }
    }
}

void worker::thermalization(double temp,int *shared_neighbours)
{
    for(int therm_step = 0; therm_step < settings::sim::MCS_therm; therm_step ++)
    {
        for(int idx = 0; idx<modelo.nSpins; idx++)
    {
        int trialSite = rn_gen::rand_site();
        modelo.trialMove(trialSite,shared_neighbours);
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

void worker::sampling(int *shared_neighbours,int *distanceMatrix)
{
    for(int samp_step = 0; samp_step <settings::sim::MCS_sample; samp_step ++)
    {
        for(int wait = 0; wait<settings::sim::MCS_decorr; wait++)
        {sweep(T,shared_neighbours);}
        //std::vector<double> corr = modelo.compute_Correlations(distanceMatrix);
        //corr_timeseries.push_back(corr);
        energies.push_back(modelo.E);
        magnetizations.push_back(modelo.M);
        t_timeseries.push_back(T);
    }
    if(ladderUpdate)
    {
    accept_timeseries.push_back(counter::swap_accepts);
    mcs_timeseries.push_back(counter::swap_trials);
    counters_t.push_back(T);  
    counter::accepts = 0;
    counter::MCS = 0;
    }
}
