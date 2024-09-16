#include "include/worker.h"
#include <iostream>

worker::worker(){}

worker::worker(int rank, int numWorkers, int *shared_neighbours)
:thermalized(false),
world_rank(rank),
world_size(numWorkers),
T_id(rank)
{
    if (numWorkers <= 0)
    {
        throw std::invalid_argument("numWorkers must be greater than 0");
    }

    if (shared_neighbours == nullptr)
    {
        throw std::invalid_argument("shared_neighbours cannot be null");
    }

    worker_dn = settings::sim::mod(rank - 1, world_size);
    worker_up = settings::sim::mod(rank + 1, world_size);

    temperatures.reserve(numWorkers);
    

    for (int idx = 0; idx < numWorkers; idx++)
    {
        temperatures.push_back(settings::sim::T_min * (pow(settings::sim::ratio, 1.0 * idx / (1 - numWorkers)))); // idx*1.0*(settings::sim::T_max-settings::sim::T_min)/(numWorkers-1)+settings::sim::T_min);
    }

    temperatures = {
        0.01, 0.131375925, 0.2192512725, 0.2897279738, 0.3408229442, 0.3796016988, 0.4089561089, 0.4348109125, 0.459629004, 0.4866419968, 0.5164514458, 0.5499549047, 0.5857143595, 0.6238921119, 0.664056418, 0.7063748626, 0.7508436232, 0.79682986, 0.843975676, 0.8925446972, 0.9432834456, 0.9968726634, 1.05371731, 1.11420416, 1.1788821600000001, 1.248377824, 1.323236132, 1.40397498, 1.49119837, 1.5855957760000001, 1.687754818, 1.797967044, 1.916501189, 2.044449512, 2.184010185, 2.33889142, 2.514784558, 2.7153148959999998, 2.935504338, 3.175375747, 3.451545256, 3.775331019, 4.17206792, 4.627818632, 5.1152590920000005, 5.697152635, 6.525809684, 7.80650979, 9.187644902999999, 10.0
    };

    compute_probabilities();
    start_counters();
    modelo = model(rank, shared_neighbours);
    std::cout<<"Worker "<<rank<<" created"<<std::endl;

    T = temperatures[rank];
    cooldown(shared_neighbours);
    thermalization(T, shared_neighbours);
    std::cout<<"Worker "<<rank<<" done thermalizing"<<std::endl;
}

void worker::compute_probabilities()
{
    for(int temperature_idx = 0; temperature_idx<world_size; temperature_idx++)
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
    counter::swap_trials        = 0;
    counter::swap_accepts       = 0;
}

void worker::cooldown(int *shared_neighbours)
{
    double highTemp = 15;
    while(highTemp > T)
    {
        thermalization(highTemp,shared_neighbours);
        highTemp-=0.5;
    }
    while(highTemp > T)
    {
        thermalization(highTemp,shared_neighbours);
        highTemp-=0.1;
    }
    thermalized = false;
}

bool worker::performTrialMove(double temp,int trialSite, int *shared_neighbours)
{
    modelo.trialMove(trialSite,shared_neighbours);
    if(modelo.delE<=0)
    {   
        modelo.acceptMove(trialSite);
        return true;
    }
    else
    {
        double rn = rn_gen::rand_double();
        double prob = probabilities_list[T_id][int(modelo.delE/2)+7];
        
        if(rn<prob)
        {
            modelo.acceptMove(trialSite);
            return true;
        }
    }
    return false;
}

void worker::sweep(double temp,int *shared_neighbours,bool counterFlag)
{
    for(int idx = 0; idx<modelo.nSpins; idx++)
    {
        counter::MCS += 1;
        int trialSite = rn_gen::rand_site();
        bool movedDone = performTrialMove(temp,trialSite,shared_neighbours);
        if(movedDone && counterFlag)
        {
            counter::accepts += 1;
        }
    }
}

void worker::thermalization(double temp,int *shared_neighbours)
{
    for(int therm_step = 0; therm_step < settings::sim::MCS_therm; therm_step ++)
    {sweep(temp,shared_neighbours,false);}
    thermalized = true;
}

void worker::storeThermodynamicData(bool storeCorrelations,int *distanceMatrix)
{
    if(storeCorrelations)
    {
        std::vector<double> corr = modelo.compute_Correlations(distanceMatrix);
        correlations.push_back(corr);
    }
    energies.push_back(modelo.E);
    magnetizations.push_back(modelo.M);
    t_timeseries.push_back(T);
}

void worker::storeCounters()
{
    accept_timeseries.push_back(counter::swap_accepts);
    mcs_timeseries.push_back(counter::swap_trials);
    counters_t.push_back(T);  
    counter::accepts = 0;
    counter::MCS = 0;
}

void worker::sampling(int *shared_neighbours,int *distanceMatrix, bool storeCorrelations, bool ladderUpdate)
{
    for(int samp_step = 0; samp_step <settings::sim::MCS_sample; samp_step ++)
    {
        for(int wait = 0; wait<settings::sim::MCS_decorr; wait++)
        {sweep(T,shared_neighbours,true);}
        storeThermodynamicData(storeCorrelations,distanceMatrix);
    }
    if(ladderUpdate)
    {
        storeCounters();
    }
}


void worker::swap_workers()
{
    settings::counter::swap_trials +=1;
    int swap;
    double E_up;
    double P_swap;
    bool myTurn = settings::sim::mod(T_id,2) == settings::sim::mod(settings::counter::swap_trials,2);
    MPI_Sendrecv(&modelo.E, 1, MPI_DOUBLE, worker_dn,100,&E_up, 1, MPI_DOUBLE, worker_up,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

 if(myTurn)
    {
        double T_up = temperatures[T_id+1];
        P_swap = exp((modelo.E-E_up)*(1/T-1/T_up));

        swap = rn_gen::rand_double() < fmin(1,P_swap) ? 1 : 0;
        if(T_id == world_size-1){swap = 0;}
        MPI_Send(&swap, 1, MPI_INT, worker_up, 110, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(&swap, 1, MPI_INT, worker_dn, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    int world_ID_dn_new = worker_dn, world_ID_up_new = worker_up;
    if(myTurn)
    {
        //Receive news from world_ID_dn about who will be your new world_ID_dn
        MPI_Recv(&world_ID_dn_new, 1, MPI_INT, worker_dn, worker_dn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
    //Send news up to world_ID_up that world_ID_dn will be at this position
        if (swap)
        {
            MPI_Send(&worker_dn, 1, MPI_INT, worker_up, world_rank, MPI_COMM_WORLD);
        } 
        else
        {
            MPI_Send(&world_rank,1, MPI_INT, worker_up,world_rank, MPI_COMM_WORLD);
        }
    }
    if(myTurn)
    {
        if (swap) 
        {
            //Send up
            MPI_Send(&world_ID_dn_new, 1, MPI_INT, worker_up, world_rank, MPI_COMM_WORLD);
            //Now I can correct (*)
            world_ID_dn_new = worker_up;
        }
        else
        {
            MPI_Send(&world_rank, 1, MPI_INT, worker_up, world_rank, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&world_ID_dn_new, 1, MPI_INT, worker_dn, worker_dn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(myTurn)
    {
        if(swap)
        {
            MPI_Send(&worker_up, 1, MPI_INT, worker_dn, world_rank, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Send(&world_rank, 1, MPI_INT, worker_dn, world_rank, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&world_ID_up_new, 1, MPI_INT, worker_up, worker_up, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }

    if(myTurn)
    {
        MPI_Recv(&world_ID_up_new, 1, MPI_INT, worker_up, worker_up, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        if (swap)
        {
            MPI_Send(&world_ID_up_new, 1, MPI_INT, worker_dn, world_rank, MPI_COMM_WORLD);
            world_ID_up_new = worker_dn;
        }
        else
        {
            MPI_Send(&world_rank,    1, MPI_INT, worker_dn, world_rank, MPI_COMM_WORLD);
        }
    }
    
    if (swap == 1)
    {
        if (myTurn == 1)
        {
            T_id = settings::sim::mod(T_id + 1, world_size);
            T = temperatures[T_id];
            MPI_Sendrecv_replace(&counter::accepts,1,MPI_INT, worker_up, 1, worker_up, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
                T_id = settings::sim::mod(T_id - 1, world_size);
                T = temperatures[T_id];
                MPI_Sendrecv_replace(&counter::accepts,1,MPI_INT, worker_dn, 1, worker_dn, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    worker_dn = world_ID_dn_new;
    worker_up = world_ID_up_new;
    
    counter::swap_accepts += swap;
    counter::swap_trials++;

}
