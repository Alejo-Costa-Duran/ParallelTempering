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
        temperatures.push_back(settings::sim::T_min);
        //settings::sim::T_min * (pow(settings::sim::ratio, 1.0 * idx / (1 - numWorkers)))); // idx*1.0*(settings::sim::T_max-settings::sim::T_min)/(numWorkers-1)+settings::sim::T_min);
    }
    /*temperatures = {0.01, 0.131376, 0.219251, 0.289728, 0.340823, 0.379602, 0.408956, 0.434811, 0.459629, 0.486642, 0.516451, 0.549955, 0.585714, 0.623892, 0.664056, 0.706375, 0.750844, 0.79683, 0.843976, 0.892545, 0.943283, 0.996873, 1.05372, 1.1142, 1.17888, 1.24838, 1.32324, 1.40397, 1.4912, 1.5856, 1.68775, 1.79797, 1.9165, 2.04445, 2.18401, 2.33889, 2.51478, 2.71531, 2.9355, 3.17538, 3.45155, 3.77533, 4.17207, 4.62782, 5.11526, 5.69715, 6.52581, 7.80651, 9.18764, 10.0
    };*/

    //compute_probabilities();
    start_counters();
    modelo = model(rank, shared_neighbours);

    T = temperatures[rank];
    //cooldown(shared_neighbours);
    lowestEnergy = modelo.E;
    lowestMagnetization= modelo.M;
    thermalization(T, shared_neighbours);

}

void worker::compute_probabilities()
{
    for(int temperature_idx = 0; temperature_idx<world_size; temperature_idx++)
    {
        std::vector<double> prob_list(16,0);
        double currentTemperature = temperatures[temperature_idx];
        for(int neighbourSum = -7; neighbourSum<8; neighbourSum+=2 )
        {
            double delE = neighbourSum-modelo.field;
            int probIdx = int((neighbourSum+7)/2);
            prob_list[probIdx] = exp(-2.0*delE/currentTemperature);
            prob_list[8+probIdx] = exp(2.0*delE/currentTemperature);
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
    double highTemp = 4.0;
    while(highTemp > T)
    {
        thermalization(highTemp,shared_neighbours);
        highTemp=highTemp*0.9;
    }
}

bool worker::performTrialMove(double temp,int trialSite, int *shared_neighbours)
{
    modelo.trialMove(trialSite,shared_neighbours);
    //int trialSpin = modelo.lattice[trialSite];
    if(modelo.delE<=0)
    {   
        modelo.acceptMove(trialSite);
        return true;
    }
    else
    {
        double rn = rn_gen::rand_double();
        //int probIdx = (1.0*modelo.delE/(-2.0*trialSpin)+modelo.field+7)/2+(trialSpin+1)*4;
        double prob = exp(-modelo.delE/temp); //probabilities_list[T_id][probIdx];
        //std::cout<<prob<<"\t"<<exp(-modelo.delE/temp)<<"\t"<<temp<<"\t"<<temperatures[T_id]<<"\n";
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
        if(modelo.E<lowestEnergy)
        {
            lowestEnergy = modelo.E;
            lowestMagnetization = modelo.M;
        }
       // if(movedDone && counterFlag)
        //{
          //  counter::accepts += 1;
            //if(modelo.E < lowestEnergy)
            //{
              //  lowestEnergy = modelo.E;
                //lowestMagnetization = modelo.M;
            //}
        //}
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
