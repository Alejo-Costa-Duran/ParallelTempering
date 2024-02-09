#include <iostream>
#include <mpi.h>
#include "./include/rngenerator.h"
#include "./include/worker.h"
#include "./include/counters.h"


void swap_workers(worker &work, std::vector<worker> &workers)
{
    settings::counter::swap_trials +=1;
    int swap;
    double E_up;
    double P_swap;
    bool myTurn = settings::sim::mod(work.T_id,2) == settings::sim::mod(settings::counter::swap_trials,2); 
    MPI_Sendrecv(&work.modelo.E, 1, MPI_DOUBLE, work.worker_dn,100,&E_up, 1, MPI_DOUBLE, work.worker_up,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::cout<<E_up<<"\t"<<workers[work.worker_up].modelo.E<<"\t"<<myTurn<<"\n";
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    std::vector<worker> workers(world_size);
    for(int w_id = 0; w_id<world_size; w_id++)
    {
        workers[w_id] = worker(w_id,world_size);
    }
    for(int iteration = 0; iteration<settings::sim::MCS_total; iteration++)
    {
        workers[world_rank].sweep();
    }
    std::cout<<"Worker number "<<world_rank<<" has accepted "<<counter::accepts<<" steps \n";

    // Get the name of the processor

    // Print off a hello world message

    // Finalize the MPI environment
    MPI_Finalize();
}