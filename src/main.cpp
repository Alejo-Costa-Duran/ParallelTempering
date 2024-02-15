#include <iostream>
#include <mpi.h>
#include "./include/rngenerator.h"
#include "./include/worker.h"
#include "./include/counters.h"
#include <fstream>
#include <string>


void swap_workers(worker &work)
{
    settings::counter::swap_trials +=1;
    int swap;
    double E_up;
    double P_swap;
    bool myTurn = settings::sim::mod(work.T_id,2) == settings::sim::mod(settings::counter::swap_trials,2);
    MPI_Sendrecv(&work.modelo.E, 1, MPI_DOUBLE, work.worker_dn,100,&E_up, 1, MPI_DOUBLE, work.worker_up,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

 if(myTurn)
    {
        double T_up = work.temperatures[settings::sim::mod(work.worker_up,work.world_size)];
        P_swap = exp((work.modelo.E-E_up)*(1/work.T-1/T_up));
        swap = rn_gen::rand_double() < fmin(1,P_swap) ? 1 : 0;
        if(work.T_id == work.world_size-1){swap = 0;}
        MPI_Send(&swap, 1, MPI_INT, work.worker_up, 110, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(&swap, 1, MPI_INT, work.worker_dn, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    int world_ID_dn_new = work.worker_dn, world_ID_up_new = work.worker_up;
    if(myTurn)
    {
        //Receive news from world_ID_dn about who will be your new world_ID_dn
        MPI_Recv(&world_ID_dn_new, 1, MPI_INT, work.worker_dn, work.worker_dn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
    //Send news up to world_ID_up that world_ID_dn will be at this position
        if (swap)
        {
            MPI_Send(&work.worker_dn, 1, MPI_INT, work.worker_up, work.world_rank, MPI_COMM_WORLD);
        } 
        else
        {
            MPI_Send(&work.world_rank,1, MPI_INT, work.worker_up,work.world_rank, MPI_COMM_WORLD);
        }
    }
    if(myTurn)
    {
        if (swap) 
        {
            //Send up
            MPI_Send(&world_ID_dn_new, 1, MPI_INT, work.worker_up, work.world_rank, MPI_COMM_WORLD);
            //Now I can correct (*)
            world_ID_dn_new = work.worker_up;
        }
        else
        {
            MPI_Send(&work.world_rank, 1, MPI_INT, work.worker_up, work.world_rank, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&world_ID_dn_new, 1, MPI_INT, work.worker_dn, work.worker_dn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(myTurn)
    {
        //Send news dn to world_ID_dn who will be his new neighbor up
        if(swap)
        {
            MPI_Send(&work.worker_up, 1, MPI_INT, work.worker_dn, work.world_rank, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Send(&work.world_rank, 1, MPI_INT, work.worker_dn, work.world_rank, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&world_ID_up_new, 1, MPI_INT, work.worker_up, work.worker_up, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //(**)I know that if I'm swapping the variable above is not actually correct yet.
            //world_ID_up_new will actually be worker.world_ID_dn
    }

    if(myTurn)
    {
        //Receive news from world_ID_up about who will be your new world_ID_up
        MPI_Recv(&world_ID_up_new, 1, MPI_INT, work.worker_up, work.worker_up, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
            //Send news up to world_ID_dn who will be at this position
        if (swap)
        {
            MPI_Send(&world_ID_up_new, 1, MPI_INT, work.worker_dn, work.world_rank, MPI_COMM_WORLD);
                //Now I can correct (**)
            world_ID_up_new = work.worker_dn;
        }
        else
        {
            MPI_Send(&work.world_rank,    1, MPI_INT, work.worker_dn, work.world_rank, MPI_COMM_WORLD);
        }
    }
    //Now do the swapping if you got lucky
    if (swap == 1)
    {
        if (myTurn == 1)
        {
            work.T_id = settings::sim::mod(work.T_id + 1, work.world_size);
            work.T = work.temperatures[work.T_id];
            MPI_Sendrecv_replace(&counter::accepts,1,MPI_INT, work.worker_up, 1, work.worker_up, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
                work.T_id = settings::sim::mod(work.T_id - 1, work.world_size);
                work.T = work.temperatures[work.T_id];
                MPI_Sendrecv_replace(&counter::accepts,1,MPI_INT, work.worker_dn, 1, work.worker_dn, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    if(swap == 1)
    {
        work.thermalized = false;
    }

    work.worker_dn = world_ID_dn_new;
    work.worker_up = world_ID_up_new;
    
    MPI_Allgather(&work.T_id,1,MPI_INT,work.T_id_list.data(),1, MPI_INT, MPI_COMM_WORLD);

    counter::swap_accepts += swap;
    counter::swap_trials++;

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


    worker work(world_rank,world_size);
    for(int iteration = 0; iteration<settings::sim::MCS_total; iteration++)
    {
        if(iteration%settings::sim::MCS_swap == 0){swap_workers(work);}
        if(!work.thermalized){work.thermalization(work.T);}
        work.sampling();
    }
    std::cout<<"Proceso "<<world_rank<< " intentó hacer " << counter::swap_trials <<"cambios de temperatura y logró "<<counter::swap_accepts<<"\n";
    std::ofstream file("../PT-Data/DatosProcesoN15" + std::to_string(world_rank) + ".csv",std::ios::app);
    //file<<"Energia\tMagnetizacion\tTemperatura\n";
        for(int idx = 0; idx < work.e_timeseries.size(); idx++)
        {
            file<<work.e_timeseries[idx]<<"\t"<<work.magn_timeseries[idx]<<"\t"<<work.t_timeseries[idx]<<"\n";
        }
    file.close();
    // Get the name of the processor

    // Print off a hello world message

    // Finalize the MPI environment
    MPI_Finalize();
}