#include <iostream>
#include <mpi.h>
#include "./include/vecinos.h"
#include "./include/rngenerator.h"
#include "./include/worker.h"
#include "./include/settings.h"
#include "./include/counters.h"
#include <fstream>
#include <string>
#include <queue>
#include <limits> // For infinity


std::vector<int> computeDistanceMatrix(std::vector<int> neighboursMatrix, int nSpins) {
    std::vector<int> distance_matrix(nSpins*(nSpins+1)/2,std::numeric_limits<int>::max()); // Initialize with large value (infinity)

    // Perform BFS from each spin to calculate distances
    for (int start = 0; start < nSpins; ++start) {
        std::queue<std::pair<int, int>> q; // Pair (spin, distance)
        std::vector<bool> visited(nSpins, false);

        // Start BFS from 'start' spin
        q.push({start, 0});
        int idx = start*(start+1)/2+start;
        distance_matrix[idx] = 0;
        visited[start] = true;

        while (!q.empty()) {
            int current = q.front().first;
            int dist = q.front().second;
            q.pop();

            // Check all neighbors of the current spin
            for (int n = 0; n < 7; ++n) {
                int neighbor = neighboursMatrix[current * 7 + n]; // Get neighbor spin
                
                // If neighbor hasn't been visited, update its distance
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    int idx =(start > neighbor) ? start * (start + 1) / 2 + neighbor : neighbor * (neighbor + 1) / 2 + start; ;
                    distance_matrix[idx] = dist + 1;
                    q.push({neighbor, dist + 1});
                }
            }
        }
    }

    return distance_matrix;
}

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
        double T_up = work.temperatures[work.T_id+1];
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

    work.worker_dn = world_ID_dn_new;
    work.worker_up = world_ID_up_new;
    
    //MPI_Allgather(&work.T_id,1,MPI_INT,work.T_id_list.data(),1, MPI_INT, MPI_COMM_WORLD);

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
    MPI_Comm node_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, world_rank, MPI_INFO_NULL, &node_comm);

    int local_rank, local_size;
    MPI_Comm_rank(node_comm, &local_rank);
    MPI_Comm_size(node_comm, &local_size);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    std::cout << "Process " << world_rank << " (local rank " << local_rank
              << ") is running on node: " << processor_name << std::endl;

    // Define variables for the shared memory window
    int* shared_neighbours;
    int* shared_distanceMatrix;
    MPI_Win win_neighbours, win_distance;

    std::vector<int> neighbours(settings::model::nSpins*7);
    std::vector<int> distanceMatrix(settings::model::nSpins*(settings::model::nSpins+1)/2);
   // Rank 0 in each node will allocate memory for the neighbours and distance matrices
    int nSpins = settings::model::nSpins;
    int neighbours_size = nSpins * 7;
    int distanceMatrix_size = (nSpins * (nSpins + 1)) / 2;

    // Allocate shared memory window for neighbors matrix
    if (local_rank == 0) {
        // Rank 0 on each node allocates memory
        MPI_Win_allocate_shared(neighbours_size * sizeof(int), sizeof(int), MPI_INFO_NULL,
                                node_comm, &shared_neighbours, &win_neighbours);

        MPI_Win_allocate_shared(distanceMatrix_size * sizeof(int), sizeof(int), MPI_INFO_NULL,
                                node_comm, &shared_distanceMatrix, &win_distance);

        // Only rank 0 loads the data
        std::ifstream file("./src/Vecinos.csv");
        if (settings::model::neighIdx == 3) {
            int number;
            std::string line;
            int idx = 0;
            while (std::getline(file, line)) {
                number = std::stoi(line);
                shared_neighbours[idx] = number;
                idx++;
            }
        } else {
            // Assign precomputed neighbour data
            for (int i = 0; i < neighbours_size; i++) {
                shared_neighbours[i] = vecinos[settings::model::neighIdx][i];
            }
        }
        file.close();
        
        // Compute distance matrix
        std::vector<int> distanceMatrix = computeDistanceMatrix(std::vector<int>(shared_neighbours, shared_neighbours + neighbours_size), nSpins);
        std::copy(distanceMatrix.begin(), distanceMatrix.end(), shared_distanceMatrix);
        std::cout << "Rank 0 on node " << processor_name << " initialized the neighbours and distance matrices.\n";
    } else {
        // Other ranks on the node attach to the shared memory
        MPI_Aint size;
        int disp_unit;

        // Attach to the shared memory window for neighbors matrix
        MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, node_comm, &shared_neighbours, &win_neighbours);
        MPI_Win_shared_query(win_neighbours, 0, &size, &disp_unit, &shared_neighbours);

        // Attach to the shared memory window for distance matrix
        MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, node_comm, &shared_distanceMatrix, &win_distance);
        MPI_Win_shared_query(win_distance, 0, &size, &disp_unit, &shared_distanceMatrix);
    }

    // Synchronize after data load
    MPI_Barrier(MPI_COMM_WORLD);
    worker work(world_rank,world_size, shared_neighbours);

        for(int iteration = 0; iteration<settings::sim::MCS_total; iteration++)
        {
            
            if(!work.thermalized){work.thermalization(work.T,shared_neighbours);};
            if(iteration%settings::sim::MCS_swap == 0){work.swap_workers();}
            work.sampling(shared_neighbours,shared_distanceMatrix,settings::sim::storeCorrelations,settings::sim::ladderUpdate);
            if(work.world_rank==1){if(iteration%10000 == 0){std::cout<<iteration<<"\n";}}
        }

    std::cout<<"Proceso "<<world_rank<< " intentó hacer " << counter::swap_trials <<"cambios de temperatura y logró "<<counter::swap_accepts<<"\n";

    std::string fileName = "../PT-Data/0Field/DatosN"+std::to_string(settings::model::nClusters)+"Proceso"+ std::to_string(world_rank) + ".csv";
    std::ofstream file(fileName);
    file<<"Energia\tMagnetizacion\tTemperatura";
    if(settings::sim::storeCorrelations)
    {
        for(int corrIdx = 0; corrIdx < settings::model::distances+1; corrIdx++)
        {
            file<<"\tCorr"+std::to_string(corrIdx);
        }
        file<<"\n";
    }
    else
    {
        file<<"\n";
    }

    for(int idx = 0; idx < work.energies.size(); idx++)
    {
        file<<work.energies[idx]<<"\t"<<work.magnetizations[idx]<<"\t"<<work.t_timeseries[idx];
        if(settings::sim::storeCorrelations)
        {
        for(int corrIdx = 0; corrIdx<settings::model::distances+1; corrIdx++)
        {
            file<<"\t"<<work.correlations[idx][corrIdx];
        }}
        file<<"\n";
    }
    file.close();
    
    if(settings::sim::ladderUpdate)
    {
    std::string counterName = "../PT-Data/0Field/CountersProceso"+std::to_string(settings::model::nClusters) + std::to_string(world_rank) + ".csv";
    std::ofstream fileCounters(counterName);
    fileCounters<<"Accepts\tTotal\tTemperature\n";
    for(int idx=0; idx<work.accept_timeseries.size(); idx++)
    {
        fileCounters<<work.accept_timeseries[idx]<<"\t"<<work.mcs_timeseries[idx]<<"\t"<<work.counters_t[idx]<<"\n";
    }
    fileCounters.close();
    }
    MPI_Win_free(&win_neighbours);
    MPI_Win_free(&win_distance);
    MPI_Finalize();
}