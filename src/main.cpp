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
    std::cout<<"Finished computing distancesn\n";
    return distance_matrix;
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
    // Define variables for the shared memory window
    int* shared_neighbours = nullptr;
    int* shared_distanceMatrix = nullptr;
    MPI_Win win_neighbours, win_distance;

   // Rank 0 in each node will allocate memory for the neighbours and distance matrices
    int nSpins = settings::model::nSpins;
    int neighbours_size = nSpins * 7;
    int distanceMatrix_size = (nSpins * (nSpins + 1)) / 2;

    // Allocate shared memory window for neighbors matrix
    if (local_rank == 0) {
        // Rank 0 on each node allocates memory
        MPI_Win_allocate_shared(neighbours_size * sizeof(int), sizeof(int), MPI_INFO_NULL,
                                node_comm, &shared_neighbours, &win_neighbours);

        if(settings::sim::storeCorrelations){       
        MPI_Win_allocate_shared(distanceMatrix_size * sizeof(int), sizeof(int), MPI_INFO_NULL,
                                node_comm, &shared_distanceMatrix, &win_distance);
        }
        // Only rank 0 loads the data
        std::ifstream file("./src/Vecinos"+std::to_string(settings::model::neighIdx)+".csv");
        int number;
        std::string line;
        int idx = 0;
        while (std::getline(file, line)) {
                number = std::stoi(line);
                shared_neighbours[idx] = number;
                idx++;
            }
        file.close();
        
        if(settings::sim::storeCorrelations){
        // Compute distance matrix
        std::vector<int> distanceMatrix = computeDistanceMatrix(std::vector<int>(shared_neighbours, shared_neighbours + neighbours_size), nSpins);
        std::copy(distanceMatrix.begin(), distanceMatrix.end(), shared_distanceMatrix);
        std::cout << "Rank 0 on node " << processor_name << " initialized the neighbours and distance matrices.\n";
        }
        std::cout << "Rank 0 on node " << processor_name << " initialized the neighbours matrix.\n";
        MPI_Win_fence(0, win_neighbours);
    } else {
        // Other ranks on the node attach to the shared memory
        MPI_Aint size;
        int disp_unit;
        int mpi_err;
        // Attach to the shared memory window for neighbors matrix
        mpi_err = MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, node_comm, &shared_neighbours, &win_neighbours);
        if (mpi_err != MPI_SUCCESS) {
            std::cout<<mpi_err<<"Something wrong when sharing memory\n";
        }
        
        MPI_Win_shared_query(win_neighbours, 0, &size, &disp_unit, &shared_neighbours);
        if(settings::sim::storeCorrelations){
        // Attach to the shared memory window for distance matrix
        MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, node_comm, &shared_distanceMatrix, &win_distance);
        MPI_Win_shared_query(win_distance, 0, &size, &disp_unit, &shared_distanceMatrix);
        }
        MPI_Win_fence(0, win_neighbours);
    }

    // Synchronize after data load
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout<<"Rank "<<world_rank<<" on node "<<processor_name<<" has successfully loaded the data.\n";
    worker work(world_rank,world_size, shared_neighbours);

    for(int iteration = 0; iteration<settings::sim::MCS_total; iteration++)
        {           
            if(!work.thermalized){work.thermalization(work.T,shared_neighbours);};
            if(iteration%settings::sim::MCS_swap == 0){work.swap_workers();}
            if(work.world_rank==1){if(iteration%100000 == 0){std::cout<<"Warmup iteration number: "<<iteration<<"          \n";}}
        }

        for(int iteration = 0; iteration<settings::sim::MCS_total; iteration++)
        {
            if(!work.thermalized){work.thermalization(work.T,shared_neighbours);};
            if(iteration%settings::sim::MCS_swap == 0){work.swap_workers();}
            work.sampling(shared_neighbours,shared_distanceMatrix,settings::sim::storeCorrelations,settings::sim::ladderUpdate);
            if(work.world_rank==1){if(iteration%100000 == 0){std::cout<<"Sampling iteration number: "<<iteration<<"        \n";}}
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
    if(settings::sim::storeCorrelations)
    {
        MPI_Win_free(&win_distance);
    }
    MPI_Finalize();
}