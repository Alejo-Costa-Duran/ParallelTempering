#include "include/worker.h"
#include <mpi.h>

worker::worker(int rank, int numWorkers)
{
    world_rank = rank;
    world_size = numWorkers;
}