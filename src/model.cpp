#include "include/model.h"

model::model()
{
    nSpins = 10;
    lattice.resize(nSpins);
    randomize_lattice();
}

void model::randomize_lattice()
{
    for(int idx = 0; idx<nSpins; idx++)
    {
        lattice[idx] = rn_gen::rand_spin();
    }
}