#include "include/rngenerator.h"

namespace rn_gen
{
    gsl_rng *rn_mt;
    void initialize_random_number_generator(int rank) {
    // Use the rank to create a unique seed for each process
    unsigned long seed = rank+1;
    rn_mt = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rn_mt, seed);
}
    int rand_site()
    {
        return gsl_rng_uniform_int(rn_mt,24);
    }
    int rand_spin()
    {
        return gsl_rng_uniform_int(rn_mt,2)*2-1;
    }
    double rand_double()
    {
        return gsl_rng_uniform(rn_mt);
    }
}