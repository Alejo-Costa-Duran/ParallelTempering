#ifndef RN_GEN_H
#define RN_GEN_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace rn_gen
{
    extern gsl_rng *rn_mt;

    void initialize_random_number_generator(int rank);
    int rand_site();
    int rand_spin();
    double rand_double();

}
#endif