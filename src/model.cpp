#include "include/model.h"
#include "include/vecinos.h"

model::model(){}
model::model(int rank)
{   
    rn_gen::initialize_random_number_generator(rank);
    nSpins = settings::model::nSpins;
    lattice.resize(nSpins);
    randomize_lattice();
    isPeriodic = settings::model::isPeriodic;
    if(isPeriodic)
    {
        neighboursMatrix = vecinos[settings::model::nClusters-1];
    }
    set_E();
    set_M();
}

void model::randomize_lattice()
{
    for(int idx = 0; idx<nSpins; idx++)
    {
        lattice[idx] = rn_gen::rand_spin();
    }
}

void model::acceptMove(int trialSite)
{
    lattice[trialSite] = -lattice[trialSite];
    E = E_trial;
    M = M_trial;
}

void model::trialMove(int trialSite)
{
    double neighSum =0;
    for(int neighbour : neighboursMatrix[trialSite])
    {
        neighSum += lattice[neighbour];
    }
    delE = -2*lattice[trialSite]*(neighSum);
    E_trial = E+delE;
    M_trial = M-2*lattice[trialSite];
}

void model::set_E()
{
    E = compute_E();
}

void model::set_M()
{

    M = compute_M();
}

double model::compute_E()
{
    double currentEnergy = 0;
    for(int site=0; site<nSpins; site++)
    {
        currentEnergy += -2*lattice[site]*field;
        for(int vec : neighboursMatrix[site])
        {
            currentEnergy += lattice[site]*lattice[vec];
        }
    }
    return currentEnergy/2;
}

int model::compute_M()
{
    int currentM = 0;
    for(int idx = 0; idx<nSpins; idx++)
    {
        currentM += lattice[idx];
    }
    return currentM;
}