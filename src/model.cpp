#include "include/model.h"
#include <fstream>
#include <iostream>
#include <string>

model::model(){}
model::model(int rank,int *shared_neighbours)
:isPeriodic(settings::model::isPeriodic),
nSpins(settings::model::nSpins),
rank(rank),
lattice(new int[settings::model::nSpins]),
E_trial(10000),
M_trial(10000),
delE(0),
field(settings::model::field)
{   
    rn_gen::initialize_random_number_generator(rank+12);
    randomize_lattice();
    set_E(shared_neighbours);
    set_M();
}


std::vector<double> model::compute_Correlations(int *shared_distanceMatrix)
{
    std::vector<double> correlations(settings::model::distances+1,0);
    std::vector<int> counts(settings::model::distances+1,0);
    counts[0] = nSpins;
    for(int currentSpin = 0; currentSpin<settings::model::nSpins;currentSpin++)
    for(int currentNeigh = currentSpin+1; currentNeigh<settings::model::nSpins; currentNeigh++)
    {
        int distanceMatrixIdx = currentNeigh*(currentNeigh+1)/2+currentSpin;
        int distance = shared_distanceMatrix[distanceMatrixIdx];
        correlations[distance] += lattice[currentNeigh]*lattice[currentSpin];
        counts[distance]+=1;
    }
    for(int correlationsIdx = 0; correlationsIdx<settings::model::distances+1; correlationsIdx++ )
    {
        correlations[correlationsIdx] = correlations[correlationsIdx]/counts[correlationsIdx];
    }
    return correlations;
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

void model::trialMove(int trialSite,int *shared_neighbours)
{
    int neighSum = 0;
    for(int vec = 0; vec<7; vec++)
    {
        int vecino = shared_neighbours[7*trialSite+vec];
        if(vecino != -1)
        {
        neighSum += lattice[vecino];
        }
    }
    delE = -2*lattice[trialSite]*neighSum+2*field*lattice[trialSite];
    E_trial = E+delE;
    M_trial = M-2*lattice[trialSite];
}

void model::set_E(int *shared_neighbours)
{
    E = compute_E(shared_neighbours);
}

void model::set_M()
{

    M = compute_M();
}

double model::compute_E(int *shared_neighbours)
{
    double currentEnergy = 0;
    for(int site=0; site<nSpins; site++)
    {
        currentEnergy += -2*lattice[site]*field;
        for(int vec = 0; vec<7; vec++)
        {
            int vecino = shared_neighbours[site*7+vec];
            if(vecino != -1)
            {
            currentEnergy += lattice[site]*lattice[vecino];
            }
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