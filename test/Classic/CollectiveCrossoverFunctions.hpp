#pragma once

#ifndef COLLECTIVECROSSOVERFUNCTION_HPP
#define COLLECTIVECROSSOVERFUNCTION_HPP

template <typename _TYPE>
void set_config(galgo::ConfigInfo<_TYPE>& config)
{
    // override some defaults
    config.mutinfo._sigma = 1.0;
    config.mutinfo._sigma_lowest = 0;
    config.mutinfo._ratio_boundary = 1;

    config.covrate = 1;  // 0.0 if no cros-over
    config.mutrate = 0.3;
    config.recombination_ratio = 0.50;

    config.elitpop = 0;
    config.tntsize = 4;
    config.Selection = RWS; // TNT; //RWS
    config.CrossOver = CollectiveCrossover; //P1XO
    config.isMultiCrossover = false;
    config.mutinfo._type = galgo::MutationType::MutationGAM_sigma_adapting_per_mutation; //MutationSPM

    config.popsize = 100;
    config.nbgen = 30;
    config.output = false;
}

template <typename _TYPE>
void set_CollectiveCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.nbgen = 30;
    config.Selection = NoSelection;
    config.csvFileName = "CollectiveCrossover";
    config.CrossOver = CollectiveCrossover; //P1XO
    config.isMultiCrossover = true;
    config.mutinfo._type = galgo::MutationType::MutationNoMutation;
    config.mutrate = 0.3;
}

template <typename _TYPE>
void set_SinglePointCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.nbgen = 30;
    config.Selection = RWS;
    config.csvFileName = "P1XO";
    config.CrossOver = P1XO; //P1XO
    config.isMultiCrossover = false;
    config.mutinfo._type = galgo::MutationType::MutationGAM_sigma_adapting_per_mutation;
    config.mutrate = 0.1;
}

template <typename _TYPE>
void set_TwoPointCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.nbgen = 30;
    config.Selection = RWS;
    config.csvFileName = "P2XO";
    config.CrossOver = P2XO; //P1XO
    config.isMultiCrossover = false;
    config.mutinfo._type = galgo::MutationType::MutationGAM_sigma_adapting_per_mutation;
    config.mutrate = 0.1;
}

#endif