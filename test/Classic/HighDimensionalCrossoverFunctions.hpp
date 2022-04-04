#pragma once

#ifndef COLLECTIVECROSSOVERFUNCTION_HPP
#define COLLECTIVECROSSOVERFUNCTION_HPP

template <typename _TYPE>
void set_config(galgo::ConfigInfo<_TYPE>& config)
{
    // override some defaults
    config.mutinfo._sigma = 1.0;
    config.mutinfo._sigma_lowest = 0.01;
    config.mutinfo._ratio_boundary = 1;

    config.covrate = 1;  // 0.0 if no cros-over
    config.mutrate = 0.7 / galgo::HDGA_NBIT;
    config.recombination_ratio = 0.50;

    config.tntsize = 4;
    config.Selection = RWS; // TNT; //RWS
    config.CrossOver = HighDimensionalGeneticAlgorithmToolboxCrossover; //P1XO
    config.isMultiCrossover = false;
    config.mutinfo._type = galgo::MutationType::MutationSPM; //MutationSPM

    config.popsize = 50;
    
    config.elitpop = 0.05 * config.popsize;
    config.nbgen = 50;
    config.output = false;
}

template <typename _TYPE>
void set_HighDimensionalGeneticAlgorithmToolboxCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "HighDimensionalGeneticAlgorithmToolboxCrossover";
    config.CrossOver = HighDimensionalGeneticAlgorithmToolboxCrossover; //P1XO
    config.isMultiCrossover = true;
}

template <typename _TYPE>
void set_SinglePointCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "P1XO";
    config.CrossOver = P1XO; //P1XO
    config.isMultiCrossover = false;
}

#endif