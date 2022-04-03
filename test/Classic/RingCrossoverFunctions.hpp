#pragma once

#ifndef RINGCROSSOVERFUNCTION_HPP
#define RINGCROSSOVERFUNCTION_HPP

template <typename _TYPE>
void set_config(galgo::ConfigInfo<_TYPE>& config)
{
    // override some defaults
    config.mutinfo._sigma = 1.0;
    config.mutinfo._sigma_lowest = 0.01;
    config.mutinfo._ratio_boundary = 0.25;

    config.covrate = 0.8;  // 0.0 if no cros-over
    config.mutrate = 0.01;
    config.recombination_ratio = 0.5;

    config.elitpop = 0;
    config.tntsize = 4;
    config.Selection = SUS; // TNT; //RWS
    config.CrossOver = RingCrossover; //P1XO
    config.isMultiCrossover = false;
    config.mutinfo._type = galgo::MutationType::MutationGAM_sigma_adapting_per_mutation; //MutationSPM

    config.popsize = 20;
    config.nbgen = 10000 / config.popsize;
    config.output = false;
}

template <typename _TYPE>
void set_RingCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "RingCrossover";
    config.CrossOver = RingCrossover; //P1XO
}

template <typename _TYPE>
void set_SinglePointCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "P1XO";
    config.CrossOver = P1XO; //P1XO
}

template <typename _TYPE>
void set_TwoPointCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "P2XO";
    config.CrossOver = P2XO; //P1XO
}

template <typename _TYPE>
void set_HeuristicCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "HeuristicCrossover";
    config.CrossOver = HeuristicCrossover; //P1XO
}

template <typename _TYPE>
void set_IntermediateCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "IntermediateCrossover";
    config.CrossOver = IntermediateCrossover; //P1XO
}

template <typename _TYPE>
void set_ArithmeticCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "ArithmeticCrossover";
    config.CrossOver = RealValuedWholeArithmeticRecombination; //P1XO
}

#endif