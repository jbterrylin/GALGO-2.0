#pragma once

#ifndef HYBRIDCROSSOVERFUNCTION_HPP
#define HYBRIDCROSSOVERFUNCTION_HPP

template <typename _TYPE>
void set_config(galgo::ConfigInfo<_TYPE>& config)
{
    // override some defaults
    config.mutinfo._sigma = 1.0;
    config.mutinfo._sigma_lowest = 0.01;
    config.mutinfo._ratio_boundary = 0.10;

    config.covrate = 0.85;  // 0.0 if no cros-over
    config.mutrate = 0.05;
    config.recombination_ratio = 0.50;

    config.elitpop = 1;
    config.tntsize = 4;
    config.Selection = RWS; // TNT; //RWS
    config.CrossOver = HybridCrossover; //P1XO
    config.isMultiCrossover = false;
    config.mutinfo._type = galgo::MutationType::MutationSPM; //MutationSPM

    config.popsize = 50;
    config.nbgen = 500;
    config.output = false;
}

template <typename _TYPE>
void set_HybridCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "HybridCrossover";
    config.CrossOver = HybridCrossover; //P1XO
    config.isMultiCrossover = true;
}

template <typename _TYPE>
void set_SinglePointCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "P1XO";
    config.CrossOver = P1XO; //P1XO
    config.isMultiCrossover = false;
}

template <typename _TYPE>
void set_TwoPointCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "P2XO";
    config.CrossOver = P2XO; //P1XO
    config.isMultiCrossover = false;
}

template <typename _TYPE>
void set_UniformCrossover(galgo::ConfigInfo<_TYPE>& config)
{
    config.csvFileName = "UniformCrossover";
    config.CrossOver = UXO; //P1XO
    config.isMultiCrossover = false;
}

#endif