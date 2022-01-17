//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#include "Galgo.hpp"

#define TEST_INIT_POP
#ifdef TEST_INIT_POP
#include "../test/Classic/HybridCrossoverFunctions.hpp"
#endif

std::vector<galgo::Parameter<galgo::_TYPE, galgo::NBIT >> myvector;

template <typename Z>  using FuncKT = std::vector<double>(*)(const std::vector<Z>&);

template <typename _TYPE, typename T>
void runGA(galgo::ConfigInfo<_TYPE>& config, FuncKT<T> Objective, std::string benchmarkName) {
    bool resultToCsv = true;
    std::cout << std::endl;
    std::cout << config.csvFileName << "->" << benchmarkName << std::endl;

    config.Objective = Objective;

    for(int i=0; i < 100; i++) {
        galgo::GeneticAlgorithm<galgo::_TYPE> my_ga(config, myvector);
        my_ga.resultToCsv = resultToCsv;
        my_ga.csvFileName += "+" + benchmarkName;
        my_ga.times = i;
        my_ga.run();
    }
}

template <typename _TYPE>
void runFuntions(galgo::ConfigInfo<_TYPE>& config) {
    myvector.clear();
    for (int z = 0; z < 2; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-5.12, (galgo::_TYPE)5.12 }));
    runGA(config, SphereObjective<galgo::_TYPE>::Objective, "SphereObjective");

    myvector.clear();
    for (int z = 0; z < 2; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-5.12, (galgo::_TYPE)5.12 }));
    runGA(config, GeneralizedRastriginObjective<galgo::_TYPE>::Objective, "GeneralizedRastriginObjective");
    
    myvector.clear();
    for (int z = 0; z < 2; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));
    runGA(config, SchaffersF6Objective<galgo::_TYPE>::Objective, "SchaffersF6Objective");
    
    myvector.clear();
    for (int z = 0; z < 2; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-600, (galgo::_TYPE)600 })); 
    runGA(config, GriewangksObjective<galgo::_TYPE>::Objective, "GriewangksObjective");
    
    myvector.clear();
    for (int z = 0; z < 2; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));
    runGA(config, HansenObjective<galgo::_TYPE>::Objective, "HansenObjective");
    
    myvector.clear();
    for (int z = 0; z < 2; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)0, (galgo::_TYPE)PI }));
    runGA(config, MichalewiczObjective<galgo::_TYPE>::Objective, "MichalewiczObjective");
}

int main()
{
#ifdef TEST_INIT_POP
    // Test init initial population
    {
        // CONFIG
        galgo::ConfigInfo<galgo::_TYPE> config;        // A new instance of config get initial defaults
        set_config<galgo::_TYPE>(config);      // Override some defaults
        set_HybridCrossover<galgo::_TYPE>(config);
        runFuntions(config);
        set_SinglePointCrossover<galgo::_TYPE>(config);
        runFuntions(config);
        set_TwoPointCrossover<galgo::_TYPE>(config);
        runFuntions(config);
    }
#endif

#ifdef _WIN32
	system("pause");
#endif

    return 0;
}

