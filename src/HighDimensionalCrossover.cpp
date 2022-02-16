//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#include "Galgo.hpp"

#define TEST_INIT_POP
#ifdef TEST_INIT_POP
#include "../test/Classic/HighDimensionalCrossoverFunctions.hpp"
#endif

std::vector<galgo::Parameter<galgo::_TYPE, galgo::HDGA_NBIT >> HDGAvector;

template <typename Z>  using FuncKT = std::vector<double>(*)(const std::vector<Z>&);

template <typename _TYPE, typename T>
void runGA(galgo::ConfigInfo<_TYPE>& config, FuncKT<T> Objective, std::string benchmarkName) {
    bool resultToCsv = true;
    std::cout << std::endl;
    std::cout << config.csvFileName << "->" << benchmarkName << std::endl;

    config.Objective = Objective;

    for(int i=0; i < 100; i++) {
        galgo::GeneticAlgorithm<galgo::_TYPE> my_ga(config, HDGAvector);
        my_ga.resultToCsv = resultToCsv;
        my_ga.csvFileName += "+" + benchmarkName;
        my_ga.times = i;
        my_ga.run();
    }
}

template <typename _TYPE>
void runFuntions(galgo::ConfigInfo<_TYPE>& config) {
    HDGAvector.clear();
    for (int z = 0; z < 2; z++) HDGAvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::HDGA_NBIT > ({ (galgo::_TYPE)-10, (galgo::_TYPE)10 }));    
    runGA(config, ShubertObjective<galgo::_TYPE>::Objective, "ShubertObjective");

    HDGAvector.clear();
    for (int z = 0; z < 10; z++) HDGAvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::HDGA_NBIT > ({ (galgo::_TYPE)-5.12, (galgo::_TYPE)5.12 }));    
    runGA(config, SphereObjective<galgo::_TYPE>::Objective, "SphereObjective");

    HDGAvector.clear();
    for (int z = 0; z < 10; z++) HDGAvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::HDGA_NBIT > ({ (galgo::_TYPE)-5, (galgo::_TYPE)10 }));    
    runGA(config, ZakharovObjective<galgo::_TYPE>::Objective, "ZakharovObjective");

    HDGAvector.clear();
    for (int z = 0; z < 10; z++) HDGAvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::HDGA_NBIT > ({ (galgo::_TYPE)-10, (galgo::_TYPE)10 }));    
    runGA(config, DixonPriceObjective<galgo::_TYPE>::Objective, "DixonPriceObjective");
}

int main()
{
#ifdef TEST_INIT_POP
    // Test init initial population
    {
        // CONFIG
        galgo::ConfigInfo<galgo::_TYPE> config;        // A new instance of config get initial defaults
        set_config<galgo::_TYPE>(config);      // Override some defaults
        set_HighDimensionalGeneticAlgorithmToolboxCrossover<galgo::_TYPE>(config);
        runFuntions(config);
        set_SinglePointCrossover<galgo::_TYPE>(config);
        runFuntions(config);
    }
#endif

#ifdef _WIN32
	system("pause");
#endif

    return 0;
}

