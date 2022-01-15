//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#include "Galgo.hpp"

#define TEST_INIT_POP
#ifdef TEST_INIT_POP
#include "../test/Classic/RingCrossoverFunctions.hpp"
#endif

// const int NBIT = 30;        // Has to remain between 1 and 64
std::vector<galgo::Parameter<galgo::_TYPE, galgo::NBIT >> myvector;

template <typename Z>  using FuncKT = std::vector<double>(*)(const std::vector<Z>&);

template <typename _TYPE, typename T>
void runGA(galgo::ConfigInfo<_TYPE>& config, FuncKT<T> Objective, std::string benchmarkName) {
    bool resultToCsv = true;
    std::cout << std::endl;
    std::cout << config.csvFileName << "->" << benchmarkName;

    config.Objective = Objective;

    for(int i=0; i < 30; i++) {
        galgo::GeneticAlgorithm<galgo::_TYPE> my_ga(config, myvector);
        my_ga.resultToCsv = resultToCsv;
        my_ga.csvFileName += "+" + benchmarkName;
        my_ga.times = i;
        my_ga.run();
    }
}

template <typename _TYPE>
void runFuntions(galgo::ConfigInfo<_TYPE>& config, int dimension) {
    myvector.clear();
    for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-5.12, (galgo::_TYPE)5.12 }));    
    runGA(config, SphereObjective<galgo::_TYPE>::Objective, "SphereObjective");

    myvector.clear();
    for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-5.12, (galgo::_TYPE)5.12 }));    
    runGA(config, AxisParallelHyperEllipsoidObjective<galgo::_TYPE>::Objective, "AxisParallelHyperEllipsoidObjective");

    myvector.clear();
    for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-65.536, (galgo::_TYPE)65.536 }));    
    runGA(config, RotatedHyperEllipsoidObjective<galgo::_TYPE>::Objective, "RotatedHyperEllipsoidObjective");

    myvector.clear();
    for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-500, (galgo::_TYPE)500 }));    
    runGA(config, NormalizedSchwefelObjective<galgo::_TYPE>::Objective, "NormalizedSchwefelObjective");

    myvector.clear();
    for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-5.12, (galgo::_TYPE)5.12 }));    
    runGA(config, GeneralizedRastriginObjective<galgo::_TYPE>::Objective, "GeneralizedRastriginObjective");

    myvector.clear();
    for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-2.048, (galgo::_TYPE)2.048 }));    
    runGA(config, RosenbrocksValleyObjective<galgo::_TYPE>::Objective, "RosenbrocksValleyObjective");
}

int main()
{
#ifdef TEST_INIT_POP
    // Test init initial population
    {
        // CONFIG
        galgo::ConfigInfo<galgo::_TYPE> config;        // A new instance of config get initial defaults
        set_config<galgo::_TYPE>(config);      // Override some defaults
        set_RingCrossover<galgo::_TYPE>(config);
        runFuntions(config, 30);
        set_SinglePointCrossover<galgo::_TYPE>(config);
        runFuntions(config, 30);
        set_TwoPointCrossover<galgo::_TYPE>(config);
        runFuntions(config, 30);
        set_HeuristicCrossover<galgo::_TYPE>(config);
        runFuntions(config, 30);
        set_IntermediateCrossover<galgo::_TYPE>(config);
        runFuntions(config, 30);
        set_ArithmeticCrossover<galgo::_TYPE>(config);
        runFuntions(config, 30);
    }
#endif

#ifdef _WIN32
	system("pause");
#endif

    return 0;
}

