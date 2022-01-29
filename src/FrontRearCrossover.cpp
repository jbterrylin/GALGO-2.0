//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#include "Galgo.hpp"

#define TEST_INIT_POP
#ifdef TEST_INIT_POP
#include "../test/Classic/FrontRearCrossoverFunctions.hpp"
#endif

std::vector<galgo::Parameter<galgo::_TYPE, galgo::NBIT >> myvector;

template <typename Z>  using FuncKT = std::vector<double>(*)(const std::vector<Z>&);

template <typename _TYPE, typename T>
void runGA(galgo::ConfigInfo<_TYPE>& config, FuncKT<T> Objective, std::string benchmarkName, int dimension) {
    bool resultToCsv = true;
    std::cout << std::endl;
    std::cout << config.csvFileName << "->" << benchmarkName << std::endl;

    if(benchmarkName == "TrapFiveObjective") {
        config.nbgen = 1000;
    } else {
        config.nbgen = 500;
    }

    if(benchmarkName == "TrapFiveObjective") {
        config.popsize = 80;
    } else {
        config.popsize = 30;
    }

    config.Objective = Objective;

    for(int i=0; i < 100; i++) {
        myvector.clear();
        for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-(pow(2,galgo::NBIT-1)-1), (galgo::_TYPE)(pow(2,galgo::NBIT-1)-1) }));
        
        galgo::GeneticAlgorithm<galgo::_TYPE> my_ga(config, myvector);
        my_ga.resultToCsv = resultToCsv;
        my_ga.csvFileName += std::string("-") + std::to_string(dimension) + "+" + benchmarkName;
        my_ga.times = i;
        my_ga.showBits = true;
        my_ga.showNFitnessEvaluateReachTarget = true;
        my_ga.targetFitness = galgo::NBIT * dimension;
        my_ga.run();
    }
}

template <typename _TYPE>
void runFuntions(galgo::ConfigInfo<_TYPE>& config, int dimension) {
    runGA(config, OneMaxObjective<galgo::_TYPE>::Objective, "OneMaxObjective", dimension);
    runGA(config, ZeroMaxObjective<galgo::_TYPE>::Objective, "ZeroMaxObjective", dimension);
    runGA(config, TrapThreeObjective<galgo::_TYPE>::Objective, "TrapThreeObjective", dimension);
    runGA(config, TrapFiveObjective<galgo::_TYPE>::Objective, "TrapFiveObjective", dimension);
}

int main()
{
#ifdef TEST_INIT_POP
    // Test init initial population
    {
        // CONFIG
        galgo::ConfigInfo<galgo::_TYPE> config;        // A new instance of config get initial defaults
        set_config<galgo::_TYPE>(config);      // Override some defaults
        if( 30 % galgo::NBIT != 0 || 
            60 % galgo::NBIT != 0 || 
            120 % galgo::NBIT != 0) {
            std::cout << "NBIT is not suitable to make chromosome length become 30, 60 or 120!" << std::endl;
        } 
        std::vector<int> dimension = {30/galgo::NBIT, 60/galgo::NBIT, 120/galgo::NBIT};
        for(int i=0; i< dimension.size(); i++) {
            set_FrontRearCrossover<galgo::_TYPE>(config);   
            runFuntions(config, dimension[i]);
            set_SinglePointCrossover<galgo::_TYPE>(config);
            runFuntions(config, dimension[i]);
            set_TwoPointCrossover<galgo::_TYPE>(config);
            runFuntions(config, dimension[i]);
            set_UniformCrossover<galgo::_TYPE>(config);
            runFuntions(config, dimension[i]);
            set_RingCrossover<galgo::_TYPE>(config);
            runFuntions(config, dimension[i]);
        }
    }
#endif

#ifdef _WIN32
	system("pause");
#endif

    return 0;
}

