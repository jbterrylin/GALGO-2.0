//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#include "Galgo.hpp"

#define TEST_INIT_POP
#ifdef TEST_INIT_POP
#include "../test/Classic/CollectiveCrossoverFunctions.hpp"
#endif

// const int NBIT = 30;        // Has to remain between 1 and 64
std::vector<galgo::Parameter<galgo::_TYPE, galgo::NBIT >> myvector;

template <typename Z>  using FuncKT = std::vector<double>(*)(const std::vector<Z>&);

template <typename _TYPE, typename T>
void runGA(galgo::ConfigInfo<_TYPE> config, FuncKT<T> Objective, std::string benchmarkName, int dimension) {
    bool resultToCsv = true;
    std::cout << std::endl;
    std::cout << config.csvFileName << "->" << benchmarkName << std::endl;

    config.Objective = Objective;
    config.nbgen *= dimension;

    for(int i=0; i < 30; i++) {
        galgo::GeneticAlgorithm<galgo::_TYPE> my_ga(config, myvector);
        my_ga.resultToCsv = resultToCsv;
        my_ga.csvFileName += "-" + std::to_string(dimension) + "+" + benchmarkName;
        my_ga.times = i;
        my_ga.run();
    }
}

template <typename _TYPE>
void runFuntions(galgo::ConfigInfo<_TYPE>& config, int dimension) {
    myvector.clear();
    for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));    
    runGA(config, ShiftedandRotatedBentCigarObjective<galgo::_TYPE>::Objective, "ShiftedandRotatedBentCigarObjective", dimension);

    myvector.clear();
    for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));    
    runGA(config, ShiftedandRotatedZakharovObjective<galgo::_TYPE>::Objective, "ShiftedandRotatedZakharovObjective", dimension);

    myvector.clear();
    for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));    
    runGA(config, ShiftedandRotatedRosenbrockObjective<galgo::_TYPE>::Objective, "ShiftedandRotatedRosenbrockObjective", dimension);

    // myvector.clear();
    // for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));    
    // runGA(config, ShiftedandRotatedRastriginObjective<galgo::_TYPE>::Objective, "ShiftedandRotatedRastriginObjective", dimension);

    // myvector.clear();
    // for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));    
    // runGA(config, ShiftedandRotatedLunacekBi_RastriginObjective<galgo::_TYPE>::Objective, "ShiftedandRotatedLunacekBi_RastriginObjective", dimension);

    // myvector.clear();
    // for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));    
    // runGA(config, ShiftedandRotatedLevyObjective<galgo::_TYPE>::Objective, "ShiftedandRotatedLevyObjective", dimension);

    // myvector.clear();
    // for (int z = 0; z < dimension; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));    
    // runGA(config, ShiftedandRotatedSchwefelObjective<galgo::_TYPE>::Objective, "ShiftedandRotatedSchwefelObjective", dimension);
}

int main()
{
#ifdef TEST_INIT_POP
    // Test init initial population
    {
        // CONFIG
        galgo::ConfigInfo<galgo::_TYPE> config;        // A new instance of config get initial defaults
        set_config<galgo::_TYPE>(config);      // Override some defaults
        std::vector<int> dimension = {30, 50, 100};
        for(int i=0; i< dimension.size(); i++) {
            set_CollectiveCrossover<galgo::_TYPE>(config);
            runFuntions(config, dimension[i]);
            set_SinglePointCrossover<galgo::_TYPE>(config);
            runFuntions(config, dimension[i]);
            set_TwoPointCrossover<galgo::_TYPE>(config);
            runFuntions(config, dimension[i]);
        }
    
        free(y);
        free(M);
        free(OShift);
        free(x_bound);
    }
#endif

#ifdef _WIN32
	system("pause");
#endif

    return 0;
}

