//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#include "Galgo.hpp"

//------------------------------------------------------------------------------
// Uncomment #define TEST_ALL_TYPE to test compiling of all types
// Uncomment #define TEST_BINAIRO to test GA for Binairos
// Uncomment #define TEST_CLASSIC_FUNCTIONS to test GA for classics functions
// Uncomment #define TEST_INIT_POP to test by initializing initial population
//------------------------------------------------------------------------------
//#define TEST_ALL_TYPE
//#define TEST_BINAIRO
//#define TEST_CLASSIC_FUNCTIONS
#define TEST_INIT_POP

#ifdef TEST_CLASSIC_FUNCTIONS
#include "../test/Classic/Functions.hpp"
#endif
#ifdef TEST_INIT_POP
#include "../test/Classic/Functions.hpp"
#endif

#ifdef TEST_ALL_TYPE
#include "../test/Types/TestTypes.hpp"
#endif

#ifdef TEST_BINAIRO
#include "../test/Binairo/GA.h"
#endif

std::vector<galgo::Parameter<galgo::_TYPE, galgo::NBIT >> myvector;

int main()
{
#ifdef TEST_ALL_TYPE
    TEST_all_types();
#endif

#ifdef TEST_BINAIRO
    GA_Binairo::test_ga_binairo(4);     // 0=resolve one free cell(hard), 1=resolve 4 free cells(very hard), 2=resolve 7 free cells(diabolical), 3 , 4==generate new grid
#endif

#ifdef TEST_CLASSIC_FUNCTIONS
    test_classic();
#endif

#ifdef TEST_INIT_POP
    // Test init initial population
    {
        bool resultToCsv = false;
        // CONFIG
        galgo::ConfigInfo<galgo::_TYPE> config;        // A new instance of config get initial defaults
        set_classic_config<galgo::_TYPE>(config);      // Override some defaults
        {
            std::cout << std::endl;
            std::cout << "Rastrigin function";

            config.Objective = OneMaxObjective<galgo::_TYPE>::Objective;
            
            for (int z = 0; z < 5; z++) myvector.push_back(galgo::Parameter<galgo::_TYPE, galgo::NBIT > ({ (galgo::_TYPE)-100, (galgo::_TYPE)100 }));
            
            galgo::GeneticAlgorithm<galgo::_TYPE> my_ga(config, myvector);
            my_ga.resultToCsv = resultToCsv;
            my_ga.csvFileName = "P1XO + Rastrigin function";

            my_ga.run();
            
            free(y);
            free(M);
            free(OShift);
            free(x_bound);
        }
    }
#endif

#ifdef _WIN32
	system("pause");
#endif

    return 0;
}

