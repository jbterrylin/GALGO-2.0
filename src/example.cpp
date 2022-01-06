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
        using _TYPE = double;       // Suppport float, double, char, int, long, ... for parameters
        const int NBIT = 64;        // Has to remain between 1 and 64
        bool resultToCsv = false;

        // CONFIG
        galgo::ConfigInfo<_TYPE> config;        // A new instance of config get initial defaults
        set_classic_config<_TYPE>(config);      // Override some defaults
        {
            std::cout << std::endl;
            std::cout << "Rastrigin function";
            // galgo::Parameter<_TYPE, NBIT > par1({ (_TYPE)-4.0,(_TYPE)5.0 });
            // galgo::Parameter<_TYPE, NBIT > par2({ (_TYPE)-4.0,(_TYPE)5.0 });
            // galgo::Parameter<_TYPE, NBIT > par3({ (_TYPE)-4.0,(_TYPE)5.0 });

            config.Objective = ShiftedandRotatedRosenbrockObjective<_TYPE>::Objective;
            
            std::vector<galgo::Parameter<_TYPE, NBIT >> myvector {};
            for (int z = 0; z < 30; z++) myvector.push_back(galgo::Parameter<_TYPE, NBIT > ({ (_TYPE)-100,(_TYPE)100 }));
            
            // std::vector<_TYPE> v;
            // for (int z = 0; z < 3 * config.popsize; z++) v.push_back( (_TYPE) (-4.0 + z *0.01) );
            // galgo::GeneticAlgorithm<_TYPE> my_ga(config, v, par1, par2, par3);
            
            galgo::GeneticAlgorithm<_TYPE> my_ga(config, myvector);
            my_ga.resultToCsv = resultToCsv;
            my_ga.csvFileName = "P1XO + Rastrigin function";

            my_ga.run();
            
            free(y);
            // free(z);
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

