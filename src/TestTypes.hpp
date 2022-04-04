//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

template <typename T>
using CROSS = void (*)(const galgo::Population<T>&, galgo::CHR<T>&, galgo::CHR<T>&);

template <typename T>
using SELECT = void(*)(galgo::Population<T>&);

template <typename _TYPE>
void set_test_config(galgo::ConfigInfo<_TYPE>& config)
{
    // override some defaults
    config.mutinfo._sigma = 1.0;
    config.mutinfo._sigma_lowest = 0.01;
    config.mutinfo._ratio_boundary = 0.10;

    config.mutrate = 0.05;
    config.recombination_ratio = 0,60;

    config.tntsize = 2;
    config.Selection = TNT;
    config.CrossOver = RealValuedSimpleArithmeticRecombination;
    config.mutinfo._type = galgo::MutationType::MutationGAM_UncorrelatedNStepSizeBoundary;

    config.popsize = 5;
    config.nbgen = 2;
    config.output = true;
}

//---------------------------------------------------
// TEST templates compiling
// Generate all templates (for a parameter type) to see if compiling/running ok 
//---------------------------------------------------
template <typename _TYPE>
void TEST_TYPE()
{
    galgo::ConfigInfo<_TYPE> config;
    set_test_config<_TYPE>(config);

    {
        std::vector<galgo::MutationType> mutcases = {
            galgo::MutationType::MutationSPM,
            galgo::MutationType::MutationBDM,
            galgo::MutationType::MutationUNM,
            galgo::MutationType::MutationGAM_UncorrelatedOneStepSizeFixed,
            galgo::MutationType::MutationGAM_UncorrelatedOneStepSizeBoundary,
            galgo::MutationType::MutationGAM_UncorrelatedNStepSize,
            galgo::MutationType::MutationGAM_UncorrelatedNStepSizeBoundary,
            galgo::MutationType::MutationGAM_sigma_adapting_per_generation,
            galgo::MutationType::MutationGAM_sigma_adapting_per_mutation
        };

        std::vector<CROSS<_TYPE>> crosscases = {
            P1XO<_TYPE>,
            P2XO<_TYPE>,
            UXO<_TYPE>,
            RealValuedSimpleArithmeticRecombination<_TYPE>,
            RealValuedSingleArithmeticRecombination<_TYPE>,
            RealValuedWholeArithmeticRecombination<_TYPE>
        };

        std::vector<SELECT<_TYPE>> selectcases = {
            RWS<_TYPE>,
            SUS<_TYPE>,
            RNK<_TYPE>,
            RSP<_TYPE>,
            TNT<_TYPE>,
            TRS<_TYPE>
        };

        const int NBIT = 32;

        for (size_t s = 0; s < selectcases.size(); s++)
        {
            for (size_t m = 0; m < mutcases.size(); m++)
            {
                config.mutinfo._type = mutcases[m];
                for (size_t c = 0; c < crosscases.size(); c++)
                {
                    std::cout << s << " - " << m << " - " << c << std::endl;
                    std::cout << "SumSameAsPrd function 2x2 = 2+2";
                    galgo::Parameter<_TYPE, NBIT> par1({ (_TYPE)1.5, (_TYPE)3, 3 }); // an initial value can be added inside the initializer list after the upper bound
                    galgo::Parameter<_TYPE, NBIT> par2({ (_TYPE)1.5, (_TYPE)3, 3 });

                    config.Objective = SumSameAsPrdObjective<_TYPE>::Objective;
                    config.Selection = selectcases[s];
                    config.CrossOver = crosscases[c];

                    galgo::GeneticAlgorithm<_TYPE> ga(config, par1, par2);
                    ga.run();
                }
            }
        }
        std::cout << " Done"  << std::endl;
    }
}

void TEST_all_types()
{
    TEST_TYPE<double>();
    TEST_TYPE<float>();
    TEST_TYPE<char>();
    TEST_TYPE<short>();
    TEST_TYPE<int>();
    TEST_TYPE<unsigned char>();
    TEST_TYPE<unsigned short>();
    TEST_TYPE<unsigned int>();
    TEST_TYPE<long>();
    TEST_TYPE<long long>();
    TEST_TYPE<unsigned long>();
    TEST_TYPE<unsigned long long>();
}
