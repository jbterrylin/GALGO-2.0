//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#include "Galgo.hpp"
//#define TEST_ALL_TYPE
//#define TEST_BINAIRO

#ifdef TEST_BINAIRO
#include "..\test\Binairo\GA.h"
#endif

template <typename T>
using CROSS = void (*)(const galgo::Population<T>&, galgo::CHR<T>&, galgo::CHR<T>&);

template <typename T>
using SELECT = void(*)(galgo::Population<T>&);

template <typename T>
using FIXEDPARAM = void(*)(galgo::Population<T>&);

template <typename T>
double Ackley(T x, T y);

// Rosenbrock objective class example
template <typename T>
class RosenbrockObjective
{
public:
   // objective function example : Rosenbrock function
   // minimizing f(x,y) = (1 - x)^2 + 100 * (y - x^2)^2
   static std::vector<double> Objective(const std::vector<T>& x)
   {
        double obj =  -(pow(1-x[0],2)+100*pow(x[1]-x[0]*x[0],2));
        return {obj};
   }
   // NB: GALGO maximize by default so we will maximize -f(x,y)
};

template <typename T>
class AckleyObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        double obj = -Ackley<T>(x[0], x[1]);
        return { obj };
    }
};

// constraints example:
template <typename T>
std::vector<double> MyConstraint(const std::vector<T>& x)
{
   return 
   {
       //x[0]*x[1]+x[0]-x[1]+1.5,   // 1) x * y + x - y + 1.5 <= 0
       //10-x[0]*x[1]               // 2) 10 - x * y <= 0
       x[0] - 2,    // x0 <= 2
       x[1] - 2     // x1 <= 2
   };
}

//
// this example is from
// https://en.wikipedia.org/wiki/Test_functions_for_optimization
// Ackley's function:
// f(x) = -20*exp(-0.2*sqrt(0.5*(x^2 + y^2))) - exp(0.5*(cos(2*pi*x) + cos(2*pi*y))) + exp(1) + 20
// solution is: (0,0)
//
template <typename T>
double Ackley(T x, T y)
{
    const double pi = 3.14159265358979323846;

    double obj_val =  -20 * std::exp(-0.2*std::sqrt(0.5*(x*x + y*y))) - std::exp(0.5*(std::cos(2 * pi*x) + std::cos(2 * pi*y))) + std::exp(1) + 20;
    return obj_val;
}

/*
Rastrigin Function
*/
template <typename T>
double pso_rastrigin(std::vector< T > particle)
{
    double result(10. * static_cast<T> (particle.size())), A(10.), PI(3.14159);
    for (auto dim : particle) {
        result += pow(dim, 2.) - (A * cos(2. * PI * dim));
    }
    return (result);
}

template <typename T>
class rastriginObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        double obj = -pso_rastrigin(x);
        return { obj };
    }
};

/*
Griewank Function
*/
template <typename T>
double pso_griewank(std::vector< T > particle) {
    double sum(0.), product(1.);
    for (int i = 0; i < particle.size(); i++) {
        sum += pow(particle[i], 2.);
        product *= cos(particle[i] / sqrt(i + 1));
    }
    return (1. + (sum / 4000.) - product);
}

template <typename T>
class GriewankObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        double obj = -pso_griewank(x);
        return { obj };
    }
};

/*
Styblinski-Tang Function
Min = (-2.903534,...,--2.903534)
*/
template <typename T>
double pso_styb_tang(std::vector< T > particle)
{
    double result(0.);
    for (auto dim : particle) {
        result += pow(dim, 4.0) - (16. * pow(dim, 2.)) + (5. * dim);
    }
    return (result / 2.);
}

template <typename T>
class StyblinskiTangObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        double obj = -pso_styb_tang(x);
        return { obj };
    }
};

template <typename T>
class SumSameAsPrdObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        int ix = (int)x[0];
        int iy = (int)x[1];
        double sum = ix + iy;
        double prd = ix * iy;
        double diff = std::fabs(sum - prd);

        double err = 1000 * diff * diff;;
        err += (100 * std::fabs(x[0] - ix)* std::fabs(x[0] - ix) + 100 * std::fabs(x[1] - iy)* std::fabs(x[1] - iy));

        double obj = -(diff + err);
        return { obj };
    }
};

//---------------------------------------------------
// TEST templates compiling
// Generate all templates (for a parameter type) to see if compiling/running ok 
//---------------------------------------------------
template <typename _TYPE>
void TEST_TYPE()
{
    galgo::MutationInfo<_TYPE> mutinfo;     // Changes mutation info as desired
    mutinfo._sigma = 1.0;
    mutinfo._sigma_lowest = 0.01;
    mutinfo._ratio_boundary = 0.10;
    mutinfo._type = galgo::MutationType::MutationGAM_UncorrelatedNStepSizeBoundary;

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

        const int       POPUL = 5;
        const int       N = 5;
        const double    MUTRATE = 0.05;
        const int       NBIT = 32;
        const double    RecombinationRatio = 0.60;

        for (size_t s = 0; s < selectcases.size(); s++)
        {
            for (size_t m = 0; m < mutcases.size(); m++)
            {
                mutinfo._type = mutcases[m];
                for (size_t c = 0; c < crosscases.size(); c++)
                {
                    std::cout << std::endl;
                    std::cout << "SumSameAsPrd function 2x2 = 2+2";
                    galgo::Parameter<_TYPE, NBIT> par1({ (_TYPE)1, (_TYPE)100, 60 });
                    galgo::Parameter<_TYPE, NBIT> par2({ (_TYPE)1, (_TYPE)100, 40 });
                    galgo::GeneticAlgorithm<_TYPE> ga(SumSameAsPrdObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2);
                    ga.mutrate = MUTRATE;
                    ga.recombination_ratio = RecombinationRatio;
                    ga.Selection = selectcases[s];
                    ga.CrossOver = crosscases[c];
                    ga.run();
                }
            }
        }
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


int main()
{
#ifdef TEST_ALL_TYPE
    TEST_all_types();
#endif

#ifdef TEST_BINAIRO
    test_ga_binairo();
#endif

    using _TYPE = float;                    // Suppport float, double, char, int, long, ... for parameters
    galgo::MutationInfo<_TYPE> mutinfo;     // Changes mutation info as desired
    mutinfo._sigma          = 1.0;
    mutinfo._sigma_lowest   = 0.01;
    mutinfo._ratio_boundary = 0.10;
    mutinfo._type           = galgo::MutationType::MutationGAM_UncorrelatedNStepSizeBoundary;

    const int       POPUL   = 50;
    const int       N       = 400;  // Number of generation to produce
    const double    MUTRATE = 0.05;
    const int       NBIT    = 63;   // has to remain between 1 and 64
    const double    RecombinationRatio = 0.60;
    const int       TNT_SIZE = 2;
    CROSS<_TYPE>    CROSSType = RealValuedSimpleArithmeticRecombination;
    SELECT<_TYPE>   SELECTType  = TNT;

    {
        {
            std::cout << std::endl;
            std::cout << "SumSameAsPrd function 2x2 = 2+2";
            galgo::Parameter<_TYPE, NBIT> par1({ (_TYPE)1, (_TYPE)100, 50 }); // an initial value can be added inside the initializer list after the upper bound
            galgo::Parameter<_TYPE, NBIT> par2({ (_TYPE)1, (_TYPE)100, 50 });
            galgo::GeneticAlgorithm<_TYPE> ga(SumSameAsPrdObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2);
            ga.mutrate = MUTRATE;  
            ga.recombination_ratio = RecombinationRatio;
            ga.tntsize = TNT_SIZE;
            ga.Selection = SELECTType;
            ga.CrossOver = CROSSType;
            ga.run();
        }

        {
            std::cout << std::endl;
            std::cout << "Rosenbrock function";
            galgo::Parameter<_TYPE, NBIT > par1({ (_TYPE)-2.0,(_TYPE)2.0 });
            galgo::Parameter<_TYPE, NBIT > par2({ (_TYPE)-2.0,(_TYPE)2.0 });
            galgo::GeneticAlgorithm<_TYPE> ga(RosenbrockObjective< _TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2);
            ga.mutrate = MUTRATE;
            ga.recombination_ratio = RecombinationRatio;
            ga.tntsize = TNT_SIZE;
            ga.Selection = SELECTType;
            ga.CrossOver = CROSSType;
            ga.run();
        }

        {
            std::cout << std::endl;
            std::cout << "Ackley function";
            galgo::Parameter<_TYPE, NBIT > par1({ (_TYPE)-4.0,(_TYPE)5.0 });
            galgo::Parameter<_TYPE, NBIT > par2({ (_TYPE)-4.0,(_TYPE)5.0 });
            galgo::GeneticAlgorithm<_TYPE> ga(AckleyObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2);
            ga.mutrate = MUTRATE;
            ga.recombination_ratio = RecombinationRatio;
            ga.tntsize = TNT_SIZE;
            ga.Selection = SELECTType;
            ga.CrossOver = CROSSType;
            ga.run();
        }

        {
            std::cout << std::endl;
            std::cout << "Rastrigin function";
            galgo::Parameter<_TYPE, NBIT > par1({ (_TYPE)-4.0,(_TYPE)5.0 });
            galgo::Parameter<_TYPE, NBIT > par2({ (_TYPE)-4.0,(_TYPE)5.0 });
            galgo::Parameter<_TYPE, NBIT > par3({ (_TYPE)-4.0,(_TYPE)5.0 });
            galgo::GeneticAlgorithm<_TYPE> ga(rastriginObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2, par3);
            ga.mutrate = MUTRATE;
            ga.recombination_ratio = RecombinationRatio;
            ga.tntsize = TNT_SIZE;
            ga.Selection = SELECTType;
            ga.CrossOver = CROSSType;
            ga.run();
        }

        {
            std::cout << std::endl;
            std::cout << "StyblinskiTang function Min = (-2.903534,...,--2.903534)";
            galgo::Parameter<_TYPE, NBIT > par1({ (_TYPE)-4.0,(_TYPE)4.0 });
            galgo::Parameter<_TYPE, NBIT > par2({ (_TYPE)-4.0,(_TYPE)4.0 });
            galgo::Parameter<_TYPE, NBIT > par3({ (_TYPE)-4.0,(_TYPE)4.0 });
            galgo::GeneticAlgorithm<_TYPE> ga(StyblinskiTangObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2, par3);
            ga.mutrate = MUTRATE;
            ga.recombination_ratio = RecombinationRatio;
            ga.tntsize = TNT_SIZE;
            ga.Selection = SELECTType;
            ga.CrossOver = CROSSType;
            ga.run();
        }

        {
            std::cout << std::endl;
            std::cout << "Griewank function";
            galgo::Parameter<_TYPE, NBIT > par1({ (_TYPE)-4.0,(_TYPE)5.0 });
            galgo::Parameter<_TYPE, NBIT > par2({ (_TYPE)-4.0,(_TYPE)5.0 });
            galgo::Parameter<_TYPE, NBIT > par3({ (_TYPE)-4.0,(_TYPE)5.0 });
            galgo::GeneticAlgorithm<_TYPE> ga(GriewankObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2, par3);
            ga.mutrate = MUTRATE;
            ga.recombination_ratio = RecombinationRatio;
            ga.tntsize = TNT_SIZE;
            ga.Selection = SELECTType;
            ga.CrossOver = CROSSType;
            ga.run();
        }
    }

    system("pause");
}

