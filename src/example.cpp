//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#include "Galgo.hpp"

//#define TEST_BINAIRO
#ifdef TEST_BINAIRO
#include "..\test\Binairo\MatriceUtil.h"
#include "..\test\Binairo\Algorithm.h"

using BINAIRO_TEST_TYPE = int;
// 4x4
//static std::vector<BINAIRO_TEST_TYPE> binairo_initial  = { 1,0,-1,0, -1,-1,1,-1, -1,-1,0,0, -1,-1,-1,-1 };
//static std::vector<BINAIRO_TEST_TYPE> binairo_solution = { 1,0, 1,0,  0, 0,1, 1,  1, 1,0,0,  0, 1, 0, 1 };
// 8x8
static std::vector<BINAIRO_TEST_TYPE> binairo_initial  = {-1,-1,-1,-1,-1,-1,-1, 0, -1, 0, 0,-1,-1, 1,-1,-1,  -1, 0,-1,-1,-1, 1,-1, 0,  -1,-1, 1,-1,-1,-1,-1,-1,
                                                           0, 0,-1, 1,-1,-1, 1,-1, -1,-1,-1,-1, 1,-1,-1,-1,   1, 1,-1,-1,-1, 0,-1, 1,  -1, 1,-1,-1,-1,-1,-1, 1};
static std::vector<BINAIRO_TEST_TYPE> binairo_solution = { 0,1,1,0,1,0,1,0,  1,0,0,1,0,1,0,1,  1,0,0,1,0,1,1,0, 0,1,1,0,1,0,0,1,
                                                           0,0,1,1,0,1,1,0,  1,0,0,1,1,0,1,0,  1,1,0,0,1,0,0,1, 0,1,1,0,0,1,0,1};
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

#ifdef TEST_BINAIRO
template <typename T> class BinairoObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        size_t n = (int)pow((double)x.size(), 0.50);
        MatriceUtil<T> mat(n, n);

        bool mismatch = false;
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                mat.set(i, j, x[n*i + j]);

                if ((x[n*i + j] != binairo_initial[n*i + j]) && ((binairo_initial[n*i + j] == 0) || (binairo_initial[n*i + j] == 1)))
                {
                    // Ignore this fixed parameter - replace with right one
                    mat.set(i, j, binairo_initial[n*i + j]);
                    //mismatch = true;
                }
            }
        }

        int cnt_0 = mat.count(0);
        int cnt_1 = mat.count(1);
        int cnt_m1 = mat.count(-1);
        int cnt_other = mat.count(2)+mat.count(-2);

        int cnt_illegal_row = 0;
        int cnt_illegal_col = 0;
        int cnt_illegal_sequ = 0;
        int cnt_illegal_rowcol = 0;
        for (size_t i = 0; i < mat.size_row(); i++)
        {
            if (mat.count_row(i, 0) > (int) (mat.size_col() / 2))
                cnt_illegal_row++;

            if (mat.count_row(i, 1) > (int)(mat.size_col() / 2))
                cnt_illegal_row++;

            if (mat.row_max_sequence(i, 0) > 2)
                cnt_illegal_sequ++;

            if (mat.row_max_sequence(i, 1) > 2)
                cnt_illegal_sequ++;

            for (size_t j = 0; j < mat.size_row(); j++)
            {
                if (i < j)
                {
                    if (mat.row_same(i, j) == true)
                        cnt_illegal_rowcol++;
                }
            }
        }

        for (size_t i = 0; i < mat.size_col(); i++)
        {
            if (mat.count_col(i, 0) > (int) (mat.size_row() / 2))
                cnt_illegal_col++;

            if (mat.count_col(i, 1) > (int)(mat.size_row() / 2))
                cnt_illegal_col++;

            if (mat.col_max_sequence(i, 0) > 2)
                cnt_illegal_sequ++;

            if (mat.col_max_sequence(i, 1) > 2)
                cnt_illegal_sequ++;

            for (size_t j = 0; j < mat.size_col(); j++)
            {
                if (i < j)
                {
                    if (mat.col_same(i, j) == true)
                        cnt_illegal_rowcol++;
                }
            }
        }

        double penality = 0.0;
        if (mismatch == true) 
            penality += 200.0;
        penality += std::fabs((double)cnt_0 - (n * n / 2));
        penality += std::fabs((double)cnt_1 - (n * n / 2));
        penality += std::fabs((double)cnt_m1 - (0));
        penality += std::fabs((double)cnt_other - (0));
        penality += (double)cnt_illegal_row;
        penality += (double)cnt_illegal_col;
        penality += (double)cnt_illegal_sequ;
        penality += (double)cnt_illegal_rowcol;



        double obj = -penality;
        return { obj };
    }
};
#endif

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

        const int POPUL = 5;
        const int N = 5;
        const double MUTRATE = 0.05;
        const int NBIT = 32;

        for (size_t s = 0; s < selectcases.size(); s++)
        {
            for (size_t m = 0; m < mutcases.size(); m++)
            {
                mutinfo._type = mutcases[m];
                for (size_t c = 0; c < crosscases.size(); c++)
                {
                    std::cout << std::endl;
                    std::cout << "SumSameAsPrd function 2x2 = 2+2";
                    galgo::Parameter<_TYPE, NBIT> par1({ (_TYPE)1, (_TYPE)100, 99 });
                    galgo::Parameter<_TYPE, NBIT> par2({ (_TYPE)1, (_TYPE)100, 99 });
                    galgo::GeneticAlgorithm<_TYPE> ga(SumSameAsPrdObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2);
                    ga.mutrate = MUTRATE;
                    ga.Selection = selectcases[s];
                    ga.CrossOver = crosscases[c];
                    ga.run();
                }
            }
        }
    }
}

void TEST_ALL_TYPE()
{
    //---------------------------------------------------
    // Generate all templates to see if compiling/running ok 
    //---------------------------------------------------
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

}

int main()
{
    //---------------------------------------------------
    // Generate all templates to see if compiling/running ok 
    //---------------------------------------------------
    if (false)
    {
        TEST_ALL_TYPE();
    }

    using _TYPE = float;                    // Suppport float, double, char, int, long, ... for parameters
    galgo::MutationInfo<_TYPE> mutinfo;     // Changes mutation info as desired
    mutinfo._sigma          = 1.0;
    mutinfo._sigma_lowest   = 0.01;
    mutinfo._ratio_boundary = 0.10;
    mutinfo._type           = galgo::MutationType::MutationGAM_UncorrelatedNStepSizeBoundary;

    const int       POPUL   = 100;
    const int       N       = 200;  // Number of generation to produce
    const double    MUTRATE = 0.05;
    const int       NBIT    = 63;   // has to remain between 1 and 64
    CROSS<_TYPE>    CROSSType = RealValuedSingleArithmeticRecombination;
    SELECT<_TYPE>   SELECTType  = RWS;

    {
#ifdef TEST_BINAIRO
        {
            //using BINAIRO_TEST_TYPE = char;
            galgo::MutationInfo<BINAIRO_TEST_TYPE> mutinfo;     // Changes mutation info as desired
            mutinfo._sigma = 1.0;
            mutinfo._sigma_lowest = 0.01;
            mutinfo._ratio_boundary = 0.10;
            mutinfo._type = galgo::MutationType::MutationSPM;

            const int       POPUL = 200;
            const int       N = 200000;
            const double    MUTRATE = 0.25;
            const int       NBIT = 2;
            CROSS<BINAIRO_TEST_TYPE>    CROSSType = P1XO;
            SELECT<BINAIRO_TEST_TYPE>   SELECTType = RWS;

            std::cout << std::endl;
            std::cout << "BINAIRO grid NxN";
            int k = 0;
            const int NBinairo = 8;// 4;
            BINAIRO_TEST_TYPE low = -1;
            BINAIRO_TEST_TYPE high = 1;
            std::vector<BINAIRO_TEST_TYPE> vlow(NBinairo * NBinairo);
            std::vector<BINAIRO_TEST_TYPE> vhigh(NBinairo * NBinairo);
            std::vector<BINAIRO_TEST_TYPE> vinit(NBinairo * NBinairo);
            for (size_t i = 0; i < NBinairo * NBinairo; i++)
            {
                vlow[i] = low;
                vhigh[i] = high;
                vinit[i] = binairo_initial[i];
            }

            MatriceUtil<BINAIRO_TEST_TYPE> mat(NBinairo, NBinairo);
            bool mismatch = false;
            std::vector<bool> force_value_flag(NBinairo * NBinairo);
            std::vector<BINAIRO_TEST_TYPE> force_value(NBinairo * NBinairo);
            for (size_t i = 0; i < NBinairo; i++)
            {
                for (size_t j = 0; j < NBinairo; j++)
                {
                    mat.set(i, j, binairo_initial[NBinairo*i + j]);
                    force_value_flag[NBinairo*i + j] = false;
                    force_value[NBinairo*i + j] = -1;
                    if (binairo_initial[NBinairo*i + j] != -1)
                    {
                        force_value_flag[NBinairo*i + j] = true;
                        force_value[NBinairo*i + j] = binairo_initial[NBinairo*i + j];
                    }
                }
            }
            display_binairio<BINAIRO_TEST_TYPE>(mat, false);

            galgo::GeneticAlgorithmN<BINAIRO_TEST_TYPE, NBIT> ga(BinairoObjective<BINAIRO_TEST_TYPE>::Objective, POPUL, N, true, mutinfo, vlow, vhigh, vinit);

            ga.mutrate = MUTRATE;
            ga.Selection = SELECTType;
            ga.CrossOver = CROSSType;
            ga.genstep = 50;
            ga.precision = 2;
            ga.force_value_flag = force_value_flag;
            ga.force_value = force_value;

            //FIXEDPARAM<BINAIRO_TEST_TYPE> f = FixedParameter;
            ga.FixedValue = nullptr; // FixedParameter;
            ga.run();

            const galgo::CHR<BINAIRO_TEST_TYPE> bestParam = ga.result(); //std::shared_ptr<Chromosome<T>>;
            bool ok = true;
            std::vector<BINAIRO_TEST_TYPE> v = bestParam->getParam();
            for (size_t i = 0; i < v.size(); i++)
            {
                if (v[i] != binairo_solution[i])
                {
                    ok = false;
                    break;
                }
            }
            if (ok) std::cout << "Solution found\n";
            else std::cout << "Solution not found\n";
            system("pause");
        }
#endif

        {
            std::cout << std::endl;
            std::cout << "SumSameAsPrd function 2x2 = 2+2";
            galgo::Parameter<_TYPE, NBIT> par1({ (_TYPE)1, (_TYPE)100, 99 }); // an initial value can be added inside the initializer list after the upper bound
            galgo::Parameter<_TYPE, NBIT> par2({ (_TYPE)1, (_TYPE)100, 99 });
            galgo::GeneticAlgorithm<_TYPE> ga(SumSameAsPrdObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2);
            ga.mutrate = MUTRATE;  
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
            ga.Selection = SELECTType;
            ga.CrossOver = CROSSType;
            ga.run();
        }
    }

    system("pause");
}

