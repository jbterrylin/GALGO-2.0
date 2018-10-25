//=================================================================================================
//                    Copyright (C) 2017 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "Galgo.hpp"

template <typename T>
double Ackley(T x, T y);

// objective class example
template <typename T>
class MyObjective
{
public:
   // objective function example : Rosenbrock function
   // minimizing f(x,y) = (1 - x)^2 + 100 * (y - x^2)^2
   static std::vector<T> Objective(const std::vector<T>& x)
   {
      T obj = -(pow(1-x[0],2)+100*pow(x[1]-x[0]*x[0],2));
      return {obj};
   }
   // NB: GALGO maximize by default so we will maximize -f(x,y)
};

template <typename T>
class AckleyObjective
{
public:
    static std::vector<T> Objective(const std::vector<T>& x)
    {
        T obj = -1 * Ackley(x[0], x[1]);
        return { obj };
    }
};

// constraints example:
// 1) x * y + x - y + 1.5 <= 0
// 2) 10 - x * y <= 0
template <typename T>
std::vector<T> MyConstraint(const std::vector<T>& x)
{
   return 
   {
       //x[0]*x[1]+x[0]-x[1]+1.5,
       //10-x[0]*x[1]
       x[0] - 2,    // x0 <= 2
       x[1] - 2     // x1 <= 2
   };
}

//
// this example is from
// https://en.wikipedia.org/wiki/Test_functions_for_optimization
//
// Ackley's function:
//
// f(x) = -20*exp(-0.2*sqrt(0.5*(x^2 + y^2))) - exp(0.5*(cos(2*pi*x) + cos(2*pi*y))) + exp(1) + 20
// 
// solution is: (0,0)
//
template <typename T>
double Ackley(T x, T y)
{
    const double pi = 3.14159265358979323846;

    T obj_val = -20 * std::exp(-0.2*std::sqrt(0.5*(x*x + y*y))) - std::exp(0.5*(std::cos(2 * pi*x) + std::cos(2 * pi*y))) + std::exp(1) + 20;
    return obj_val;
}

/*
Rastrigin Function
Requires <vector>, <math.h> and C++11
Compilation: g++ fileName.cpp -std=c++11 -o progName
*/
template <typename T>
double pso_rastrigin(std::vector< T > particle) 
{
    T result(10. * static_cast< T > (particle.size())), A(10.), PI(3.14159);
    for (auto dim : particle) {
        result += pow(dim, 2.) - (A * cos(2. * PI * dim));
    }
    return (result);
}

template <typename T>
class rastriginObjective
{
public:
    static std::vector<T> Objective(const std::vector<T>& x)
    {
        T obj = -1 * pso_rastrigin(x);
        return { obj };
    }
};

/*
Griewank Function
Requires <vector>, <math.h> and C++11
*/
double pso_griewank(std::vector< double > particle, void *advanced_settings) {
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
    static std::vector<T> Objective(const std::vector<T>& x)
    {
        T obj = -1 * pso_rastrigin(x);
        return { obj };
    }
};

/*
Styblinski-Tang Function
Requires <vector>, <math> and C++11
Compilation: g++ fileName.cpp -std=c++11 -o progName
Min = (-2.903534,...,--2.903534)
*/
template <typename T>
T pso_styb_tang(std::vector< T > particle)
{
    T result(0.);
    for (auto dim : particle) {
        result += pow(dim, 4.0) - (16. * pow(dim, 2.)) + (5. * dim);
    }
    return (result / 2.);
}

template <typename T>
class StyblinskiTangObjective
{
public:
    static std::vector<T> Objective(const std::vector<T>& x)
    {
        T obj = -1 * pso_styb_tang(x);
        return { obj };
    }
};


template <typename T>
class SumSameAsPrdObjective
{
public:
    static std::vector<T> Objective(const std::vector<T>& x)
    {
        int ix = (int)x[0];
        int iy = (int)x[1];
        T sum = ix + iy;
        T prd = ix * iy;
        T diff = std::fabs(sum - prd);

        T err = 1000000*diff ;
        err +=  1000 * std::fabs(x[0] - (T)ix)   + 1000 * std::fabs(x[1] - (T)iy);

        T obj = -1 * (diff + err);
        return { obj };
    }
};

int main()
{
    // initializing parameters lower and upper bounds
    // an initial value can be added inside the initializer list after the upper bound
    //galgo::Parameter<double> par1({ -2.0,2.0 });
    //galgo::Parameter<double> par2({ -2.0,2.0 });
    // here both parameter will be encoded using 16 bits the default value inside the template declaration
    // this value can be modified but has to remain between 1 and 64
    // initiliazing genetic algorithm
    //galgo::GeneticAlgorithm<double> ga(MyObjective<double>::Objective, 100, 200, true, par1, par2);
    //ga.Constraint = MyConstraint;

    using _TYPE = float;
    galgo::MutationInfo<_TYPE> mutinfo;

    //mutinfo._type = galgo::MutationType::MutationSPM;
    //mutinfo._type = galgo::MutationType::MutationGAM_UncorrelatedOneStepSizeFixed;
    //mutinfo._type = galgo::MutationType::MutationGAM_UncorrelatedOneStepSizeBoundary;
    //mutinfo._type = galgo::MutationType::MutationGAM_UncorrelatedNStepSize;
    //mutinfo._type = galgo::MutationType::MutationGAM_sigma_adapting_per_generation;
    mutinfo._type = galgo::MutationType::MutationGAM_sigma_adapting_per_mutation;

    const int POPUL     = 100;
    const int N         = 200;
    const _TYPE MUTRATE = 0.05;
    const int NBIT      = 64;

    {
        std::cout << std::endl;
        std::cout << "SumSameAsPrd function 2x2 = 2+2";
        galgo::Parameter<_TYPE, NBIT> par1({ 1, 10 });
        galgo::Parameter<_TYPE, NBIT> par2({ 1, 10 });
        galgo::GeneticAlgorithm<_TYPE> ga(SumSameAsPrdObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2);
        ga.mutrate = MUTRATE;
        ga.run();
    }
    
    {
        std::cout << std::endl;
        std::cout << "Rosenbrock function";
        galgo::Parameter<_TYPE, NBIT > par1({ -2.0,2.0 });
        galgo::Parameter<_TYPE, NBIT > par2({ -2.0,2.0 });
        galgo::GeneticAlgorithm<_TYPE> ga(MyObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2);
        ga.mutrate = MUTRATE;
        ga.run();
    }

    {
        std::cout << std::endl;
        std::cout << "Ackley function";
        galgo::Parameter<_TYPE, NBIT > par1({ -4.0,5.0 });
        galgo::Parameter<_TYPE, NBIT > par2({ -4.0,5.0 });
        galgo::GeneticAlgorithm<_TYPE> ga(AckleyObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2);
        ga.mutrate = MUTRATE;
        ga.run();
    }
    
    {
        std::cout << std::endl;
        std::cout << "Rastrigin function";
        galgo::Parameter<_TYPE, NBIT > par1({ -4.0,5.0 });
        galgo::Parameter<_TYPE, NBIT > par2({ -4.0,5.0 });
        galgo::Parameter<_TYPE, NBIT > par3({ -4.0,5.0 });
        galgo::GeneticAlgorithm<_TYPE> ga(rastriginObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2, par3);
        ga.mutrate = MUTRATE;
        ga.run();
    }

    {
        std::cout << std::endl;
        std::cout << "StyblinskiTang function Min = (-2.903534,...,--2.903534)";
        galgo::Parameter<_TYPE, NBIT > par1({ -4.0,4.0 });
        galgo::Parameter<_TYPE, NBIT > par2({ -4.0,4.0 });
        galgo::Parameter<_TYPE, NBIT > par3({ -4.0,4.0 });
        galgo::GeneticAlgorithm<_TYPE> ga(StyblinskiTangObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2, par3);
        ga.mutrate = MUTRATE;
        ga.run();
    }

    {
        std::cout << std::endl;
        std::cout << "Griewank function";
        galgo::Parameter<_TYPE, NBIT > par1({ -4.0,5.0 });
        galgo::Parameter<_TYPE, NBIT > par2({ -4.0,5.0 });
        galgo::Parameter<_TYPE, NBIT > par3({ -4.0,5.0 });
        galgo::GeneticAlgorithm<_TYPE> ga(GriewankObjective<_TYPE>::Objective, POPUL, N, true, mutinfo, par1, par2, par3);
        ga.mutrate = MUTRATE;
        ga.run();
    }

    system("pause");
}
