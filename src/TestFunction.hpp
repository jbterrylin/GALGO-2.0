//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#ifndef TESTFUNCTION_HPP
#define TESTFUNCTION_HPP

//--------------------------
// constraints example:
//--------------------------
template <typename T>
std::vector<double> MyConstraint(const std::vector<T>& x)
{
    double x0 = (double)x[0];
    double x1 = (double)x[1];
    return 
    {
       //x[0]*x[1]+x[0]-x[1]+1.5,   // 1) x * y + x - y + 1.5 <= 0
       //10-x[0]*x[1]               // 2) 10 - x * y <= 0
       x0 - 2,    // x0 <= 2
       x1 - 2     // x1 <= 2
    };
}

//--------------------------
// functions
//--------------------------
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
    double result(10. * static_cast<T> (particle.size())), A(10.);
    for (auto dim : particle) {
        result += pow(dim, 2.) - (A * cos(2. * PI * dim));
    }
    return (result);
}

/*
Griewank Function
*/
template <typename T>
double pso_griewank(std::vector< T > particle) 
{
    double sum(0.), product(1.);
    for (int i = 0; i < particle.size(); i++) {
        sum += pow(particle[i], 2.);
        product *= cos(particle[i] / sqrt(i + 1));
    }
    return (1. + (sum / 4000.) - product);
}

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

//--------------------------
// Objectives
//--------------------------
template <typename T>
class rastriginObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = -pso_rastrigin<double>(xd);
        return { obj };
    }
};

template <typename T>
class GriewankObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = -pso_griewank<double>(xd);
        return { obj };
    }
};

template <typename T>
class StyblinskiTangObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for(size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = -pso_styb_tang<double>(xd);
        return { obj };
    }
};

template <typename T>
class RosenbrockObjective
{
public:
    // objective function example : Rosenbrock function
    // minimizing f(x,y) = (1 - x)^2 + 100 * (y - x^2)^2
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        double x0 = (double)x[0];
        double x1 = (double)x[1];
        double obj = -(pow(1.0 - x0, 2.0) + 100 * pow(x1 - x0*x0, 2.0));
        return { obj };
    }
    // NB: GALGO maximize by default so we will maximize -f(x,y)
};

template <typename T>
class AckleyObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        double x0 = (double)x[0];
        double x1 = (double)x[1];
        double obj = -Ackley<double>(x0, x1);
        return { obj };
    }
};

template <typename T>
class SumSameAsPrdObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        double x0 = (double)x[0];
        double x1 = (double)x[1];

        int ix = (int)x0;
        int iy = (int)x1;
        double sum = ix + iy;
        double prd = ix * iy;
        double diff  = std::fabs(sum - prd);

        double err = 1000 * diff * diff;;
        err += (100 * std::fabs(x0 - ix)* std::fabs(x0 - ix) + 100 * std::fabs(x1 - iy)* std::fabs(x1 - iy));

        double obj = -(diff + err);
        return { obj };
    }
};


#endif