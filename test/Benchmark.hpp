#pragma once

#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

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
double Aaa(T x, T y)
{
    const double pi = 3.14159265358979323846;

    double obj_val =  -20 * std::exp(-0.2*std::sqrt(0.5*(x*x + y*y))) - std::exp(0.5*(std::cos(2 * pi*x) + std::cos(2 * pi*y))) + std::exp(1) + 20;
    return obj_val;
}

template <typename T>
class AaaaObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        double x0 = (double)x[0];
        double x1 = (double)x[1];
        double obj = -Aaa<double>(x0, x1);
        return { obj };
    }
};

#endif