#pragma once

#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

template <typename T>
double sphere(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 0; i < particle.size(); i++) {
        sum += pow(particle[i], 2.);
    }
    return sum;
}

template <typename T>
class SphereObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = sphere<double>(xd);
        return { obj };
    }
};

template <typename T>
double axisParallelHyperEllipsoid(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 1; i <= particle.size(); i++) {
        sum += (i * pow(particle[i], 2.));
    }
    return sum;
}

template <typename T>
class AxisParallelHyperEllipsoidObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = axisParallelHyperEllipsoid<double>(xd);
        return { obj };
    }
};

template <typename T>
double rotatedHyperEllipsoid(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 1; i <= particle.size(); i++) {
        double inner (0.);
        for (int j = 1; j <= i; j++)
            inner += pow(particle[j],2);
        sum += inner;
    }

    return sum;
}

template <typename T>
class RotatedHyperEllipsoidObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = rotatedHyperEllipsoid<double>(xd);
        return { obj };
    }
};

template <typename T>
double rotatedHyperEllipsoid(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 1; i <= particle.size(); i++) {
        double inner (0.);
        for (int j = 1; j <= i; j++)
            inner += pow(particle[j],2);
        sum += inner;
    }

    return sum;
}

template <typename T>
class RotatedHyperEllipsoidObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = rotatedHyperEllipsoid<double>(xd);
        return { obj };
    }
};

#endif