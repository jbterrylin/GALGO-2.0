//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#ifndef GALGO_H
#define GALGO_H

// #include <WINDOWS.H>      
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

#include <algorithm>
#include <bitset>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <climits>
#include <cmath>
#include <cstring>

#include <fstream>
#include <sstream>
#include <filesystem>

/*-------------------------------------------------------------------------------------------------*/

namespace galgo {

    // forward declarations
    template <typename T>
    struct ConfigInfo;

    template <typename T>
    struct MutationInfo;

    template <typename T>
    class BaseParameter;

    template <typename T, int N = 16>
    class Parameter;

    template <typename T>
    class GeneticAlgorithm;

    template <typename T, int PARAM_NBIT>
    class GeneticAlgorithmN;

    template <typename T>
    class Population;

    template <typename T>
    class Chromosome;

    // convenient typedefs
    template <typename T>
    using CHR = std::shared_ptr<Chromosome<T>>;

    template <typename T>
    using PAR = std::unique_ptr<BaseParameter<T>>;

    template <typename T, int...N>
    using TUP = std::tuple<const Parameter<T, N>&...>;

    using _TYPE = double;       // Suppport float, double, char, int, long, ... for parameters

    const int NBIT = 15;
}

    /*-------------------------------------------------------------------------------------------------*/

// #ifdef _OPENMP 
// #include <omp.h>
// // getting maximum number of threads available
//     static const int MAX_THREADS = omp_get_max_threads();
// #endif

    /*-------------------------------------------------------------------------------------------------*/

#include "Randomize.hpp"
#include "Converter.hpp"
#include "Parameter.hpp"
#include "../test/OtherFuctions.hpp"
#include "Evolution.hpp"
#include "Chromosome.hpp"
#include "Population.hpp"
#include "Config.hpp"
#include "GeneticAlgorithm.hpp"
#include "../test/Crossover.hpp"
#include "../test/Cec17.hpp"
#include "../test/Benchmark.hpp"

//================================================================================================= 

#endif

