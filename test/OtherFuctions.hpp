#ifndef OTHERFUNCTIONS_HPP
#define OTHERFUNCTIONS_HPP

// Review Paper of Various Selection Methods in Genetic Algorithm
// actually this algorithm contain some problem (not sure my guess is correct or not)
//      because add in start point, last pointer will bigger than fitsum
//      so last pointer should take x(0), but in this algo not seem like that
//      but it cannot be if fitnessTillI > fitsum then take x(0) because it maybe x(1), x(2)...
template <typename T>
void StochasticUniformSelection(galgo::Population<T>& x)
{
    // adjusting all fitness to positive values
    x.adjustFitness();

    // computing fitness sum
    double fitsum = x.getSumFitness();
    auto matsize = x.matsize();

    auto distance = fitsum / matsize;

    double start = galgo::uniform<double>(0.0, distance);

    std::vector<double> pointers = {};

    for(int i=0; i<x.matsize(); i++) {
        pointers.push_back(start + i * distance);
    }


    int i=0;
    double fitnessTillI = x(i)->fitness;
    for(int p=0; p<pointers.size(); p++ ) {
        while(fitnessTillI < pointers[p]) {
            i++;
            fitnessTillI += x( i )->fitness;
        }
        x.select( i );
    }
}

template <typename T>
void NoSelection(galgo::Population<T>& x)
{
    for (int i = 0; i < x.matsize(); i++)
    {
        x.select(i);
    }
}

template <typename T>
void NoMutation(galgo::CHR<T>& chr)
{
    return;
}

#endif
