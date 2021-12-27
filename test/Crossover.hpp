#ifndef CROSSOVER_HPP
#define CROSSOVER_HPP

template <typename T>
using MULTI_CROSS = void (*)(const galgo::Population<T>&);

template <typename _TYPE>
std::vector<MULTI_CROSS<_TYPE>> crosscases = {
    RealValuedWholeArithmeticRecombination<_TYPE>
};

// These alr in Evolution.hpp
    // Single point crossover 
    // Two point crossover
    // Arithmetic crossover

// Ring Crossover
template <typename T>
void RingCrossover(const galgo::Population<T>& x, galgo::CHR<T>& chr1, galgo::CHR<T>& chr2)
{
    // choosing randomly 2 chromosomes from mating population
    int idx1 = galgo::uniform<int>(0, x.matsize());
    int idx2 = galgo::uniform<int>(0, x.matsize());
    if (x.matsize() >= 2)
    {
        while (idx1 == idx2) { idx2 = galgo::uniform<int>(0, x.matsize()); } // find not unique parents
    }

    // choosing randomly a position for cross-over
    int pos = galgo::uniform<int>(0, chr1->size());

    auto reverseBits = [](std::string bits) 
    { 
        int n = bits.length();

        // Swap character starting from two
        // corners
        for (int i = 0; i < n / 2; i++)
            std::swap(bits[i], bits[n - i - 1]);
        return bits;
    };

    chr1->setPortion(*x[idx1], 0, pos);
    chr1->chr = reverseBits(chr1->chr) + (*x[idx2]).chr.substr(0, (*x[idx1]).size() - pos - 1);

    chr2->chr = (*x[idx1]).chr.substr(pos + 1) + reverseBits((*x[idx2]).chr.substr((*x[idx1]).size() - pos - 1));

    double r = chr1->recombination_ratio();
    const galgo::Chromosome<T>& chrmat1 = *x[idx1];
    const galgo::Chromosome<T>& chrmat2 = *x[idx2];

    // Transmit sigma
    transmit_sigma<T>(r, chrmat1, chrmat2, chr1, chr2);
}

template <typename T>
void HeuristicCrossover(const galgo::Population<T>& x, galgo::CHR<T>& chr1, galgo::CHR<T>& chr2)
{
    // choosing randomly 2 chromosomes from mating population
    int idx1 = galgo::uniform<int>(0, x.matsize());
    int idx2 = galgo::uniform<int>(0, x.matsize());
    if (x.matsize() >= 2)
    {
        while (idx1 == idx2) { idx2 = galgo::uniform<int>(0, x.matsize()); } // find not unique parents
    }

    double r = 1.2;
    const galgo::Chromosome<T>& chrmat1 = *x[idx1];
    const galgo::Chromosome<T>& chrmat2 = *x[idx2];

    // std::cout << "(*x[idx1]).getTotal()[0]" << (*x[idx1]).getTotal() << std::endl;
    // std::cout << "(*x[idx2]).getTotal()[0]" << (*x[idx2]).getTotal() << std::endl;

    if(idx1 < idx2) {
        chr1->setPortion(*x[idx1], 0);
        // change idx2
        for (int i = 0; i < chr2->nbgene(); i++)
        {
            chr2->initGene(i, (T)( chrmat1.get_value(i) + (r * (chrmat1.get_value(i) - chrmat2.get_value(i))) ));
        }
    } else {
        chr2->setPortion(*x[idx2], 0);
        // change idx1
        for (int i = 0; i < chr1->nbgene(); i++)
        {
            chr1->initGene(i, (T)( chrmat2.get_value(i) + (r * (chrmat2.get_value(i) - chrmat1.get_value(i))) ));
        }
    }

    // Transmit sigma
    transmit_sigma<T>(chr1->recombination_ratio(), chrmat1, chrmat2, chr1, chr2);
}

template <typename T>
void IntermediateCrossover(const galgo::Population<T>& x, galgo::CHR<T>& chr1, galgo::CHR<T>& chr2)
{
    // choosing randomly 2 chromosomes from mating population
    int idx1 = galgo::uniform<int>(0, x.matsize());
    int idx2 = galgo::uniform<int>(0, x.matsize());
    if (x.matsize() >= 2)
    {
        while (idx1 == idx2) { idx2 = galgo::uniform<int>(0, x.matsize()); } // find not unique parents
    }

    double r = 1;
    const galgo::Chromosome<T>& chrmat1 = *x[idx1];
    const galgo::Chromosome<T>& chrmat2 = *x[idx2];
    
    for (int i = 0; i < chr2->nbgene(); i++)
    {
        chr1->initGene(i, (T)( chrmat1.get_value(i) + galgo::uniform<float>(0, 1) * r * (chrmat2.get_value(i) - chrmat1.get_value(i)) ));
    }
    for (int i = 0; i < chr2->nbgene(); i++)
    {
        chr2->initGene(i, (T)( chrmat2.get_value(i) + galgo::uniform<float>(0, 1) * r * (chrmat1.get_value(i) - chrmat2.get_value(i)) ));
    }

    // Transmit sigma
    transmit_sigma<T>(chr1->recombination_ratio(), chrmat1, chrmat2, chr1, chr2);
}

template <typename T>
void CollectiveCrossover(const galgo::Population<T>& x, galgo::CHR<T>& chr1, galgo::CHR<T>& chr2)
{
    int matingPoolSize = 0;
    for (int n = 1; n <= x.popsize(); n++) matingPoolSize += n;

    for (int i = 0; i < chr1->nbgene(); i++) {
        int targetChromosomeInPool = galgo::uniform<int>(0, matingPoolSize);
        int temp = 0;
        for (int n = 0; n < x.popsize(); n++) {
            temp += n;
            if(targetChromosomeInPool < temp){
                chr1->initGene(i, (T)( x(n-1)->get_value(i) ));
                break;
            } else if (x.popsize() == n+1) {
                chr1->initGene(i, (T)( x(n)->get_value(i) ));
            }
        }
    }

    for (int i = 0; i < chr2->nbgene(); i++) {
        int targetChromosomeInPool = galgo::uniform<int>(0, matingPoolSize);
        int temp = 0;
        for (int n = 0; n < x.popsize(); n++) {
            temp += n;
            if(targetChromosomeInPool < temp){
                chr2->initGene(i, (T)( x(n-1)->get_value(i) ));
                break;
            } else if (x.popsize() == n+1) {
                chr2->initGene(i, (T)( x(n)->get_value(i) ));
            }
        }
    }
    
    int idx1 = galgo::uniform<int>(0, x.matsize());
    int idx2 = galgo::uniform<int>(0, x.matsize());
    if (x.matsize() >= 2)
    {
        while (idx1 == idx2) { idx2 = galgo::uniform<int>(0, x.matsize()); } // find not unique parents
    }

    const galgo::Chromosome<T>& chrmat1 = *x[idx1];
    const galgo::Chromosome<T>& chrmat2 = *x[idx2];

    // Transmit sigmas
    transmit_sigma<T>(chr1->recombination_ratio(), chrmat1, chrmat2, chr1, chr2);
}

template <typename T>
void HighDimensionalGeneticAlgorithmToolboxCrossover(const galgo::Population<T>& x)
{
    int idx1 = galgo::uniform<int>(0, x.matsize());
    int idx2 = galgo::uniform<int>(0, x.matsize());
    if (x.matsize() >= 2)
    {
        while (idx1 == idx2) { idx2 = galgo::uniform<int>(0, x.matsize()); } // find not unique parents
    }

    // const galgo::Chromosome<T>& chrmat1 = *x[idx1];
    // const galgo::Chromosome<T>& chrmat2 = *x[idx2];

    // Transmit sigmas
    // transmit_sigma<T>(chr1->recombination_ratio(), chrmat1, chrmat2, chr1, chr2);
}

#endif