#ifndef CROSSOVER_HPP
#define CROSSOVER_HPP

// These alr in Evolution.hpp
    // Single point crossover 
    // Two point crossover
    // Arithmetic crossover

// Ring Crossover
template <typename T>
void RingCrossover(const galgo::Population<T>& x, std::vector< galgo::CHR<T> >& chr)
{
    // choosing randomly 2 chromosomes from mating population
    int idx1 = galgo::uniform<int>(0, x.matsize());
    int idx2 = galgo::uniform<int>(0, x.matsize());
    if (x.matsize() >= 2)
    {
        while (idx1 == idx2) { idx2 = galgo::uniform<int>(0, x.matsize()); } // find not unique parents
    }

    // std::cout << "(*x[idx1]).chr" << (*x[idx1]).chr << std::endl;
    // std::cout << "(*x[idx2]).chr" << (*x[idx2]).chr << std::endl;

    // choosing randomly a position for cross-over
    int pos = galgo::uniform<int>(0, (*x[idx1]).chr.size());
    // std::cout << "pos" << pos << std::endl;
    auto reverseBits = [](std::string bits) 
    { 
        int n = bits.length();

        // Swap character starting from two
        // corners
        for (int i = 0; i < n / 2; i++)
            std::swap(bits[i], bits[n - i - 1]);
        return bits;
    };

    chr[0]->setPortion(*x[idx1], 0, pos);
    chr[0]->chr = reverseBits(chr[0]->chr) + (*x[idx2]).chr.substr(0, (*x[idx1]).size() - pos - 1);

    chr[1]->chr = (*x[idx1]).chr.substr(pos + 1) + reverseBits((*x[idx2]).chr.substr((*x[idx1]).size() - pos - 1));

    // std::cout << "chr[0]->chr" << chr[0]->chr << std::endl;
    // std::cout << "chr[1]->chr" << chr[1]->chr << std::endl;

    // double r = chr[0]->recombination_ratio();
    // const galgo::Chromosome<T>& chrmat1 = *x[idx1];
    // const galgo::Chromosome<T>& chrmat2 = *x[idx2];

    // Transmit sigma
    transmit_sigma<T>(chr[0]->recombination_ratio(), *x[idx1], *x[idx2], chr[0], chr[1]);
}

template <typename T>
void HeuristicCrossover(const galgo::Population<T>& x, std::vector< galgo::CHR<T> >& chr)
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
        chr[0]->setPortion(*x[idx1], 0);
        // change idx2
        for (int i = 0; i < chr[1]->nbgene(); i++)
        {
            chr[1]->initGene(i, (T)( chrmat1.get_value(i) + (r * (chrmat1.get_value(i) - chrmat2.get_value(i))) ));
        }
    } else {
        chr[1]->setPortion(*x[idx2], 0);
        // change idx1
        for (int i = 0; i < chr[0]->nbgene(); i++)
        {
            chr[0]->initGene(i, (T)( chrmat2.get_value(i) + (r * (chrmat2.get_value(i) - chrmat1.get_value(i))) ));
        }
    }

    // Transmit sigma
    transmit_sigma<T>(chr[0]->recombination_ratio(), chrmat1, chrmat2, chr[0], chr[1]);
}

template <typename T>
void IntermediateCrossover(const galgo::Population<T>& x, std::vector< galgo::CHR<T> >& chr)
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
    
    for (int i = 0; i < chr[1]->nbgene(); i++)
    {
        chr[0]->initGene(i, (T)( chrmat1.get_value(i) + galgo::uniform<float>(0, 1) * r * (chrmat2.get_value(i) - chrmat1.get_value(i)) ));
    }
    for (int i = 0; i < chr[1]->nbgene(); i++)
    {
        chr[1]->initGene(i, (T)( chrmat2.get_value(i) + galgo::uniform<float>(0, 1) * r * (chrmat1.get_value(i) - chrmat2.get_value(i)) ));
    }

    // Transmit sigma
    transmit_sigma<T>(chr[0]->recombination_ratio(), chrmat1, chrmat2, chr[0], chr[1]);
}

// turn to multi
template <typename T>
void CollectiveCrossover(const galgo::Population<T>& x, std::vector< galgo::CHR<T> >& chr)
{
    // Algorithm 1(Crossover)
    std::vector< galgo::Chromosome<T> >  popc;

    for(int i=0; i < 2 * ( 0.7 *( chr.size()/2 ) ); i++) {
    // for(int i=0; i< chr.size(); i++) {
        popc.push_back( galgo::Chromosome<T>( *chr[i] ) );
    }

    int matingPoolSize = 0;
    for (int n = 1; n <= x.popsize(); n++) matingPoolSize += n;

    for(int popi=0; popi<popc.size(); popi++) {
        for (int i = 0; i < chr[0]->nbgene(); i++) {
            int targetChromosomeInPool = galgo::uniform<int>(0, matingPoolSize);
            int temp = 0;
            for (int n = 0; n < x.popsize(); n++) {
                temp += n;
                if(targetChromosomeInPool < temp){
                    popc[popi].initGene(i, (T)( x(n-1)->get_value(i) ));
                    break;
                } else if (x.popsize() == n+1) {
                    popc[popi].initGene(i, (T)( x(n)->get_value(i) ));
                    break;
                }
            }
        }
    }

    // mutate
    std::vector< galgo::Chromosome<T> >  popm;
    for(int i=0; i < 0.3 *chr.size(); i++) {
        int target = galgo::uniform<int>(0, x.popsize());
        popm.push_back( galgo::Chromosome<T>( *x(target) ) );
    }

    if (chr[0]->mutrate() != 0.0) {
        for(int popm_i=0; popm_i < popm.size(); popm_i++) {
            double mutrate = popm[popm_i].mutrate();
            const std::vector<T>& lowerBound = popm[popm_i].lowerBound();
            const std::vector<T>& upperBound = popm[popm_i].upperBound();

            std::normal_distribution<double> distribution01(0.0, 1.0);

            // int i  = galgo::uniform<int>(0, popm[popm_i].nbgene());
            // looping on number of genes
            for (int i = 0; i < popm[popm_i].nbgene(); ++i)
            {
                // generating a random probability
                if (galgo::proba(galgo::rng) <= mutrate)
                {
                    T value = popm[popm_i].get_value(i);
                    double sigma = popm[popm_i].get_sigma(i);

                    if (sigma < 0.00000000001) // never copied from parent
                    {
                        sigma = ((double)(upperBound[i] - lowerBound[i])) * popm[popm_i].mutinfo()._ratio_boundary;
                        if (sigma < popm[popm_i].mutinfo()._sigma_lowest)
                            sigma = popm[popm_i].mutinfo()._sigma_lowest;
                        popm[popm_i].sigma_update(i, sigma);
                    }

                    std::normal_distribution<double> distribution((double)value, sigma);
                    double norm = distribution(galgo::rng);
                    T new_value = (T)(std::min(std::max((T)norm, lowerBound[i]), upperBound[i])) + value;
                    popm[popm_i].initGene(i, new_value);
                }
            }
        }
    }

    std::vector< galgo::Chromosome<T> > pop = popc;
    pop.insert(pop.end(), popm.begin(), popm.end());
    for(int i=0; i<x.popsize(); i++) {
        pop.push_back(*x(i));
    }

    // get highest fitness to list
    for(int i=0; i<pop.size(); i++) {
        pop[i].evaluate();
    }
    std::sort(pop.begin(),pop.end(),[](const galgo::Chromosome<T>& chr1,const galgo::Chromosome<T>& chr2)->bool{return chr1.fitness > chr2.fitness;});
    for(int i=0; i<chr.size(); i++) {
        chr[i]->chr = pop[i].chr;
    }
}

template <typename T>
void HighDimensionalGeneticAlgorithmToolboxCrossover(const galgo::Population<T>& x, std::vector< galgo::CHR<T> >& chr)
{
    // std::vector< galgo::Chromosome<T> > matPool;
    // for(int i=0; i<x.matsize(); i++) {
    //     matPool.push_back(*x[i]);
    // }

    double p1 = 0.9, p2=0.9;
    std::vector<int> selectedGenes;
    for(int i=0; i<(*x[0]).nbgene(); i++) {
        if(galgo::proba(galgo::rng) <= p1) {
            selectedGenes.push_back(i);
        } else {
            selectedGenes.push_back(-1);
        }
    }

    for(int chri=0; chri<chr.size(); chri=chri+2) {
        for(int genei=0; genei<(*x[0]).nbgene(); genei++) {
            // if -1 means no need to change
            if(selectedGenes[genei] == -1 || (selectedGenes[genei] != -1 && galgo::proba(galgo::rng) > p2)) {
                // std::cout << "chr[chri]->chr" << std::endl;
                // std::cout << x[chri]->chr << std::endl;
                // std::cout << chr[chri]->chr << std::endl;
                chr[chri]->initGene(genei, x[chri]->get_value(genei));
                // std::cout << chr[chri]->chr << std::endl;
                // std::cout << "------------" << std::endl;
                chr[chri+1]->initGene(genei, x[chri+1]->get_value(genei));
            } else {
                int nGeneBit = x[0]->size() / x[0]->nbgene();
                int pos = galgo::uniform<int>(0, nGeneBit);

                // std::cout << "chr[chri]->chr" << std::endl;
                // std::cout << "nGeneBit" << nGeneBit << "pos" << pos << std::endl;
                // std::cout << x[chri]->chr << std::endl;
                // std::cout << chr[chri]->chr << std::endl;

                std::string geneBits1 = x[chri]->chr.substr(genei* x[0]->size()/x[0]->nbgene(), nGeneBit);
                std::string geneBits2 = x[chri+1]->chr.substr(genei* x[0]->size()/x[0]->nbgene(), nGeneBit);

                // std::cout << "geneBits1" << geneBits1 << std::endl;
                // std::cout << "geneBits2" << geneBits2 << std::endl;

                auto genes1temp = geneBits1.substr(pos, nGeneBit);
                geneBits1 = geneBits1.substr(0, pos) + geneBits2.substr(pos, nGeneBit);
                geneBits2 = geneBits2.substr(0, pos) + genes1temp;

                // std::cout << "geneBits1" << geneBits1 << std::endl;
                // std::cout << "geneBits2" << geneBits2 << std::endl;

                chr[chri]->chr += geneBits1;
                chr[chri+1]->chr += geneBits2;

                // std::cout << "chr[chri]->chr" << chr[chri]->chr << std::endl;
                // std::cout << "chr[chri+1]->chr" << chr[chri+1]->chr << std::endl;
                // std::cout << "------------" << std::endl;
            }
        }
        transmit_sigma<T>(chr[0]->recombination_ratio(), *x[chri], *x[chri+1], chr[chri], chr[chri+1]);
    }
}

template <typename T>
void FrontRearCrossover(const galgo::Population<T>& x, std::vector< galgo::CHR<T> >& chr)
{
    // choosing randomly 2 chromosomes from mating population
    int idx1 = galgo::uniform<int>(0, x.matsize());
    int idx2 = galgo::uniform<int>(0, x.matsize());
    if (x.matsize() >= 2)
    {
        while (idx1 == idx2) { idx2 = galgo::uniform<int>(0, x.matsize()); } // find not unique parents
    }

    // choosing randomly a position for cross-over
    int pos = galgo::uniform<int>(0, chr[0]->size());

    chr[0]->chr = (*x[idx2]).chr.substr(chr[0]->size() - pos);
    chr[0]->setPortion(*x[idx1], pos);
    chr[1]->setPortion(*x[idx2], 0, chr[0]->size() - pos - 1);
    chr[1]->chr += (*x[idx1]).chr.substr(0, pos);

    double r = chr[0]->recombination_ratio();
    const galgo::Chromosome<T>& chrmat1 = *x[idx1];
    const galgo::Chromosome<T>& chrmat2 = *x[idx2];

    // Transmit sigma
    transmit_sigma<T>(r, chrmat1, chrmat2, chr[0], chr[1]);
}

template <typename T>
std::vector< galgo::Chromosome<T> > hybridCrossoverAlgo1(const galgo::Population<T>& x, std::vector< galgo::Chromosome<T> > chr) {
    int child_size = 0;
    int pop_pat [chr[0].size()] = {  };

    for (int idx = 0; idx < chr.size(); idx++) {
        for (int i = 0; i < (*x[idx]).chr.length(); i++) {
            if((*x[idx]).chr[i] == '1') {
                pop_pat[i]++;
            }
        }
    }

    do {
        // choosing randomly 2 chromosomes from mating population
        int idx1 = galgo::uniform<int>(0, x.matsize());
        int idx2 = galgo::uniform<int>(0, x.matsize());
        if (x.matsize() >= 2)
        {
            while (idx1 == idx2) { idx2 = galgo::uniform<int>(0, x.matsize()); } // find not unique parents
        }

        for (int j = 0; j <  (*x[0]).chr.length(); j++) {
            if(pop_pat[j]-chr.size()/2 < 3) {
                chr[child_size].chr += (*x[idx2]).chr[j];
                chr[child_size+1].chr += (*x[idx1]).chr[j];
            } else {
                auto temp = galgo::uniform<double>(0, 1);
                if(temp < 1 - (double)pop_pat[j]/chr.size()) {
                    chr[child_size].chr += '1';
                    chr[child_size+1].chr += '1';
                } else {
                    chr[child_size].chr += '0';
                    chr[child_size+1].chr += '0';
                }
            }

            pop_pat[j] = 0;
            for (int idx = 0; idx < chr.size(); idx++) {
                if((*x[idx]).chr[j] == '1') {
                    pop_pat[j]++;
                }
            }
        }
        child_size += 2;
    } while(child_size < chr.size());
    return chr;
}

template <typename T>
std::vector< galgo::Chromosome<T> > hybridCrossoverAlgo2(const galgo::Population<T>& x, std::vector< galgo::Chromosome<T> > chr, std::vector< galgo::Chromosome<T> > algo1result) {
    for(int i=0; i< chr.size(); i = i + 2) {
        int idx1 = galgo::uniform<int>(0, algo1result.size());
        int idx2 = galgo::uniform<int>(0, algo1result.size());
        if (algo1result.size() >= 2)
        {
            while (idx1 == idx2) { idx2 = galgo::uniform<int>(0, algo1result.size()); } // find not unique parents
        }

        int size_one = galgo::uniform<int>(1, algo1result[0].chr.size()/2 + 1);
        int a_start, b_start;
        if(size_one == algo1result[0].chr.size()/2) { 
            a_start = algo1result[0].size()/2;
            b_start = 0;
        }else {
            a_start = algo1result[0].size()/2 + galgo::uniform<int>(0, algo1result[0].chr.size()/2 - size_one + 1);
            b_start = galgo::uniform<int>(0, algo1result[0].chr.size()/2 - size_one + 1) ;
        }
        chr[i].chr = algo1result[idx1].chr;
        chr[i+1].chr = algo1result[idx2].chr;

        std::string tempa = chr[i].chr.substr (a_start, size_one);
        std::string tempb = chr[i+1].chr.substr (b_start, size_one);
        chr[i].chr.replace(a_start, size_one, tempb);
        chr[i+1].chr.replace(b_start, size_one , tempa);

        int size_two = galgo::uniform<int>(1, algo1result[0].chr.size()/2);
        int c_start, d_start;
        c_start = galgo::uniform<int>(0, algo1result[0].chr.size()/2 - size_two);
        d_start = galgo::uniform<int>(0, algo1result[0].chr.size()/2 - size_two);

        std::string tempc = chr[i+1].chr.substr (c_start, size_two);
        std::string tempd = chr[i].chr.substr (d_start, size_two);
        chr[i].chr.replace(c_start, size_two, tempc);
        chr[i+1].chr.replace(d_start, size_two , tempd);
    }

    return chr;
}

// template <typename T>
// std::vector< galgo::Chromosome<T> > hybridCrossoverAlgo3(const galgo::Population<T>& x, std::vector< galgo::Chromosome<T> > chr, std::vector< galgo::Chromosome<T> > algo2result) {
//     for(int i=0; i< chr.size(); i++) {
//         chr[i].chr = algo2result[i].chr;
//         std::cout << "chr[i].chr" << chr[i].chr << std::endl;
//         int start = galgo::uniform<int>(0, algo2result[0].chr.size());
//         std::cout << "start" << start << std::endl;

//         int j = start;
//         int counter = -1;
//         do {
//             std::cout << "---------------------" << std::endl;
//             if(j == algo2result[0].chr.size()/2 || j == 0)
//                 counter = 1;
//             else
//                 if(j != algo2result[0].chr.size()) {
//                     std::cout << "chr[i].chr.substr (j, 1)" << chr[i].chr << std::endl;
//                     chr[i].evaluate();
//                     auto oldFitness = chr[i].fitness;
//                     if(chr[i].chr.substr (j, 1) == "0")
//                         chr[i].chr.replace(j, 1, "1");
//                     else
//                         chr[i].chr.replace(j, 1, "0");

//                     chr[i].evaluate();
//                     if(chr[i].fitness < oldFitness)
//                         if(chr[i].fitness < oldFitness) 
//                             if(chr[i].chr.substr (j, 1) == "0")
//                                 chr[i].chr.replace(j, 1, "1");
//                             else
//                                 chr[i].chr.replace(j, 1, "0");
//                     std::cout << "chr[i].chr.substr (j, 1)" << chr[i].chr << std::endl;
//                     std::cout << "oldFitness" << oldFitness << std::endl;
//                     std::cout << "chr[i].fitness" << chr[i].fitness << std::endl;
//                 }
//             std::cout << "---------------------" << std::endl;
//             j += counter;
//         } while(algo2result[0].chr.size() != j);
//         std::cout << "chr[i].chr" << chr[i].chr << std::endl;
//     }
//     return chr;
// }

template <typename T>
std::vector< galgo::Chromosome<T> > hybridCrossoverAlgo3(const galgo::Population<T>& x, std::vector< galgo::Chromosome<T> > chr, std::vector< galgo::Chromosome<T> > algo2result) {
    for(int i=0; i< chr.size(); i++) {
        chr[i].chr = algo2result[i].chr;
        algo2result[i].evaluate();
        // std::cout << "chr[i].chr" << chr[i].chr << std::endl;
        int start = galgo::uniform<int>(0, algo2result[0].chr.size());

        int j = start;
        int counter = -1;
        do {
            // std::cout << "---------------------" << std::endl;
            if(j == algo2result[0].chr.size()/2 || j == 0)
                counter = 1;
            else
                if(j != algo2result[0].chr.size()) {
                    // std::cout << "chr[i].chr.substr (j, 1)" << chr[i].chr << std::endl;
                    if(chr[i].chr.substr (j, 1) == "0")
                        chr[i].chr.replace(j, 1, "1");
                    else
                        chr[i].chr.replace(j, 1, "0");
                    // std::cout << "chr[i].chr.substr (j, 1)" << chr[i].chr << std::endl;
                }
            // std::cout << "---------------------" << std::endl;
            j += counter;
        } while(algo2result[0].chr.size() != j);
        chr[i].evaluate();
        
        // if result is bad then change back
        if(chr[i].fitness < algo2result[i].fitness) {
            chr[i].chr = algo2result[i].chr;
        }
        // std::cout << "chr[i].chr" << chr[i].chr << std::endl;
    }
    return chr;
}

// radius = 3
// todo
template <typename T>
void HybridCrossover(const galgo::Population<T>& x, std::vector< galgo::CHR<T> >& chr)
{
    std::vector< galgo::Chromosome<T> > temp;
    for(int i=0; i< chr.size(); i++) {
        temp.push_back( galgo::Chromosome<T>( *chr[i] ) );
    }
    auto algoresult = hybridCrossoverAlgo1(x, temp);
    algoresult = hybridCrossoverAlgo2(x, temp, algoresult);
    algoresult = hybridCrossoverAlgo3(x, temp, algoresult);

    for(int i=0; i< algoresult.size(); i++) {
        chr[i]->chr = algoresult[i].chr;
    }

    double r = chr[0]->recombination_ratio();
    for(int i=0; i< chr.size(); i = i + 2) {
        const galgo::Chromosome<T>& chrmat1 = *x[i];
        const galgo::Chromosome<T>& chrmat2 = *x[i+1];

        // Transmit sigma
        transmit_sigma<T>(r, chrmat1, chrmat2, chr[i], chr[i+1]);
    }

    // mutation
}

#endif