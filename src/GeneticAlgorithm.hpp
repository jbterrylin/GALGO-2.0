//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#ifndef GENETICALGORITHM_HPP
#define GENETICALGORITHM_HPP

namespace galgo
{
    /*-------------------------------------------------------------------------------------------------*/
    template <typename T>
    class GeneticAlgorithm
    {
        template <typename K> friend class Population;
        template <typename K> friend class Chromosome;

        template <typename Z>  using FuncKT = std::vector<double>(*)(const std::vector<Z>&);

    protected:
        Population<T> pop;             // population of chromosomes
        std::vector<PAR<T>> param;     // parameter(s) 
        std::vector<T> lowerBound;     // parameter(s) lower bound
        std::vector<T> upperBound;     // parameter(s) upper bound
        std::vector<T> initialSet;     // initial set of parameter(s)
        std::vector<int> idx;          // indexes for chromosome breakdown

    public:
        // objective function pointer
        FuncKT<T> Objective;

        // selection method initialized to roulette wheel selection                                   
        void(*Selection)(Population<T>&) = RWS;

        // cross-over method initialized to 1-point cross-over                                
        void(*CrossOver)(const Population<T>&, CHR<T>&, CHR<T>&) = P1XO;

        // mutation method initialized to single-point mutation 
        void(*Mutation)(CHR<T>&) = SPM;

        // adaptation to constraint(s) method                                        
        void(*Adaptation)(Population<T>&) = nullptr;

        // constraint(s)                               
        std::vector<double>(*Constraint)(const std::vector<T>&) = nullptr;

        MutationInfo<T> mutinfo;

        double covrate = .50;   // cross-over rate
        double mutrate = .05;   // mutation rate   
        double SP = 1.5;        // selective pressure for RSP selection method 
        double tolerance = 0.0; // terminal condition (inactive if equal to zero)
        double recombination_ratio = 0.50; // Real Valued crossover ratio

        int elitpop = 1;   // elit population size
        int matsize;       // mating pool size, set to popsize by default
        int tntsize = 10;  // tournament size
        int genstep = 10;  // generation step for outputting results
        int precision = 10; // precision for outputting results
        bool output;   // control if results must be outputted

        // Prototype to set fixed value of parameters while evolving
        void(*FixedValue)(Population<T>&, int k) = nullptr;
        std::vector<bool> force_value_flag;
        std::vector<T> force_value;

        bool(*StopCondition)(GeneticAlgorithm<T>&) = nullptr;

        // constructor
        template <int...N> GeneticAlgorithm(FuncKT<T> objective, int popsize, int nbgen, bool output, MutationInfo<T> mutinfo, const Parameter<T, N>&...args);
        template <int...N> GeneticAlgorithm(const ConfigInfo<T>& config, const Parameter<T, N>&...args);
        template <int...N> GeneticAlgorithm(const ConfigInfo<T>& config, std::vector<T> init_values, const Parameter<T, N>&...args);

        // run genetic algorithm
        void run();

        // return best chromosome 
        const CHR<T>& result() const;

        // print results for each new generation
        void print(bool force = false) const;

        Population<T>& get_pop() { return pop; }

        int get_genstep() const { return genstep; }

    protected:
        void setMutation(const MutationInfo<T>& mt)
        {
            mutinfo = mt;
            if (mt._type == MutationType::MutationSPM) { Mutation = SPM; }
            else if (mt._type == MutationType::MutationBDM) { Mutation = BDM; }
            else if (mt._type == MutationType::MutationUNM) { Mutation = UNM; }
            else if (mt._type == MutationType::MutationGAM_UncorrelatedOneStepSizeFixed) { Mutation = GAM_UncorrelatedOneStepSizeFixed; }
            else if (mt._type == MutationType::MutationGAM_UncorrelatedOneStepSizeBoundary) { Mutation = GAM_UncorrelatedOneStepSizeBoundary; }
            else if (mt._type == MutationType::MutationGAM_UncorrelatedNStepSize) { Mutation = GAM_UncorrelatedNStepSize; }
            else if (mt._type == MutationType::MutationGAM_UncorrelatedNStepSizeBoundary) { Mutation = GAM_UncorrelatedNStepSizeBoundary; }
            else if (mt._type == MutationType::MutationGAM_sigma_adapting_per_generation) { Mutation = GAM_sigma_adapting_per_generation; }
            else if (mt._type == MutationType::MutationGAM_sigma_adapting_per_mutation) { Mutation = GAM_sigma_adapting_per_mutation; }
            else Mutation = SPM;
        }

        GeneticAlgorithm(const ConfigInfo<T>& config); // No parameters set

        int nbbit;     // total number of bits per chromosome
        int nbgen;     // number of generations
        int nogen = 0; // numero of generation
        int nbparam;   // number of parameters to be estimated
        int popsize;   // population size

        // end of recursion for initializing parameter(s) data
        template <int I = 0, int...N>
        typename std::enable_if<I == sizeof...(N), void>::type init(const TUP<T, N...>&);

        // recursion for initializing parameter(s) data
        template <int I = 0, int...N>
        typename std::enable_if < I < sizeof...(N), void>::type init(const TUP<T, N...>&);

        // check inputs validity
        void check() const;

        void init_from_config(const ConfigInfo<T>& config);

        std::vector<T> _init_values;
    };

    template <typename T> template <int...N>
    GeneticAlgorithm<T>::GeneticAlgorithm(const ConfigInfo<T>& config, std::vector<T> init_values, const Parameter<T, N>&...args)
        : _init_values(init_values)
    {
        init_from_config(config);

        nbbit = sum(N...);
        nbparam = sizeof...(N);
        TUP<T, N...> tp(args...);
        init(tp);
    }

    /*-------------------------------------------------------------------------------------------------*/
    template <typename T, int PARAM_NBIT>
    class GeneticAlgorithmN : public GeneticAlgorithm<T>
    {
    public:
        // constructor
        GeneticAlgorithmN(const ConfigInfo<T>& config,
            std::vector<T>& _lowerBound,
            std::vector<T>& _upperBound,
            std::vector<T>& _initialSet);
    };

    /*-------------------------------------------------------------------------------------------------*/
    template <typename T>
    void GeneticAlgorithm<T>::init_from_config(const ConfigInfo<T>& config)
    {
        setMutation(config.mutinfo); // Mutation is set here

        Objective = config.Objective;
        Selection = config.Selection;
        CrossOver = config.CrossOver;
        Adaptation = config.Adaptation;
        Constraint = config.Constraint;
        FixedValue = config.FixedValue;
        StopCondition = config.StopCondition;

        covrate = config.covrate;
        mutrate = config.mutrate;
        SP = config.SP;
        tolerance = config.tolerance;
        recombination_ratio = config.recombination_ratio;

        elitpop = config.elitpop;
        tntsize = config.tntsize;
        genstep = config.genstep;
        precision = config.precision;

        force_value_flag = config.force_value_flag;
        force_value = config.force_value;

        nbgen = config.nbgen;
        popsize = config.popsize;
        matsize = popsize;          // matsize default to popsize
        output = config.output;

        nogen = 0;
    }

    template <typename T> template <int...N>
    GeneticAlgorithm<T>::GeneticAlgorithm(const ConfigInfo<T>& config, const Parameter<T, N>&...args)
    {
        init_from_config(config);

        nbbit = sum(N...);
        nbparam = sizeof...(N);
        TUP<T, N...> tp(args...);
        init(tp);
    }

    template <typename T>
    GeneticAlgorithm<T>::GeneticAlgorithm(const ConfigInfo<T>& config)
    {
        init_from_config(config);
    }

    // constructor
    template <typename T> template <int...N>
    GeneticAlgorithm<T>::GeneticAlgorithm(FuncKT<T> objective, int popsize, int nbgen, bool output, MutationInfo<T> mutinfo, const Parameter<T, N>&...args)
    {
        setMutation(mutinfo);
        Objective = objective;

        // getting total number of bits per chromosome
        nbbit = sum(N...);
        nbgen = nbgen;

        // getting number of parameters in the pack
        nbparam = sizeof...(N);

        popsize = popsize;
        matsize = popsize;
        output = output;

        // unpacking parameter pack in tuple
        TUP<T, N...> tp(args...);

        // initializing parameter(s) data
        init(tp);
    }

    // constructor
    template <typename T, int PARAM_NBIT>
    GeneticAlgorithmN<T, PARAM_NBIT>::GeneticAlgorithmN(const ConfigInfo<T>& _config, std::vector<T>& _lowerBound, std::vector<T>& _upperBound, std::vector<T>& _initialSet)
        : GeneticAlgorithm<T>(_config)
    {
        for (int i = 0; i < _lowerBound.size(); i++)
        {
            std::vector<T> w;
            w.push_back(_lowerBound[i]); w.push_back(_upperBound[i]); w.push_back(_initialSet[i]);
            Parameter<T, PARAM_NBIT> p(w);
            param.emplace_back(new decltype(p)(p));

            if (i == 0) idx.push_back(0);
            else idx.push_back(idx[i - 1] + PARAM_NBIT);
        }
        lowerBound = _lowerBound;
        upperBound = _upperBound;
        initialSet = _initialSet;

        nbbit = (int)_lowerBound.size()*PARAM_NBIT;
        nbparam = (int)_lowerBound.size();
    };

    /*-------------------------------------------------------------------------------------------------*/

    // end of recursion for initializing parameter(s) data
    template <typename T> template <int I, int...N>
    inline typename std::enable_if<I == sizeof...(N), void>::type
        GeneticAlgorithm<T>::init(const TUP<T, N...>& tp) {}

    // recursion for initializing parameter(s) data
    template <typename T> template <int I, int...N>
    inline typename std::enable_if < I < sizeof...(N), void>::type
        GeneticAlgorithm<T>::init(const TUP<T, N...>& tp)
    {
        // getting Ith parameter in tuple
        auto par = std::get<I>(tp);

        // getting Ith parameter initial data
        const std::vector<T>& data = par.getData();

        // copying parameter data
        param.emplace_back(new decltype(par)(par));

        lowerBound.push_back(data[0]);
        upperBound.push_back(data[1]);

        // if parameter has initial value
        if (data.size() > 2) {
            initialSet.push_back(data[2]);
        }
        // setting indexes for chromosome breakdown
        if (I == 0) {
            idx.push_back(0);
        }
        else {
            idx.push_back(idx[I - 1] + par.size());
        }
        // recursing
        init<I + 1>(tp);
    }

    /*-------------------------------------------------------------------------------------------------*/

    // check inputs validity
    template <typename T>
    void GeneticAlgorithm<T>::check() const
    {
        if (!initialSet.empty()) {
            for (int i = 0; i < nbparam; ++i) {
                if (initialSet[i] < lowerBound[i] || initialSet[i] > upperBound[i]) {
                    throw std::invalid_argument("Error: in class galgo::Parameter<T,N>, initial parameter value cannot be outside the parameter boundaries, please choose a value between its lower and upper bounds.");
                }
            }
            if (initialSet.size() != (unsigned)nbparam) {
                throw std::invalid_argument("Error: in class galgo::GeneticAlgorithm<T>, initial set of parameters does not have the same dimension than the number of parameters, please adjust.");
            }
        }
        if (SP < 1.0 || SP > 2.0) {
            throw std::invalid_argument("Error: in class galgo::GeneticAlgorithm<T>, selective pressure (SP) cannot be outside [1.0,2.0], please choose a real value within this interval.");
        }
        if (elitpop > popsize || elitpop < 0) {
            throw std::invalid_argument("Error: in class galgo::GeneticAlgorithm<T>, elit population (elitpop) cannot outside [0,popsize], please choose an integral value within this interval.");
        }
        if (covrate < 0.0 || covrate > 1.0) {
            throw std::invalid_argument("Error: in class galgo::GeneticAlgorithm<T>, cross-over rate (covrate) cannot outside [0.0,1.0], please choose a real value within this interval.");
        }
        if (genstep <= 0) {
            throw std::invalid_argument("Error: in class galgo::GeneticAlgorithm<T>, generation step (genstep) cannot be <= 0, please choose an integral value > 0.");
        }
    }

    /*-------------------------------------------------------------------------------------------------*/

    // run genetic algorithm
    template <typename T>
    void GeneticAlgorithm<T>::run()
    {
        // checking inputs validity
        check();

        // setting adaptation method to default if needed
        if (Constraint != nullptr && Adaptation == nullptr) {
            Adaptation = DAC;
        }

        // initializing population
        pop = Population<T>(*this);

        if (output) {
            std::cout << "\n Running Genetic Algorithm...\n";
            std::cout << " ----------------------------\n";
        }

        // creating population
        if (_init_values.size() > 0)
            pop.creation(_init_values);
        else
            pop.creation();

        // initializing best result and previous best result
        double bestResult = pop(0)->getTotal();
        double prevBestResult = bestResult;

        // outputting results 
        if (output) print();

        // starting population evolution
        for (nogen = 1; nogen <= nbgen; ++nogen)
        {
            // evolving population
            pop.evolution();

            // getting best current result
            bestResult = pop(0)->getTotal();

            // outputting results
            if (output) print();

            // checking convergence
            if (tolerance != 0.0)
            {
                if (fabs(bestResult - prevBestResult) < fabs(tolerance))
                {
                    break;
                }
                prevBestResult = bestResult;
            }

            if (StopCondition != nullptr)
            {
                if (StopCondition(*this) == true)
                {
                    break;
                }
            }
        }

        // outputting contraint value
        if (Constraint != nullptr)
        {
            // getting best parameter(s) constraint value(s)
            std::vector<double> cst = pop(0)->getConstraint();
            if (output) {
                std::cout << "\n Constraint(s)\n";
                std::cout << " -------------\n";
                for (unsigned i = 0; i < cst.size(); ++i) {
                    std::cout << " C";

                    if (nbparam > 1) {
                        std::cout << std::to_string(i + 1);
                    }
                    std::cout << "(x) = " << std::setw(6) << std::fixed << std::setprecision(precision) << cst[i] << "\n";
                }
                std::cout << "\n";
            }
        }
    }

    /*-------------------------------------------------------------------------------------------------*/

    // return best chromosome
    template <typename T>
    inline const CHR<T>& GeneticAlgorithm<T>::result() const
    {
        return pop(0);
    }

    /*-------------------------------------------------------------------------------------------------*/

    // print results for each new generation
    template <typename T>
    void GeneticAlgorithm<T>::print(bool force) const
    {
        // getting best parameter(s) from best chromosome
        std::vector<T> bestParam = pop(0)->getParam();
        std::vector<double> bestResult = pop(0)->getResult();

        if ((nogen % genstep == 0) || force)
        {
            std::cout << " Generation = " << std::setw(std::to_string(nbgen).size()) << nogen << " |";
            for (int i = 0; i < nbparam; ++i)
            {
                std::cout << " X";
                if (nbparam > 1) {
                    std::cout << std::to_string(i + 1);
                }
                std::cout << " = " << std::setw(2 + precision) << std::fixed << std::setprecision(precision) << bestParam[i] << " |";
            }
            for (unsigned i = 0; i < bestResult.size(); ++i)
            {
                std::cout << " F";
                if (bestResult.size() > 1) {
                    std::cout << std::to_string(i + 1);
                }
                std::cout << "(x) = " << std::setw(12) << std::fixed << std::setprecision(precision) << bestResult[i];
                if (i < bestResult.size() - 1) {
                    std::cout << " |";
                }
                else {
                    std::cout << "\n";
                }
            }

        }
    }


}
#endif