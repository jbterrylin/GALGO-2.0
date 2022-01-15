//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================
#pragma once

#ifndef CONFIG_HPP
#define CONFIG_HPP

namespace galgo 
{
    enum class MutationType 
    {
        MutationSPM,
        MutationBDM,
        MutationUNM,
        MutationGAM_UncorrelatedOneStepSizeFixed,
        MutationGAM_UncorrelatedOneStepSizeBoundary,
        MutationGAM_UncorrelatedNStepSize,
        MutationGAM_UncorrelatedNStepSizeBoundary,
        MutationGAM_sigma_adapting_per_generation,
        MutationGAM_sigma_adapting_per_mutation,
        MutationNoMutation,
    };

    template <typename T>
    struct MutationInfo
    {
        MutationInfo<T>() :
            _type(MutationType::MutationSPM),
            _sigma(1.0),
            _ratio_boundary(1.0 / 6.0),
            _sigma_lowest(0.0001)
        {
        }

        MutationType _type;
        double _sigma;
        double _ratio_boundary;
        double _sigma_lowest;
    };

    template <typename ParamTYPE>
    struct ConfigInfo
    {
        ConfigInfo() : mutinfo()
        {
            // DEFAULTS
            covrate = .50;
            mutrate = .05;
            SP = 1.5;
            tolerance = 0.0;
            recombination_ratio = 0.50;

            elitpop = 1;
            tntsize = 10;
            genstep = 10;
            precision = 10;

            Objective = nullptr;
            Selection = RWS;
            CrossOver = P1XO;
            isMultiCrossover = false;
            //Mutation = SPM; // derived from by mutinfo._type
            Adaptation = nullptr;
            Constraint = nullptr;
            FixedValue = nullptr;
            StopCondition = nullptr;

            nbgen = 10;
            popsize = 10;
            output = false;
        }

        MutationInfo<ParamTYPE> mutinfo;

        double covrate;
        double mutrate;
        double SP;
        double tolerance;
        double recombination_ratio;

        int elitpop;
        //int matsize; // set to popsize when ga is constructed, maybe change by ga.matsize = ... after constructor and before ga.run()
        int tntsize;
        int genstep;
        int precision;

        std::vector<double>(*Objective)(const std::vector<ParamTYPE>&);
        void(*Selection)(Population<ParamTYPE>&);
        // void(*CrossOver)(const Population<ParamTYPE>&, CHR<ParamTYPE>&, CHR<ParamTYPE>&);
        void(*CrossOver)(const Population<ParamTYPE>&, std::vector< CHR<ParamTYPE> >&);
        bool isMultiCrossover = false;
        void(*Mutation)(CHR<ParamTYPE>&);
        void(*Adaptation)(Population<ParamTYPE>&) = nullptr;
        std::vector<double>(*Constraint)(const std::vector<ParamTYPE>&);
        void(*FixedValue)(Population<ParamTYPE>&, int k);
        bool(*StopCondition)(galgo::GeneticAlgorithm<ParamTYPE>&);

        std::vector<bool>       force_value_flag;
        std::vector<ParamTYPE>  force_value;
        std::string csvFileName;

        int nbgen;
        int popsize;
        bool output;
    };
}
#endif
