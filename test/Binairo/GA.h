//--------------------------------------------------
//  Fichier: GA.h
//
//  Copyright (C) 2018 Alain Lanthier, septembre 2018
//  License: MIT License 
//--------------------------------------------------
#pragma once

#include "MatriceUtil.h"
#include "Algorithm.h"

using BINAIRO_TEST_TYPE = int;
// 4x4
//static std::vector<BINAIRO_TEST_TYPE> binairo_initial  = { 1,0,-1,0, -1,-1,1,-1, -1,-1,0,0, -1,-1,-1,-1 };
//static std::vector<BINAIRO_TEST_TYPE> binairo_solution = { 1,0, 1,0,  0, 0,1, 1,  1, 1,0,0,  0, 1, 0, 1 };
// 8x8
//static std::vector<BINAIRO_TEST_TYPE> binairo_initial  = {-1,-1,-1,-1,-1,-1,-1, 0, -1, 0, 0,-1,-1, 1,-1,-1,  -1, 0,-1,-1,-1, 1,-1, 0,  -1,-1, 1,-1,-1,-1,-1,-1,
//                                                           0, 0,-1, 1,-1,-1, 1,-1, -1,-1,-1,-1, 1,-1,-1,-1,   1, 1,-1,-1,-1, 0,-1, 1,  -1, 1,-1,-1,-1,-1,-1, 1};
//static std::vector<BINAIRO_TEST_TYPE> binairo_solution = { 0,1,1,0,1,0,1,0,  1,0,0,1,0,1,0,1,  1,0,0,1,0,1,1,0, 0,1,1,0,1,0,0,1,
//                                                           0,0,1,1,0,1,1,0,  1,0,0,1,1,0,1,0,  1,1,0,0,1,0,0,1, 0,1,1,0,0,1,0,1};
static std::vector<BINAIRO_TEST_TYPE> binairo_initial;
void make_binairo(int no = 0)
{
    std::string s;
    if (no == 2)
    {
        // hard3.txt
        s =
        std::string("*1***0**1*") +
        std::string("0*0*******") +
        std::string("******1***") +
        std::string("**1**0****") +
        std::string("0*********") +
        std::string("*******0**") +
        std::string("*******1*1") +
        std::string("**0***0***") +
        std::string("******0***") +
        std::string("****0*****");
    }

    else if (no == 1)
    {
        // hard2.txt
        s =
        std::string("*1***0**1*") +
        std::string("0*0*******") +
        std::string("******1***") +
        std::string("**1**0****") +
        std::string("0*********") +
        std::string("*******00*") +
        std::string("1******1*1") +
        std::string("**0***0***") +
        std::string("******0*1*") +
        std::string("****0*****");
    }

    else
    {
        // hard.txt
        s =
        std::string("*1*1*0**1*") +
        std::string("0*0*******") +
        std::string("******11**") +
        std::string("**1**0****") +
        std::string("0*********") +
        std::string("*******00*") +
        std::string("1******1*1") +
        std::string("**0***0***") +
        std::string("******0*1*") +
        std::string("****0**0**");
}

    //--------------------------------------------
    // 0 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0 |
    //--------------------------------------------
    // 1 | 0 | 1 | 0 | * | * | 1 | 0 | 0 | 1 | 1 |
    //--------------------------------------------
    // 2 | 1 | 0 | 1 | * | * | 0 | 1 | 1 | 0 | 0 |
    //--------------------------------------------
    // 3 | * | * | 1 | * | * | 0 | 1 | 0 | 0 | 1 |
    //--------------------------------------------
    // 4 | 0 | 1 | 0 | 0 | 1 | 1 | 0 | 1 | 1 | 0 |
    //--------------------------------------------
    // 5 | * | * | 1 | * | 0 | 1 | 1 | 0 | 0 | 1 |
    //--------------------------------------------
    // 6 | 1 | 0 | 0 | 1 | 0 | 0 | 1 | 1 | 0 | 1 |
    //--------------------------------------------
    // 7 | * | * | 0 | * | 1 | * | 0 | * | 1 | 0 |
    //--------------------------------------------
    // 8 | * | * | 1 | * | * | * | 0 | * | 1 | 0 |
    //--------------------------------------------
    // 9 | * | * | 1 | * | 0 | * | 1 | 0 | 0 | 1 |
    //--------------------------------------------

        //*1*1*0**1*
        //0*0*******
        //******11**
        //**1**0****
        //0*********
        //*******00*
        //1******1*1
        //**0***0***
        //******0*1*
        //****0**0**

        //0   1   0   1   1   0   0   1   1   0   358
        //0   1   0   1   0   1   0   0   1   1   339
        //1   0   1   0   1   0   1   1   0   0   684
        //1   0   1   1   0   0   1   0   0   1   713
        //0   1   0   0   1   1   0   1   1   0   310
        //0   1   1   0   0   1   1   0   0   1   409
        //1   0   0   1   0   0   1   1   0   1   589
        //1   0   0   1   1   0   0   1   1   0   614
        //0   1   1   0   1   1   0   0   1   0   434
        //1   0   1   0   0   1   1   0   0   1   665
        //205 818 211 844 678 307 217 684 806 345

    MatriceUtil<int> mat(10, 10);
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            if (s.at(10 * i + j) == '*') mat.set(i, j, -1);
            else if (s.at(10 * i + j) == '0') mat.set(i, j, 0);
            else if (s.at(10 * i + j) == '1') mat.set(i, j, 1);
        }
    }

    // Forced values
    bool is_valid;
    if (try_resolve_binairio(mat, is_valid) == true)
    {
        // SOLVED
    }

    binairo_initial = std::vector<BINAIRO_TEST_TYPE>(100);
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            binairo_initial[10 * i + j] = mat[i][j];
        }
    }
}

template <typename T>
void FixedParameterBinairo(galgo::Population<T>& x, int k)
{
    std::vector<galgo::CHR<T>>& np = x.get_newpop();

    {
        for (int j = 0; j < x.ga_algo()->force_value_flag.size(); j++)
        {
            if (x.ga_algo()->force_value_flag[j])
            {
                np[k]->initGene(j, x.ga_algo()->force_value[j]);
                if (np[k]->get_value(j) != x.ga_algo()->force_value[j])
                {
                    std::cout << "ERROR - Invalid decode/encode desired_value:" << x.ga_algo()->force_value[j] << " set_value: " << np[k]->get_value(j) << "\n";
                }
            }
        }

        size_t n = (int)pow((double)np[k]->nbgene(), 0.50);
        MatriceUtil<T> mat(n, n);

        int v;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                v = np[k]->get_value((int)n*i + j);
                mat.set(i, j, v);

                if ((v != binairo_initial[n*i + j]) && ((binairo_initial[n*i + j] == 0) || (binairo_initial[n*i + j] == 1)))
                {
                    // Put fixed parameter back into matrice
                    mat.set(i, j, binairo_initial[n*i + j]);
                    np[k]->initGene((int)n*i + j, binairo_initial[n*i + j]);
                }
            }
        }

        // Put forced parameter in matrice
        bool is_valid;
        if (try_resolve_binairio(mat, is_valid) == true)
        {
            // SOLVED
        }

        if (is_valid == true)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    np[k]->initGene((int)n*i + j, mat[i][j]);
                }
            }
        }

    }
}

template <typename T> class BinairoObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        int n = (int)pow((double)x.size(), 0.50);
        MatriceUtil<T> mat(n, n);

        bool mismatch = false;
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                mat.set(i, j, x[n*i + j]);

                if ((x[n*i + j] != binairo_initial[n*i + j]) && ((binairo_initial[n*i + j] == 0) || (binairo_initial[n*i + j] == 1)))
                {
                    mismatch = true;
                }
            }
        }

        int cnt_0 = mat.count(0);
        int cnt_1 = mat.count(1);
        int cnt_m1 = mat.count(-1);
        int cnt_other = (n*n) - (cnt_0 + cnt_1 + cnt_m1);

        int cnt_illegal_row = 0;
        int cnt_illegal_col = 0;
        int cnt_illegal_sequ = 0;
        int cnt_illegal_rowcol = 0;
        for (size_t i = 0; i < mat.size_row(); i++)
        {
            if (mat.count_row(i, 0) >(int) (n / 2))
                cnt_illegal_row++;

            if (mat.count_row(i, 1) > (int)(n/ 2))
                cnt_illegal_row++;

            if (mat.row_max_sequence(i, 0) > 2)
                cnt_illegal_sequ++;

            if (mat.row_max_sequence(i, 1) > 2)
                cnt_illegal_sequ++;

            for (size_t j = 0; j < mat.size_row(); j++)
            {
                if (i < j)
                {
                    if (mat.row_same(i, j) == true)
                        cnt_illegal_rowcol++;
                }
            }
        }

        for (size_t i = 0; i < mat.size_col(); i++)
        {
            if (mat.count_col(i, 0) >(int) (n / 2))
                cnt_illegal_col++;

            if (mat.count_col(i, 1) > (int)(n / 2))
                cnt_illegal_col++;

            if (mat.col_max_sequence(i, 0) > 2)
                cnt_illegal_sequ++;

            if (mat.col_max_sequence(i, 1) > 2)
                cnt_illegal_sequ++;

            for (size_t j = 0; j < mat.size_col(); j++)
            {
                if (i < j)
                {
                    if (mat.col_same(i, j) == true)
                        cnt_illegal_rowcol++;
                }
            }
        }

        double penality = 0.0;
        if (mismatch == true)  penality += 20000.0;
        penality += std::fabs((double)cnt_other - (0));
        penality += 50 * (double)cnt_illegal_row;
        penality += 50 * (double)cnt_illegal_col;
        penality += 50 * (double)cnt_illegal_sequ;
        penality += 50 * (double)cnt_illegal_rowcol;
        penality += 10 * std::fabs((double)cnt_m1 - (0));

        double obj = -penality;
        return { obj };
    }
};

template <typename T>
bool StopGABinairo(galgo::GeneticAlgorithm<T>& ga)
{
    const galgo::CHR<BINAIRO_TEST_TYPE> bestParam = ga.result(); //std::shared_ptr<Chromosome<T>>;
    if (bestParam->fitness >= 0.0)
    {
        std::cout << "SOLVED\n";
        if (ga.output)  ga.print(true);
        //galgo::Population<BINAIRO_TEST_TYPE>& x = ga.get_pop();
        //std::vector<galgo::CHR<T>>& np = x.get_newpop();

        size_t n = (int)pow((double)bestParam->nbgene(), 0.50);
        MatriceUtil<T> mat(n, n);

        int v;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                v = bestParam->get_value((int)n*i + j);
                mat.set(i, j, v);
            }
        }

        bool is_valid;
        if (try_resolve_binairio(mat, is_valid) == true)
        {
            // SOLVED
            display_binairio<BINAIRO_TEST_TYPE>(mat, false);
        }

        return true;
    }
    return false;
}

void test_ga_binairo(int no = 0)
{
    {
        const int NBIT = 2;
        make_binairo(no);

        std::cout << std::endl;
        std::cout << "BINAIRO grid NxN\n";
        int k = 0;
        const int NBinairo = 10;
        //[-1,0,1] 3 values/4 (2 bits) = -1=-1+00, 0=-1+01, 1=-1+10
        BINAIRO_TEST_TYPE low = -1;
        BINAIRO_TEST_TYPE high = 1;
        std::vector<BINAIRO_TEST_TYPE> vlow(NBinairo * NBinairo);
        std::vector<BINAIRO_TEST_TYPE> vhigh(NBinairo * NBinairo);
        std::vector<BINAIRO_TEST_TYPE> vinit(NBinairo * NBinairo);
        for (size_t i = 0; i < NBinairo * NBinairo; i++)
        {
            vlow[i] = low;
            vhigh[i] = high;
            vinit[i] = binairo_initial[i];
        }

        MatriceUtil<BINAIRO_TEST_TYPE> mat(NBinairo, NBinairo);
        bool mismatch = false;
        std::vector<bool> force_value_flag(NBinairo * NBinairo);
        std::vector<BINAIRO_TEST_TYPE> force_value(NBinairo * NBinairo);
        for (size_t i = 0; i < NBinairo; i++)
        {
            for (size_t j = 0; j < NBinairo; j++)
            {
                mat.set(i, j, binairo_initial[NBinairo*i + j]);
                force_value_flag[NBinairo*i + j] = false;
                force_value[NBinairo*i + j] = -1;
                if (binairo_initial[NBinairo*i + j] != -1)
                {
                    force_value_flag[NBinairo*i + j] = true;
                    force_value[NBinairo*i + j] = binairo_initial[NBinairo*i + j];
                }
            }
        }
        display_binairio<BINAIRO_TEST_TYPE>(mat, false);

        galgo::ConfigInfo<BINAIRO_TEST_TYPE> config;
        config.Objective = BinairoObjective<BINAIRO_TEST_TYPE>::Objective;

        config.mutinfo._sigma = 1.0;
        config.mutinfo._sigma_lowest = 0.01;
        config.mutinfo._ratio_boundary = 0.10;

        config.popsize  = 400;
        config.nbgen    = 10000000;
        config.output   = true;

        config.covrate = 0.10;
        config.elitpop = 50; // Keep enough single unmodified individuals
        config.recombination_ratio = 0.50;
        config.mutrate = 0.05;

        config.Selection = RWS;
        config.CrossOver = RealValuedSingleArithmeticRecombination;
        config.mutinfo._type = galgo::MutationType::MutationSPM;
        config.genstep = 50;
        config.precision = 2;

        config.force_value_flag = force_value_flag;
        config.force_value = force_value;
        config.FixedValue = FixedParameterBinairo;  // nullptr;
        config.StopCondition = StopGABinairo;

        galgo::GeneticAlgorithmN<BINAIRO_TEST_TYPE, NBIT> ga(config, vlow, vhigh, vinit);
        ga.run();

        //const galgo::CHR<BINAIRO_TEST_TYPE> bestParam = ga.result(); //std::shared_ptr<Chromosome<T>>;
        //std::vector<BINAIRO_TEST_TYPE> v = bestParam->getParam();

        system("pause");
    }
}