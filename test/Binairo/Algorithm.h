//--------------------------------------------------
//  Fichier: Algorithm.h
//  Contient fonctions utilitaires pour grille Binairo
//
//  Copyright (C) 2018 Samuel Lanthier, septembre 2018
//  License: MIT License 
//--------------------------------------------------
#pragma once

#include "MatriceUtil.h"
#include <cstdlib> 
#include <iostream> 
#include <ctime>
#include <time.h>
#include <vector>
#include <string>
#include <fstream>

// Description: Fonction pour afficher une grille Binairo (matrice)
// Parametre m : Matrice a utiliser
// Parametre show_verification : Si true afficher les nombres de verification
//
template<class T>
void display_binairio(MatriceUtil<T>& m, bool show_verification = false)
{
    std::vector<size_t> n_verif_col(m.size_col(), 0);

    std::cout << "    ";
    for (size_t j = 0; j < m.size_col(); j++)
    {
        std::cout << " " << j << " |";
    }

    if (show_verification)
        std::cout << " Verification";
    std::cout << std::endl;

    std::cout << "----";
    for (size_t j = 0; j < m.size_col(); j++)
        std::cout << "----";
    std::cout << std::endl;

    for (size_t i = 0; i < m.size_row(); i++)
    {
        std::cout << " " << i << " |";

        size_t n_verif_row = 0;
        for (size_t j = 0; j < m.size_col(); j++)
        {
            T v = m.get(i, j);
            if (v == -1)
            {
                std::cout << " * ";
            }
            else if (v == 0)
            {
                std::cout << " 0 ";
                n_verif_row = (2 * n_verif_row) + 0;
                n_verif_col[j] = (2 * n_verif_col[j]) + 0;
            }
            else if (v == 1)
            {
                std::cout << " 1 ";
                n_verif_row = (2 * n_verif_row) + 1;
                n_verif_col[j] = (2 * n_verif_col[j]) + 1;
            }
            std::cout << "|";
        }

        if (show_verification)
            std::cout << " " << n_verif_row;
        std::cout << std::endl;

        std::cout << "----";
        for (size_t j = 0; j < m.size_col(); j++)
            std::cout << "----";
        std::cout << std::endl;
    }

    if (show_verification)
    {
        std::cout << "   ";
        for (size_t j = 0; j < m.size_col(); j++)
        {
            std::cout << " " << n_verif_col[j];
        }
    }

    std::cout << std::endl;
}

// Description: Fonction pour valider si les cellules vides d'une grille peuvent etre 0 ou 1
// Parametre m_work : Matrice a utiliser
// Retourne bool - true si grille demeure valide en essayant des 0 ou 1 dans cellules vides, false sinon
//
template<class T>
bool valid_empty_cells(MatriceUtil<T> m_work)
{
    // Check if empty cell can still get 0 or 1
    bool ok = true;
    MatriceUtil<T> m_test_0(m_work.size_row(), m_work.size_col());
    MatriceUtil<T> m_test_1(m_work.size_row(), m_work.size_col());

    for (size_t i = 0; i < m_work.size_row(); i++)
    {
        if (ok == false)
            break;

        for (size_t j = 0; j < m_work.size_col(); j++)
        {
            if (m_work.get(i, j) == -1)
            {
                m_test_0 = m_work;
                m_test_1 = m_work;

                m_test_0.set(i, j, 0);
                if (valid_binairio<T>(m_test_0, false) == true)
                {
                    // ok
                }
                else
                {
                    m_test_1.set(i, j, 1);
                    if (valid_binairio<T>(m_test_1, false) == true)
                    {
                        // ok 
                    }
                    else
                    {
                        // invalid
                        ok = false;
                        break;
                    }
                }
            }
        }
    }
    return ok;
}

// Description: Fonction pour valider une grille Binairo
// Parametre m : Matrice a utiliser
// Parametre check_row_col_same : Si false, il s' agit d' une grille incomplete a valider, sinon grille complete a valider
// Retourne bool - true si grille valide, false sinon
//
template<class T>
bool valid_binairio(MatriceUtil<T> m, bool check_row_col_same = true)
{
    for (size_t i = 0; i < m.size_row(); i++)
    {
        if (m.count_row(i, 0) > (int) (m.size_col() / 2))
            return false;

        if (m.count_row(i, 1) > (int)(m.size_col() / 2))
            return false;

        if (check_row_col_same)
        {
            if (m.count_row(i, 0) != m.count_row(i, 1))
                return false;
        }

        if (m.row_max_sequence(i, 0) > 2)
            return false;

        if (m.row_max_sequence(i, 1) > 2)
            return false;

        if (check_row_col_same)
        {
            for (size_t j = 0; j < m.size_row(); j++)
            {
                if (i < j)
                {
                    if (m.row_same(i, j) == true)
                        return false;
                }
            }
        }
    }

    for (size_t i = 0; i < m.size_col(); i++)
    {
        if (m.count_col(i, 0) >(int) (m.size_row() / 2))
            return false;

        if (m.count_col(i, 1) > (int)(m.size_row() / 2))
            return false;

        if (check_row_col_same)
        {
            if (m.count_col(i, 0) != m.count_col(i, 1))
                return false;
        }

        if (m.col_max_sequence(i, 0) > 2)
            return false;

        if (m.col_max_sequence(i, 1) > 2)
            return false;

        if (check_row_col_same)
        {
            for (size_t j = 0; j < m.size_col(); j++)
            {
                if (i < j)
                {
                    if (m.col_same(i, j) == true)
                        return false;
                }
            }
        }
    }

    return true;
}

// Description: Fonction pour resoudre une grille Binairo en essayant seulement des 0 ou 1
// Parametre(input/output) m : Matrice a utiliser et a remplir avec valeurs (0 ou 1) valides
// Parametre (output) is_valid : Si true, la grille est encore valide apres avoir mis des 0 ou 1, sinon false
// Retourne bool - true si grille valide et complete (solutionner), false sinon
//
template<class T>
bool try_resolve_binairio(MatriceUtil<T>& m, bool& is_valid)
{
    is_valid = false;
    MatriceUtil<T> m_work = m;
    bool filling_atleast_one_cell = true;   // Indicateur pour continuer a placer des 0 ou 1

    while (filling_atleast_one_cell)
    {
        // Placer valeurs forcer
        m_work.fill_forced();

        filling_atleast_one_cell = false;
        for (size_t i = 0; i < m_work.size_row(); i++)
        {
            for (size_t j = 0; j < m_work.size_col(); j++)
            {
                int force_value = -1;

                // Verifier si 0 ou 1 est forcer sur cellule vide
                if (m_work.get(i, j) == -1)
                {
                    MatriceUtil<T> m_test_0 = m_work;
                    MatriceUtil<T> m_test_1 = m_work;

                    // Placer 0 et valeurs forcer
                    m_test_0.set(i, j, 0);
                    m_test_0.fill_forced();

                    if (valid_binairio<T>(m_test_0, false) == true)
                    {
                        // Placer 1
                        m_test_1.set(i, j, 1);
                        m_test_1.fill_forced();
                        if (valid_binairio<T>(m_test_1, false) == false)
                        {
                            // ok - seulement 0 est possible
                            force_value = 0;
                            filling_atleast_one_cell = true; // continuer a placer
                        }
                    }
                    else
                    {
                        // 0 est pas possible
                        // Placer 1
                        m_test_1.set(i, j, 1);
                        m_test_1.fill_forced();
                        if (valid_binairio<T>(m_test_1, false) == true)
                        {
                            // ok - seulement 1 est possible
                            force_value = 1;
                            filling_atleast_one_cell = true; // continuer a placer
                        }
                    }

                    if (force_value == 0) m_work = m_test_0;
                    else if (force_value == 1) m_work = m_test_1;
                }
            }
        }
    }

    bool is_full = true;
    for (size_t i = 0; i < m_work.size_row(); i++)
    {
        if (is_full == false)
            break;

        for (size_t j = 0; j < m_work.size_col(); j++)
        {
            if (m_work.get(i, j) == -1)
            {
                is_full = false;
                break;
            }
        }
    }

    if (is_full == true)
    {
        if (valid_binairio<T>(m_work) == true)
        {
            m = m_work;
            return true;    // La grille est complete et valide (donc solutionner)
        }
        else
        {
            // La grille est complete mais pas valide
        }
    }
    else
    {
        if (valid_binairio<T>(m_work, false) == true)
        {
            // check if empty cell can still get 0 or 1
            bool ok = valid_empty_cells<T>(m_work);

            if (ok == true)
            {
                is_valid = true;
            }
        }
    }

    if (is_valid)
    {
        m = m_work;
    }

    return false;
}

// Description: Fonction pour resoudre completement une grille Binairo
// Parametre(input/output) m : Matrice a utiliser et a remplir avec valeurs (0 ou 1) valides
// Retourne bool - true si grille solutionner completement, false sinon
//
template<class T>
bool resolve_binairio(MatriceUtil<T>& m)
{
    // Esayer de resoudre seulement en mettant des 0 et 1 dans les cellules vides
    bool is_valid;
    if (try_resolve_binairio(m, is_valid) == true)
    {
        // On a trouver la solution
        return true;
    }

    std::cout << "Grille courante:" << std::endl;
    display_binairio<T>(m, false);

    // La grille est trop complexe pour etre resolu seulement en mettant des 0 et 1 dans les cellules vides
    // Il faut essayer 2 hypotheses (0 ou 1) par cellule vide et voir si la grille se resout alors

    // cnt_empty_cells: Nombre de cellules vides restantes
    size_t cnt_empty_cells = m.count(-1);

    // Essayer 0 ou 1 dans les cellules vides et voir si on peux resoudre la grille alors
    if ((is_valid == true) && (cnt_empty_cells > 0))
    {
        std::vector<std::pair<int, int>> v_cell = m.get_cells(-1);
        std::vector< MatriceUtil<T>> v_one;

        std::cout << "Grille complexe. Creation de scenarios pour " << v_cell.size() << " cellules vides." << std::endl;

        // Essayer 0 dans une cellule vide puis tenter de resoudre
        for (size_t c = 0; c < v_cell.size(); c++)
        {
            MatriceUtil<T> m_work = m;
            m_work.set(v_cell[c].first, v_cell[c].second, 0);
            if (try_resolve_binairio(m_work, is_valid) == true)
            {
                m = m_work;
                return true;
            }

            if (is_valid)
            {
                v_one.push_back(m_work);
            }
        }

        // Essayer 1 dans une cellule vide puis tenter de resoudre
        for (size_t c = 0; c < v_cell.size(); c++)
        {
            MatriceUtil<T> m_work = m;
            m_work.set(v_cell[c].first, v_cell[c].second, 1);
            if (try_resolve_binairio(m_work, is_valid) == true)
            {
                m = m_work;
                return true;
            }

            if (is_valid)
            {
                v_one.push_back(m_work);
            }
        }

        // On a maintenant v_one.size() hypotheses (grilles) a verifier
        // Chaque hypothese de v_one est une grille partiellement completer (et valide) avec une cellule en hypothese
        // Valider chacune des ces hypotheses 

        std::cout << "Nombre de scenarios (avec 1 cellule en hypothese) en cours de resolution: " << v_one.size() << std::endl;

        std::vector< MatriceUtil<T>> v_two;
        for (size_t k = 0; k < v_one.size(); k++)
        {
            std::vector<std::pair<int, int>> v_cell = v_one[k].get_cells(-1);

            // Essayer 0 dans une cellule vide puis tenter de resoudre
            for (size_t c = 0; c < v_cell.size(); c++)
            {
                MatriceUtil<T> m_work = v_one[k];
                m_work.set(v_cell[c].first, v_cell[c].second, 0);
                if (try_resolve_binairio(m_work, is_valid) == true)
                {
                    m = m_work;
                    return true;
                }

                if (is_valid)
                {
                    v_two.push_back(m_work);
                }
            }

            // Essayer 1 dans une cellule vide puis tenter de resoudre
            for (size_t c = 0; c < v_cell.size(); c++)
            {
                MatriceUtil<T> m_work = v_one[k];
                m_work.set(v_cell[c].first, v_cell[c].second, 1);
                if (try_resolve_binairio(m_work, is_valid) == true)
                {
                    m = m_work;
                    return true;
                }

                if (is_valid)
                {
                    v_two.push_back(m_work);
                }
            }
        }

        // On a maintenant v_two.size() hypotheses (grilles) a verifier
        // Chaque hypothese de v_two est une grille partiellement completer (et valide) avec deux cellules en hypothese
        // Valider chacune des ces hypotheses 

        std::cout << "Nombre de scenarios (avec 2 cellules en hypothese) en cours de resolution: " << v_two.size() << std::endl;

        std::vector< MatriceUtil<T>> v_three;
        for (size_t k = 0; k < v_two.size(); k++)
        {
            if (k % 100 == 0)
            {
                std::cout << "Scenario: " << k+1 <<  " en cours de resolution..." << std::endl;
            }
            std::vector<std::pair<int, int>> v_cell = v_two[k].get_cells(-1);

            // Essayer 0 dans une cellule vide puis tenter de resoudre
            for (size_t c = 0; c < v_cell.size(); c++)
            {
                MatriceUtil<T> m_work = v_two[k];
                m_work.set(v_cell[c].first, v_cell[c].second, 0);
                if (try_resolve_binairio(m_work, is_valid) == true)
                {
                    m = m_work;
                    return true;
                }

                if (is_valid)
                {
                    v_three.push_back(m_work);
                }
            }

            // Essayer 1 dans une cellule vide puis tenter de resoudre
            for (size_t c = 0; c < v_cell.size(); c++)
            {
                MatriceUtil<T> m_work = v_two[k];
                m_work.set(v_cell[c].first, v_cell[c].second, 1);
                if (try_resolve_binairio(m_work, is_valid) == true)
                {
                    m = m_work;
                    return true;
                }

                if (is_valid)
                {
                    v_three.push_back(m_work);
                }
            }
        }

        // On a maintenant v_three.size() hypotheses (grilles) a verifier
        // Chaque hypothese de v_three est une grille partiellement completer (et valide) avec trois cellules en hypothese

        std::cout << "Nombre de scenarios (3 cellules) en cours de resolution: " << v_three.size() << std::endl;

        // Abandonner - grille trop complexe a verifier (trop d'hypotheses a valider) ou invalide
        std::cout << "Abandon - grille trop complexe a verifier (trop d'hypotheses a valider) ou invalide" << std::endl;
    }

    return false;
}

