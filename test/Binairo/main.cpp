//--------------------------------------------------
//  Binairio
//
//  Copyright (C) 2018 Samuel Lanthier, septembre 2018
//  License: MIT License 
//--------------------------------------------------
#include "MatriceUtil.h"
#include "Algorithm.h"
#include <conio.h>

// Description Grille Binairo:
//  Chacune des lignes et chacune des colonnes doit contenir autant de 0 que de 1. 
//  Ceci signifie donc qu’il y aura 5 cases contenant la valeur 0 et 5 cases contenant la valeur 1 pour chaque ligne et chaque colonne.
//  Il est interdit d’avoir une succession de plus de 2 chiffres identiques, que ce soit à l’horizontale ou à la verticale.
//  On ne peut donc pas avoir 3 cases consécutives contenant la valeur 0 ni trois cases consécutives contenant la valeur 1.
//  Dans la grille, chaque ligne doit être unique, i.e.deux lignes ne peuvent pas avoir la même séquence de 1 et de 0.
//  Dans la grille, chaque colonne doit être unique, i.e.deux colonnes ne peuvent pas avoir la même séquence de 1 et de 0.
//  Chaque grille n’a qu’une seule solution

// N: Nombre de colonne/rangee pour une grille Binairo
const int N = 10;

const std::string str_menu_continuer = "Presser une touche pour continuer (x pour terminer) ";

// Fonction pour continuer ou terminer le menu
bool menu_continuer()
{
    std::cout << str_menu_continuer;
    int ch = _getch();
    ch = toupper(ch);
    if (ch == 'X') return true;
    else return false;
}

// Fonction pour lire une grille Binairo (a resoudre) dans une matrice
template<class T>
void LireFichier(std::ifstream& in, MatriceUtil<T>& m)
{
    int x = 0;
    while (!in.eof())
    {
        int y = 0;
        std::string line;
        std::getline(in, line);
        for (char ch : line)
        {
            if (ch == '0') m.set(x, y, 0);
            else if (ch == '1') m.set(x, y, 1);
            else m.set(x, y, -1);
            ++y;
        }
        ++x;
    }
    in.close();
}

// Fonction pour ecrire la solution d'une grille Binairo dans un fichier
template<class T>
void EcrireFichier(MatriceUtil<T>& m, std::string filename)
{
    std::vector<size_t> n_verif_col(m.size_col(), 0);

    std::ofstream output(std::string("Solution-") + filename);
    if (output.is_open() == true)
    {
        for (int i = 0; i < (int)m.size_row(); ++i) 
        {
            size_t n_verif_row = 0;

            for (int j = 0; j < (int)m.size_col(); ++j)
            {
                int v = m.get(i, j);
                if (v == 0)
                {
                    output << " 0  ";
                    n_verif_row = (2 * n_verif_row) + 0;
                    n_verif_col[j] = (2 * n_verif_col[j]) + 0;
                }
                else if (v == 1)
                {
                    output << " 1  ";
                    n_verif_row = (2 * n_verif_row) + 1;
                    n_verif_col[j] = (2 * n_verif_col[j]) + 1;
                }
            }
            output << " " << n_verif_row << std::endl;
        }

        for (size_t j = 0; j < m.size_col(); j++)
        {
            output << " " << n_verif_col[j];
        }
        std::cout << std::endl;

        output.close();
    }
}

// Fonction pour resoudre en batch toutes les grilles
void resoudre_batch()
{
    std::vector<std::string> v = {
        "ct01.txt",
        "ct02.txt",
        "ct03.txt",
        "ct04.txt",
        "ct05.txt",
        "ct06.txt",
        "ct07.txt",
        "ct08.txt",
        "ct09.txt",
        "ct10.txt",
        "ct11.txt",
        "ct12.txt",
        "ct13.txt",
        "ct14.txt"/*,
        "hard.txt"*/
    };

    for (auto& f : v)
    {
        MatriceUtil<int> mat(N, N);

        std::ifstream in(f);

        // Valider fichier existe
        if (in.is_open() == false)
        {
            continue;
        }

        // Lire fichier (grille)
        LireFichier(in, mat);

        if (resolve_binairio<int>(mat) == true)
        {
            // Afficher solution grille
            std::cout << "Solution grille :" << f  << std::endl;
            display_binairio<int>(mat, true);

            // Ecrire solution grille
            EcrireFichier<int>(mat, f);
        }
        else
        {
            std::cout << "Pas de solution trouvée pour :" << f << std::endl;
        }
    }
}

// Point d' entrer du programme principale
void main(int argc, char* argv[])
{
    // Mode batch
    if (argc == 2)
    {
        resoudre_batch();
        return;
    }

    // mat: Grille Binairo a resoudre
    MatriceUtil<int> mat(N, N);

    // Done est indicateur de fin de programme
    bool done = false;

    // Afficher menu tant usager n' a pas terminer
    while (done == false)
    {
        system("CLS");
        std::cout << "Entrez le nom du fichier contenant la grille (Ex: ct01.txt) : ";

        // filename contient nom du fichier (grille) a lire
        std::string filename;

        // Saisir le nom du fichier a lire
        std::cin >> filename;
        std::cout << std::endl;

        // Valider le nom du fichier
        if (filename.empty())
        {
            std::cout << "Nom fichier vide - recommencer" << std::endl;
            done = menu_continuer();
            continue;
        }

        // Ouvrir fichier
        std::ifstream in(filename);

        // Valider fichier existe
        if (in.is_open() == false)
        {
            std::cout << "Fichier inexistant - recommencer" << std::endl;
            done = menu_continuer();
            continue;
        }

        // Lire fichier (grille)
        LireFichier(in, mat);

        // Valider grille lue
        if (mat.size_row() != mat.size_col())
        {
            std::cout << "Nombre de colonnes pas identique au nombre de lignes  - recommencer" << std::endl;
            done = menu_continuer();
            continue;
        }

        // Valider grille lue
        if (mat.size_row() != N)
        {
            std::cout << "Nombre de lignes doit etre " << N << " - recommencer" << std::endl;
            done = menu_continuer();
            continue;
        }

        // Afficher grille
        std::cout << std::endl;
        std::cout << "Grille initiale:" << std::endl;
        display_binairio<int>(mat);

        // Resoudre grille
        if (resolve_binairio<int>(mat) == true)
        {
            // Afficher solution grille
            std::cout << "Solution grille:" << std::endl;
            display_binairio<int>(mat, true);

            // Ecrire solution grille
            EcrireFichier<int>(mat, filename);
        }
        else
        {
            std::cout << "Pas de solution trouvée" << std::endl;
        }

        // Demander a l'usager de continuer ou non
        done = menu_continuer();
    }

}