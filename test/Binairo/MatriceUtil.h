//--------------------------------------------------
// Classe MatriceUtil deriver de Matrice
// Ajout de fonctions utilitaires
//
// Copyright (C) 2018 Samuel Lanthier, septembre 2018
// License: MIT License 
//--------------------------------------------------
#pragma once
#include "Matrice.h"

template <class T>
class MatriceUtil : public Matrice<T>
{
public:
    MatriceUtil() : Matrice<T>() { n = 0; m = 0; }
    MatriceUtil(size_t r, size_t c) : Matrice<T>(r, c) {}
    MatriceUtil(const MatriceUtil & ref) : Matrice<T>(ref) {}

    //-------------------------------------------------------------
    // Ajout de fonctions utilitaires
    //-------------------------------------------------------------
    inline void checkIndexes(size_t i, size_t j) const
    {
        if (i > size_row() || j > size_col()) {
            throw new std::runtime_error("The indexes is outside the valid range");
        }
    }

    void set(T x, T y, T value) { _v_.at(x)[y] = value; }

    size_t size_row()  const  { return (size_t)n; }
    size_t size_col()  const  { return (size_t)m; }

    void fill(T value)
    {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                _v_[i][j] = value;
            }
        }
    }

    int count(T value)
    {
        int sum = 0;
        for (size_t i = 0; i < size_row(); ++i) {
            for (size_t j = 0; j < size_col(); ++j) {
                if (_v_[i][j] == value)
                    sum++;
            }
        }
        return sum;
    }

    std::vector<std::pair<int, int>> get_cells(T value)
    {
        std::vector<std::pair<int, int>> v;
        for (size_t i = 0; i < size_row(); i++)
        {
            for (size_t j = 0; j < size_col(); ++j)
            {
                if (_v_[i][j] == value)
                    v.push_back(std::pair<int, int>{(int)i, (int)j});
            }
        }
        return v;
    }

    void set(size_t i, size_t j, const T value)
    {
        checkIndexes(i, j);
        _v_[i][j] = value;
    }

    T get(size_t i, size_t j) const
    {
        checkIndexes(i, j);
        return _v_[i][j];
    }

    // fill other 0/1 when row/col has all 0 or 1 already needed
    void fill_forced()
    {
        int n_half = (int)(size_row() / 2);
        for (size_t r = 0; r < n; r++) // per row
        {
            int cnt_0 = count_row(r, 0);
            if (cnt_0 == n_half)
            {
                for (size_t c = 0; c < m; c++) // per col
                {
                    if (_v_[r][c] == -1)
                    {
                        _v_[r][c] = 1; // forced 1
                    }
                }
            }

            int cnt_1 = count_row(r, 1);
            if (cnt_1 == n_half)
            {
                for (size_t c = 0; c < m; c++) // per col
                {
                    if (_v_[r][c] == -1)
                    {
                        _v_[r][c] = 0; // forced 0
                    }
                }
            }
        }

        for (size_t c = 0; c < m; c++) // per col
        {
            int cnt_0 = count_col(c, 0);
            if (cnt_0 == n_half)
            {
                for (size_t r = 0; r < n; r++) // per row
                {
                    if (_v_[r][c] == -1)
                    {
                        _v_[r][c] = 1; // forced 1
                    }
                }
            }

            int cnt_1 = count_col(c, 1);
            if (cnt_1 == n_half)
            {
                for (size_t r = 0; r < n; r++) // per row
                {
                    if (_v_[r][c] == -1)
                    {
                        _v_[r][c] = 0; // forced 0
                    }
                }
            }
        }

    }

    int count_row(size_t row, T value)
    {
        if (row >= n)
        {
            throw new std::runtime_error("The indexes is outside the valid range");
        }

        int sum = 0;
        for (size_t j = 0; j < m; j++)
        {
            if (_v_[row][j] == value)
                sum++;
        }
        return sum;
    }

    int count_col(size_t col, T value)
    {
        if (col >= m)
        {
            throw new std::runtime_error("The indexes is outside the valid range");
        }

        int sum = 0;
        for (size_t i = 0; i < n; i++)
        {
            if (_v_[i][col] == value)
                sum++;
        }
        return sum;
    }

    int row_max_sequence(size_t row, T value)
    {
        if (row >= n)
        {
            throw new std::runtime_error("The indexes is outside the valid range");
        }

        int sum = 0;
        int max = 0;
        for (size_t j = 0; j < m; j++)
        {
            if (_v_[row][j] == value)
                sum++;
            else
                sum = 0;

            if (sum > max)
                max = sum;
        }

        return max;
    }

    int col_max_sequence(size_t col, T value)
    {
        if (col >= m)
        {
            throw new std::runtime_error("The indexes is outside the valid range");
        }

        int sum = 0;
        int max = 0;
        for (size_t i = 0; i < n; i++)
        {
            if (_v_[i][col] == value)
                sum++;
            else
                sum = 0;

            if (sum > max)
                max = sum;
        }
        return max;
    }

    bool row_same(size_t row1, size_t row2)
    {
        if (row1 >= n || row2 >= n)
        {
            throw new std::runtime_error("The indexes is outside the valid range");
        }

        for (size_t j = 0; j < m; j++)
        {
            if (_v_[row1][j] != _v_[row2][j])
                return false;
        }

        return true;
    }

    bool col_same(size_t col1, size_t col2)
    {
        if (col1 >= n || col2 >= n)
        {
            throw new std::runtime_error("The indexes is outside the valid range");
        }

        for (size_t i = 0; i < n; i++)
        {
            if (_v_[i][col1] != _v_[i][col2])
                return false;
        }

        return true;
    }
};
