//--------------------------------------------------
//  Binairio
//
//  Copyright (C) 2018 Samuel Lanthier, septembre 2018
//  License: MIT License 
//--------------------------------------------------
#pragma once
#include <vector>
#include <exception>

template <class T>
class Matrice
{
public:
    Matrice()  : n(0),  m(0) {}

	Matrice(size_t r, size_t c ) : n(r), m(c), _v_(r)
	{
		for(size_t i = 0; i < r; i++ )
			_v_[ i ].resize( c );
	}

    Matrice(const Matrice & ref) : n(ref.n), m(ref.m), _v_(ref._v_) {}

	const std::vector<T> & operator[]( int r ) const { return _v_[ r ]; }
    std::vector<T> & operator[]( int r )             { return _v_[ r ]; }
	const std::vector<T> & at (int r) const          {	return _v_.at( r );	}
	std::vector<T> & at (int r)                      {	return _v_.at( r );	}

protected:
    size_t n;   // row size
    size_t m;   // col size
    std::vector< std::vector<T> > _v_;
};
