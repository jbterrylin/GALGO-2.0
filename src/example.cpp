//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#include "Galgo.hpp"

//------------------------------------------------------------------------------
// Uncomment #define TEST_ALL_TYPE to test compiling of all types
// Uncomment #define TEST_BINAIRO to test GA for Binairos
// Uncomment #define TEST_CLASSIC_FUNCTIONS to test GA for classics functions
//------------------------------------------------------------------------------
//#define TEST_ALL_TYPE
//#define TEST_BINAIRO
#define TEST_CLASSIC_FUNCTIONS

#ifdef TEST_CLASSIC_FUNCTIONS
#include "..\test\Classic\Functions.hpp"
#endif

#ifdef TEST_ALL_TYPE
#include "..\test\Types\TestTypes.hpp"
#endif

#ifdef TEST_BINAIRO
#include "..\test\Binairo\GA.h"
#endif


int main()
{  
#ifdef TEST_ALL_TYPE
    TEST_all_types();
#endif

#ifdef TEST_BINAIRO
    GA_Binairo::test_ga_binairo(4);     // 0=resolve one free cell(hard), 1=resolve 4 free cells(very hard), 2=resolve 7 free cells(diabolical), 3 , 4==generate new grid
#endif

#ifdef TEST_CLASSIC_FUNCTIONS
    test_classic();
#endif

    system("pause");
    return 0;
}

