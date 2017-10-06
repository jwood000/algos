/*! 
 *  ******* Original factor.cc comments (pertaining to factorR) *******
 *  \file factor.cc
 *  \brief C function used for factorization
 *
 *  \version 1
 *
 *  \date Created: 04/12/04
 *  \date Last modified: Time-stamp: <2013-06-08 22:47:18 antoine>
 *
 *  \author Antoine Lucas (help from Immanuel Scholz) (R adaptation)
 *          Original C code from libgmp.
 *
 *  \note Licence: GPL
 *  ************** End original comments *****************************
 *  
 *  \date Last modified: Time-stamp: <2017-10-02 17:19:50 jwood000>
 *  
 *  \author Joseph Wood
 *  
 */

#include "Rgmp.h"
#include "quadraticsieve.h"
#include "rsafactorize.h"

using namespace std;

SEXP QuadraticSieveContainer (SEXP n) {
    bigvec nBig = bigintegerR::create_bignum(n);
    mpz_t nmpz;
    mpz_init_set(nmpz, nBig[0].value.getValue());
    bigvec result;

    quadraticSieve (nmpz, 0, 0, 40000);

    return bigintegerR::create_SEXP(result);
}