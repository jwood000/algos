/*! 
 *  \file rsafactorize.cc
 *  \brief C function that transfers input from R to 
 *          quadraticSieve function for factoring large
 *            numbers and returning result to R console
 *
 *  \version 1
 *
 *  \date Created: 10/06/17
 *  \date Last modified: Time-stamp: <2017-10-06 12:33:33 EDT jwood000>
 *
 *  \author Joseph Wood
 *
 *  \note Licence: GPL (>=) 2  
 */

#include <Rcpp.h>
#include "Rgmp.h"
#include "quadraticsieve.h"
#include "rsafactorize.h"

using namespace Rcpp;

typedef std::vector<signed long int> v1d;
typedef std::vector<v1d> v2d;

SEXP QuadraticSieveContainer (SEXP m, SEXP m2, SEXP n, SEXP FB) {
    bigvec nBig = bigintegerR::create_bignum(n);
    bigvec m2Big = bigintegerR::create_bignum(m2);
    mpz_t nmpz;
    unsigned long int i, mySize = m2Big.size();
    v1d myMat, freeV = Rcpp::as<v1d>(FB);
    v2d mat = Rcpp::as<v2d>(m);
    mpz_t big2mpz[mySize];
    for (i = 0; i < mySize; i++) {mpz_init_set(big2mpz[i], m2Big[i].value.getValue());}
    
    mpz_init_set(nmpz, nBig[0].value.getValue());
    bigvec result;
    
    quadraticSieve (nmpz, 0, 0, 0, result);
    return bigintegerR::create_SEXP(result);
}