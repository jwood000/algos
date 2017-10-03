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
#include "factorize.h"
#include "factor.h"

using namespace std;

// Exact same function obtained from the gmp package
SEXP factorR (SEXP n)
{
    bigvec v = bigintegerR::create_bignum(n), result;
    if(v.size() > 0) {
        mpz_t val;
        mpz_init(val);
        mpz_t_sentry val_s(val);
        mpz_set(val,v[0].value.getValueTemp());

        int sgn = mpz_sgn(val);
        if(sgn == 0)
            error(_("Cannot factorize 0"));
        if(sgn<0)
        {
            mpz_abs(val,val);
            result.value.push_back(biginteger(-1));
        }
        //
        // function from gmplib, in demo/factorize.c
        //
        factor(val,result);
    }
    return bigintegerR::create_SEXP(result);
}

SEXP getDivisorsC (SEXP n)
{
  bigvec v = bigintegerR::create_bignum(n);
  bigvec result;
  if(v.size() > 0) {
    if (mpz_cmp_ui(v[0].value.getValueTemp(), 1) == 0) {
        result.push_back(1);
    } else {
        mpz_t val;
        mpz_init(val);
        mpz_t_sentry val_s(val);
        mpz_set(val,v[0].value.getValueTemp());
        
        int sgn = mpz_sgn(val);
        if(sgn == 0)
            error(_("Cannot factorize 0"));
        if(sgn<0)
        {
            mpz_abs(val,val);
            result.value.push_back(biginteger(-1));
        }
        allDivisors (val,result);
    }
  }
  return bigintegerR::create_SEXP(result);
}

SEXP QuadraticSieveContainer (SEXP n) {
    bigvec nBig = bigintegerR::create_bignum(n);
    mpz_t nmpz;
    mpz_init_set(nmpz, nBig[0].value.getValue());
    bigvec result;

    QuadraticSieve (nmpz, 0, 0, 40000);

    return bigintegerR::create_SEXP(result);
}

SEXP QuadraticResidueContainer (SEXP n, SEXP p) {
    mpz_t nmpz, pmpz;
    bigvec v1 = bigintegerR::create_bignum(n);
    bigvec v2 = bigintegerR::create_bignum(p);
    mpz_init(nmpz);
    mpz_init(pmpz);
    mpz_set(nmpz,v1[0].value.getValueTemp());
    mpz_set(pmpz,v2[0].value.getValueTemp());
    bigvec result;
    
    TonelliShanksC(nmpz, pmpz, result);
    
    return bigintegerR::create_SEXP(result);
}
