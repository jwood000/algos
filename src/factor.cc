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
  bigvec result, myFacs, bigFacs, tempBig;
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

    std::vector<int> lengths;
    std::vector<biginteger>::iterator it;
    biginteger prev = result.value[0];
    int n = result.size();
    lengths.reserve(n);
    bigFacs.value.reserve(n);
    int numUni = 0;
    bigFacs.value.push_back(prev);
    lengths.push_back(1);
    for(it = result.value.begin() + 1; it < result.value.end(); it++) {
        if (prev == *it) {
            lengths[numUni]++;
        } else {
            bigFacs.value.push_back(*it);
            lengths.push_back(1);
            numUni++;
            prev = *it;
        }
    }

    long i, j, k;
    long mySize, facSize;
    myFacs.value.reserve(lengths[0]+1);
    mpz_t temp;
    mpz_init (temp);

    for (i = 0; i <= lengths[0]; ++i) {
        mpz_pow_ui(temp, bigFacs.value[0].getValue(), i);
        myFacs.value.push_back(temp);
    }

    if (numUni > 0) {
        for (j = 1; j <= numUni; ++j) {
            facSize = myFacs.size();
            mySize = facSize * lengths[j];
            tempBig.value.reserve(mySize);
            for (i = 1; i <= lengths[j]; ++i) {
                for (k = 0; k < facSize; k++) {
                    mpz_pow_ui(temp, bigFacs.value[j].getValue(), i);
                    mpz_mul(temp, temp, myFacs[k].value.getValue());
                    tempBig.push_back(temp);
                }
            }
            myFacs.value.reserve(myFacs.value.size() + mySize);
            myFacs.value.insert(myFacs.value.end(),
                                tempBig.value.begin(),tempBig.value.end());
            tempBig.clear();
        }
    }

    std::sort(myFacs.value.begin(), myFacs.value.end());
  }
  return bigintegerR::create_SEXP(myFacs);
}

SEXP QuadraticSieveContainer (SEXP n) {
    bigvec nBig = bigintegerR::create_bignum(n);
    mpz_t nmpz;
    mpz_init_set(nmpz, nBig[0].value.getValue());
    bigvec result;

    QuadraticSieve (nmpz, 0, 0, 40000);

    return bigintegerR::create_SEXP(result);
}
