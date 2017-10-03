/*! \file factor.h
 *  \brief header for factor functions set
 *
 *  \version 1
 *
 *  \date Created: 2005
 *  \date Last modified: Time-stamp: <2008-02-17 21:41:29 antoine>
 *
 *
 *  \note Licence: GPL
 */

#ifndef GMP_R_FACTOR_HEADER_
#define GMP_R_FACTOR_HEADER_ 1

#include "bigintegerR.h"


extern "C"
{

  SEXP getDivisorsC (SEXP n);
  SEXP factorR (SEXP n);
  SEXP QuadraticSieveContainer (SEXP n);
  SEXP QuadraticResidueContainer (SEXP n, SEXP p);

}

void TonelliShanksC (mpz_t a, mpz_t p, bigvec & result);

void QuadraticSieve (mpz_t a, double e1, double e2, unsigned long lb);

/** \brief Function ued to test factorization with small numbers
 */
void factor_using_division (mpz_t t, unsigned int limit,  bigvec & result) ;

/** \brief Function used for factorization
 */
void factor_using_division_2kp (mpz_t t, unsigned int limit, unsigned long p,  bigvec & result) ;

/** \brief Pollard Rho method for factorization
 */
void factor_using_pollard_rho (mpz_t n, int a_int, unsigned long p, bigvec & result);

/** \brief Function that call an algorithm for factorization
 */
void factor (mpz_t t, bigvec & result);

void allDivisors (mpz_t t, bigvec & result);

#endif
