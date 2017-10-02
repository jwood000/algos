/*! \file solve.cc
 *  \brief functions to solve matrix
 *
 *  \version 1
 *
 *  \date Created: 25/05/06
 *  \date Last modified: Time-stamp: <2006-05-25 23:05:20 antoine>
 *
 *  \author A. Lucas
 *
 *  \note Licence: GPL
 */


#include "solve.h"

#include "bigrationalR.h"
#include "bigintegerR.h"

// inverse a rational matrix
SEXP inverse_q(SEXP A)
{

  bigvec_q a = bigrationalR::create_bignum(A);

  return(solve_gmp_R::inverse_q(a));
}

SEXP solve_gmp_R::inverse_q(bigvec_q a)
{
  if(a.nrow * a.nrow != a.size())
    error(_("Argument 1 must be a square matrix"));
  bigvec_q b (a.size());
  b.nrow = a.nrow;

  // initialize b to identity
  for(unsigned int i=0; i<b.nrow ; ++i)
    for(unsigned int j=0; j<b.nrow ; ++j)
	b.value[i+j*b.nrow].setValue((i == j) ? 1 : 0);

  solve_gmp_R::solve(a,b);

  return(bigrationalR::create_SEXP(b));
}


SEXP inverse_z (SEXP A)
{
  bigvec a = bigintegerR::create_bignum(A);
  if(a.modulus.size() == 1 &&  !a.modulus[0].isNA()) {
    bigvec b (a.size() );
    b.nrow = a.nrow;
    if(a.nrow * a.nrow != a.size())
      error(_("Argument 1 must be a square matrix"));

    b.modulus.push_back(a.modulus[0]);
    // initialize b to identity
    for(unsigned int i=0; i<b.nrow ; ++i)
      for(unsigned int j=0; j<b.nrow ; ++j)
	b.value[i+j*b.nrow].setValue((i == j) ? 1 : 0);

    solve_gmp_R::solve(a,b);

    return(bigintegerR::create_SEXP(b));
  }
  else {
    bigvec_q aq (a);
    return(solve_gmp_R::inverse_q(aq));
  }
}

// solve AX=B
SEXP solve_z(SEXP A,SEXP B)
{
  bool common_modulus=true;
  bigvec a = bigintegerR::create_bignum(A);
  bigvec b = bigintegerR::create_bignum(B);
  if(a.modulus.size() == 1 )
    if(!a.modulus[0].isNA())
      {
	if(b.modulus.size()>1)
	  common_modulus = false; // -> solve with rational
	if(b.modulus.size() == 1)
	  {
	    if(b.modulus[0] == a.modulus[0])
	      common_modulus = false; // -> solve with rational
	  }
	else
	  b.modulus.push_back(a.modulus[0]);
	// solve in Z/nZ
	if(common_modulus)
	  {
	    // case: b a vector
	    if(b.nrow<1)
	      b.nrow = b.size();

	    if(a.nrow * a.nrow != a.size())
	      error(_("Argument 1 must be a square matrix"));

	    if(a.nrow != b.nrow)
	      error(_("Dimensions do not match"));

	    solve_gmp_R::solve(a,b);

	    return(bigintegerR::create_SEXP(b));
	  }
      }

  bigvec_q aq (a);
  bigvec_q bq (b);
  return(solve_gmp_R::solve_q(aq,bq));
}


// solve AX=B
SEXP solve_q(SEXP A,SEXP B)
{
  bigvec_q a = bigrationalR::create_bignum(A);
  bigvec_q b = bigrationalR::create_bignum(B);

  return(solve_gmp_R::solve_q(a,b));
}

// solve AX = B
SEXP solve_gmp_R::solve_q(bigvec_q a, bigvec_q b)
{
  if(a.nrow * a.nrow != a.size())
    error(_("Argument 1 must be a square matrix"));

  // case: b a vector
  if(b.nrow<0)
    b.nrow = b.size();

  if(a.nrow != b.nrow)
    error(_("Dimensions do not match"));

  solve_gmp_R::solve(a,b);

  return(bigrationalR::create_SEXP(b));
}
