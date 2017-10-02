
#include <R.h>
#include <Rinternals.h>
// and one thing from Rdefines.h :
#define NEW_LIST(n) allocVector(VECSXP,n)

#include "apply.h"
#include "bigintegerR.h"
#include "bigrationalR.h"


// X a matrix, or a bigz integer
// line: true = we return a list of all lines
//
SEXP gmpMatToListZ(SEXP X, SEXP line)
{

  SEXP ans;
  // dangerous... no check => use with care in R
  int lines =  INTEGER(line)[0];

  bigvec matrix = bigintegerR::create_bignum(X);

  unsigned int ncol = matrix.size() / matrix.nrow;
  unsigned int nrow =  matrix.nrow;

  if(lines == 1)
    {
      // RETURN a list of all lines
      PROTECT (ans = NEW_LIST(matrix.nrow) );
      for(unsigned int i = 0; i < nrow; ++i)
	{
	  bigvec oneLine ;
	  for(unsigned int j = 0; j < ncol; ++j)
	    {
	      oneLine.value.push_back(matrix.value[i+j*nrow]);

	      // modulus, if one by cell
	      if(matrix.modulus.size() ==matrix.value.size() )
		oneLine.modulus.push_back(matrix.modulus[i+j*nrow]);

	    }

	  // modulus, if one by line
	  if(((matrix.modulus.size() == nrow ) || (matrix.modulus.size() == 1) ) && (matrix.modulus.size() !=matrix.value.size()) )
	    oneLine.modulus.push_back(matrix.modulus[i % matrix.modulus.size() ]);


	  SET_VECTOR_ELT(ans, i,bigintegerR::create_SEXP(oneLine));

	}
      UNPROTECT(1);
    }
  else
    {
      // RETURN a list of all rows !
      PROTECT (ans = NEW_LIST(ncol) );
      for(unsigned int j = 0; j < ncol; ++j)
	{
	  bigvec oneLine ;
	  for(unsigned int i = 0; i < nrow; ++i)
	    {
	      oneLine.value.push_back(matrix.value[i+j*nrow]);

	      // modulus, if one by cell
	      if(matrix.modulus.size() ==matrix.value.size() )
		oneLine.modulus.push_back(matrix.modulus[i+j*nrow]);

	    }

	  // modulus, if one by line
	  if( (matrix.modulus.size() == 1)  && (matrix.modulus.size() !=matrix.value.size()) )
	    oneLine.modulus.push_back(matrix.modulus[0 ]);


	  SET_VECTOR_ELT(ans, j,bigintegerR::create_SEXP(oneLine));

	}
      UNPROTECT(1);

    }

  return(ans);
}


// X a matrix, or a bigq rational
// line: true = we return a list of all lines
//
SEXP gmpMatToListQ(SEXP X, SEXP line)
{

  SEXP ans;
  // dangerous... no check => use with care in R
  bool lines =  INTEGER(line)[0];

  bigvec_q matrix = bigrationalR::create_bignum(X);

  unsigned int ncol = matrix.size() / matrix.nrow;
  unsigned int nrow =  matrix.nrow;

  if(lines)
    {
      // RETURN a list of all lines
      PROTECT (ans = NEW_LIST(matrix.nrow) );
      for(unsigned int i = 0; i < nrow; ++i)
	{
	  bigvec_q oneLine ;
	  for(unsigned int j = 0; j < ncol; ++j)
	    {
	      oneLine.value.push_back(matrix.value[i+j*nrow]);
	    }
	  SET_VECTOR_ELT(ans, i,bigrationalR::create_SEXP(oneLine));

	}
      UNPROTECT(1);
    }
  else
    {
      // RETURN a list of all rows !
      PROTECT (ans = NEW_LIST(ncol) );
      for(unsigned int j = 0; j < ncol; ++j)
	{
	  bigvec_q oneLine ;
	  for(unsigned int i = 0; i < nrow; ++i)
	    {
	      oneLine.value.push_back(matrix.value[i+j*nrow]);
	    }

	  SET_VECTOR_ELT(ans, j,bigrationalR::create_SEXP(oneLine));

	}
      UNPROTECT(1);

    }

  return(ans);
}


