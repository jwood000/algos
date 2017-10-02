/************************************************************/
/*! \file matrix.cc
 *  \brief C++ function to add matrix support
 *
 *  \version 1
 *
 *  \date Created: 19/02/06
 *  \date Last modified: Time-stamp: <2006-06-17 23:10:44 antoine>
 *
 *  \author A. Lucas
 *
 * \note
 *  as usually, matrix x[i,j] (n x p) is represented by a vector
 *              x[i + j x n]  (i=0..n-1 ; j=0..p-1)
 *
 *  \note Licence: GPL
 */

#include <vector>

using namespace std;

#include "Rgmp.h"

#include "biginteger.h"
#include "bigintegerR.h"
#include "matrix.h"
// need to call matrix_mul_q()
#include "matrixq.h"

// given that x is "bigz" or "bigq",
// return TRUE if x is a bigz/q *matrix*: R's is.matrixZQ(.)
SEXP is_matrix_zq(SEXP x) {
    return Rf_ScalarLogical(Rf_getAttrib(x, Rf_mkString("nrow")) != R_NilValue);
}

// C++ side of R function matrix.bigz()
SEXP as_matrixz (SEXP x, SEXP nrR, SEXP ncR, SEXP byrowR, SEXP mod)
{
  int
    nc=INTEGER(ncR)[0],
    nr=INTEGER(nrR)[0],
    byrow=INTEGER(byrowR)[0];

  // get "bigz" vector, this makes all conversion int to bigz etc...
  bigvec mat = bigintegerR::create_bignum(x);
  int lendat = mat.value.size();
  // int sizemod = mat.modulus.size();
  // when modulus specified
  bigvec modulus = bigintegerR::create_bignum(mod);
  if(modulus.value.size()>0) // should be allways the case
    if(!modulus.value[0].isNA()) {
      mat.modulus.resize(modulus.size());
      for (unsigned int i = 0; i < modulus.size(); i++)
	mat.modulus[i] = modulus.value[i];
      // sizemod = modulus.size();
    }

  // A copy from R project to get correct dimension
  // all warnings...
  //
  if (nr == NA_INTEGER) // This is < 0
    error(_("matrix: invalid 'nrow' value (too large or NA)"));
  if (nr < 0)
    error(_("matrix: invalid 'nrow' value (< 0)"));
  if (nc < 0)
    error(_("matrix: invalid 'ncol' value (< 0)"));
  if (nc == NA_INTEGER)
    error(_("matrix: invalid 'ncol' value (too large or NA)"));

  if(lendat > 0 ) {
    if (lendat > 1 && (nr * nc) % lendat != 0) {
      if (((lendat > nr) && (lendat / nr) * nr != lendat) ||
	  ((lendat < nr) && (nr / lendat) * lendat != nr))
	warning(_("data length [%d] is not a sub-multiple or multiple of the number of rows [%d] in matrix"),
		lendat, nr);
      else if (((lendat > nc) && (lendat / nc) * nc != lendat) ||
	       ((lendat < nc) && (nc / lendat) * lendat != nc))
	warning(_("data length [%d] is not a sub-multiple or multiple of the number of columns [%d] in matrix"),
		lendat, nc);
    }
    else if ((lendat > 1) && (nr * nc == 0)){
      warning(_("data length exceeds size of matrix"));
    }
  }

  // update dimension parameters
  if(nr == 1)
    nr = (int)ceil(lendat / (double) nc);
  if(nc == 1)
    nc = (int)ceil(lendat / (double)nr);

  // when we extend "x"
  if(nc*nr > lendat)
    {
      mat.value.resize(nr*nc);
      for(int i = lendat; i < nr*nc; i++)
	mat.value[i] = mat.value[i % lendat];
    }
  mat.nrow = nr;
  if(byrow)
    {
      bigvec mat2 = matrixz::bigint_transpose (mat, nc,nr);
      mat2.nrow = nr;// FIXME - needed ??
      return( bigintegerR::create_SEXP (mat2));
    }

  return( bigintegerR::create_SEXP (mat));
}


/*
 * Transposition
 */
SEXP bigint_transposeR(SEXP x)
{
  SEXP dimAttr = Rf_getAttrib(x, Rf_mkString("nrow"));
  bigvec mat = bigintegerR::create_bignum(x);
  int nr, n = mat.size();

  if (dimAttr == R_NilValue) { // vector
    nr = n;
  } else if (TYPEOF(dimAttr) == INTSXP) {
    nr = INTEGER(dimAttr)[0];
  } else { nr = -1;// -Wall
    error(_("argument must be a matrix of class \"bigz\""));
  }
  int nc = (int) n / nr;
  // Rprintf(" o bigI_tr(<%d x %d>) ..\n", nr,nc);
  return( bigintegerR::create_SEXP(matrixz::bigint_transpose(mat, nr,nc)));
}


/* \brief  matrix cross product
 *
 * \param a matrix (n x p)
 * \param trans if(trans), compute tcrossprod(), else crossprod()
 * \return  crossprod(a) := t(a) %*% a  [p x p]  or
 *         tcrossprod(a) := a %*% t(a)  [n x n]
 */
SEXP matrix_crossp_z (SEXP a, SEXP trans)
{
  bool useMod = FALSE,
    tr = (bool)Rf_asLogical(trans);
  bigvec mat_a = bigintegerR::create_bignum(a);
  int sizemod = mat_a.modulus.size(),
    a_nrow = mat_a.nrow,
    a_len = mat_a.size();

  // in case of a vector; crossprod() returns scalar product,
  // whereas             tcrossprod() gives  n x n matrix.
  if(a_nrow < 0)
    a_nrow = a_len;
  int a_ncol = a_len / a_nrow;

  // Result R is R[1..m, 1..m] -- and R_{ij} = sum_{k=1}^p  A.. B..
  int m, p;
  if(tr) { // tcrossprod()
    m= a_nrow; p= a_ncol;
  } else { //  crossprod()
    m= a_ncol; p= a_nrow;
  }
  bigvec res(m*m);
  res.nrow= m;

  mpz_t R_ij, tt;
  mpz_init(R_ij);
  mpz_init(tt);
  mpz_t common_modulus; mpz_init(common_modulus);

  if(sizemod <= 1) { // maybe 'useMod' i.e., can use common modulus:
    if(sizemod == 1) {
      mpz_set(common_modulus, mat_a.modulus[0].getValueTemp());
      if(!mat_a.modulus[0].isNA())
	useMod = TRUE;
    }
  }

  // here the computation:
  for(int i=0; i < m; i++)
    for(int j=0; j < m; j++) {
      mpz_set_ui(R_ij, 0);
      bool isna = false;
#define K_LOOP								\
      for(int k=0; k < p; k++) {					\
	/* R_ij = \sum_{k=1}^p  a_{ik} b_{kj} */			\
	if( !(A_I_K.isNA() || B_K_J.isNA())) {				\
	  mpz_mul(tt, A_I_K.getValueTemp(), B_K_J.getValueTemp());	\
	  mpz_add(R_ij, tt,R_ij);					\
	}								\
	else {								\
	  isna = true; break;						\
	}								\
      }

      if(tr) {//------------- tcrossprod ---------------------------

#define A_I_K  mat_a.value [ i + k *a_nrow]
#define B_K_J  mat_a.value [ j + k *a_nrow]
	K_LOOP
#undef A_I_K
#undef B_K_J

       } else {//------------- crossprod ---------------------------

#define A_I_K  mat_a.value [ k + i *a_nrow]
#define B_K_J  mat_a.value [ k + j *a_nrow]
	K_LOOP
#undef A_I_K
#undef B_K_J

       }

      if(isna) {
	res.value[i + j*m].setValue(0);
	res.value[i + j*m].NA(true);
      }
      else
	res.value[i + j*m].setValue(R_ij);

    } // for(i ..)  for(j ..)

  if(useMod)
    res.modulus.push_back(biginteger(common_modulus));

  mpz_clear(R_ij);
  mpz_clear(tt);
  mpz_clear(common_modulus);

  return( bigintegerR::create_SEXP (res));
} // matrix_crossp_z()
#undef K_LOOP

/* \brief  matrix multiplication
 *
 * returns matrix multiplication  T(a) %*% b  or  b %*% T(a)
 * \param a matrix
 * \param b matrix
 * \param op operation code: 0: %*%,  1: crossprod,  2: tcrossprod
 *     (same codes as in R's do_matprod() in src/main/array.c )
 */
SEXP matrix_mul_z (SEXP a, SEXP b, SEXP op)
{
  if(!strcmp(class_P(b), "bigq")) { // b  "bigq",  --> use q arithm:
      return(matrix_mul_q(a, b, op));
  }
  // FIXME: we may know that a is 'bigz' - but we don't know at all about b !!
  // -----  create_bignum(.) should be much more careful (better: have a careful option!)
  bool useMod = FALSE;// if(useMod)  use a *common* modulus
  int o_ = Rf_asInteger(op); // INTEGER(op)[0]
  bigvec mat_a = bigintegerR::create_bignum(a),
         mat_b = bigintegerR::create_bignum(b);

  int sizemod_a = mat_a.modulus.size(),
      sizemod_b = mat_b.modulus.size();

  int
    a_nrow = mat_a.nrow, a_len = mat_a.size(),
    b_nrow = mat_b.nrow, b_len = mat_b.size(),
    a_ncol = -1, b_ncol = -1;// -Wall

  // distinguish cases of vectors / matrices ---------------------
  if(a_nrow < 0) {
    if(b_nrow < 0) { // *both* are vectors
      if(o_ == 0) {
	a_nrow = 1;
	a_ncol = a_len;
      } else {
	a_nrow = a_len;
	a_ncol = 1;
      }
      b_nrow = b_len;
      b_ncol = 1;

    } else { // a : vector,   b : matrix
      b_ncol = b_len / b_nrow;
      if(o_ == 0) {
	if (a_len == b_nrow) {	/* x as row vector */
	  a_nrow = 1;
	  a_ncol = b_nrow; /* == a_len */
	}
	else if (b_nrow == 1) {	/* x as col vector */
	  a_nrow = a_len;
	  a_ncol = 1;
	}
      } else if(o_ == 1) { /* crossprod() */
	if (a_len == b_nrow) {	/* x is a col vector */
	  a_nrow = b_nrow; /* == a_len */
	  a_ncol = 1;
	}
	/* else if (b_nrow == 1) ... not being too tolerant
	   to treat x as row vector, as t(x) *is* row vector */
      } else { // o_ == 2 -- tcrossprod()
	if (a_len == b_ncol) {	/* x as row vector */
	  a_nrow = 1;
	  a_ncol = b_ncol; /* == a_len */
	}
	else if (b_ncol == 1) {	/* x as col vector */
	  a_nrow = a_len;
	  a_ncol = 1;
	}
      }
    }
  }
  else if (b_nrow < 0) { // a : matrix,   b : vector
    a_ncol = a_len / a_nrow;
    if (o_ == 0) {
      if (b_len == a_ncol) {	/* y as col vector */
	b_nrow = a_ncol;
	b_ncol = 1;
      }
      else if (a_ncol == 1) {	/* y as row vector */
	b_nrow = 1;
	b_ncol = b_len;
      }
    }
    else if (o_ == 1) { /* crossprod() */
      if (b_len == a_nrow) {	/* y is a col vector */
	b_nrow = a_nrow;
	b_ncol = 1;
      }
    }
    else { /* tcrossprod --	   y is a col vector */
      b_nrow = b_len;
      b_ncol = 1;
    }

  } else { // a, b  *both* matrices
    a_ncol = a_len / a_nrow;
    b_ncol = b_len / b_nrow;
  }

  if(((o_ == 0) && (a_ncol != b_nrow)) ||
     ((o_ == 1) && (a_nrow != b_nrow)) || // crossprod()
     ((o_ == 2) && (a_ncol != b_ncol))    // tcrossprod()
     )
    error(_("Matrix dimensions do not match"));

  // Result R is R[1..n, 1..m] -- and R_{ij} = sum_{k=1} ^ p  A.. B..
  int n,m, p;
  if(o_ == 0) {
    n= a_nrow; m= b_ncol; p= a_ncol;// = b_nrow
  }else if (o_ == 1) {
    n= a_ncol; m= b_ncol; p= a_nrow;// = b_nrow
  }else if (o_ == 2) {
    n= a_nrow; m= b_nrow; p= a_ncol;// = b_ncol
  }else {
    error(_("invalid 'op' code in matrix_mul_z()"));
    n = m = p = -1;// -Wall
  }

  bigvec res(n*m);
  res.nrow=n;

  mpz_t common_modulus, tt;
  mpz_init(tt);
  mpz_init(common_modulus);

  /* modulus when modulus are "global" (i.e. of size 1) and
   * either are the same, or only one of a or b is specified
   */
  if( !(sizemod_a > 1 || sizemod_b > 1)) {
    if((sizemod_a == 1) && (sizemod_b == 1)) {
      mpz_set(common_modulus, mat_a.modulus[0].getValueTemp());
      if(mpz_cmp(common_modulus, mat_b.modulus[0].getValueTemp()) == 0
	 && !mat_a.modulus[0].isNA())
	useMod = TRUE;
    }
    else { // at least one of the sizemod_*  is > 1 :
      if ((sizemod_a == 1) && !mat_a.modulus[0].isNA()) {
	mpz_set(common_modulus, mat_a[0].modulus.getValueTemp());
	useMod = TRUE;
      } else if ((sizemod_b == 1) && !mat_b.modulus[0].isNA()) {
	mpz_set(common_modulus, mat_b[0].modulus.getValueTemp());
	useMod = TRUE;
      }
    }
  }
  // bigmod tmp;

  // here the computation:
  for(int i=0; i < n; i++)
    for(int j=0; j < m; j++)
      {
#define	R_IJ res.value[ i + j*n]
#define K_LOOP								\
	for(int k=0; k < p; k++)					\
	    {								\
	      if(A_I_K.isNA() || B_K_J.isNA()) {			\
		R_IJ.setValue(0); R_IJ.NA(true);			\
		break;							\
	      }								\
	      /* Z = A_I_K * B_K_J */					\
	      mpz_mul(tt, A_I_K.getValueTemp(), B_K_J.getValueTemp());	\
	      /* R_IJ = R_IJ + A_I_K * B_K_J  */			\
	      mpz_add(tt, tt, R_IJ.getValueTemp());			\
	      if(useMod)						\
		mpz_mod(tt,tt,common_modulus);				\
	      R_IJ.setValue(tt);					\
	    }

	R_IJ.setValue(0);

	if(o_ == 0) { //------------- %*% ---------------------------

#define A_I_K  mat_a.value [ i + k *a_nrow]
#define B_K_J  mat_b.value [ k + j *b_nrow]
	  K_LOOP
#undef A_I_K
#undef B_K_J

	} else if(o_ == 1){//------------- crossprod ---------------------------

#define A_I_K  mat_a.value [ k + i *a_nrow]
#define B_K_J  mat_b.value [ k + j *b_nrow]
	  K_LOOP
#undef A_I_K
#undef B_K_J

	} else {//(o_ == 2) ------------- tcrossprod ---------------------------

#define A_I_K  mat_a.value [ i + k *a_nrow]
#define B_K_J  mat_b.value [ j + k *b_nrow]
	  K_LOOP
#undef A_I_K
#undef B_K_J

	}
      }
  if(useMod)
    res.modulus.push_back(biginteger(common_modulus));

  mpz_clear(tt);
  mpz_clear(common_modulus);

  return( bigintegerR::create_SEXP (res));
} // matrix_mul_z()
#undef R_IJ
#undef K_LOOP


SEXP biginteger_rbind(SEXP args)
{
  int i=0,j=0;
  bigvec result;
  bigvec v;

  result = bigintegerR::create_bignum(VECTOR_ELT(args,0));
  if(result.nrow==0)
    result.nrow = result.size();

  result = matrixz::bigint_transpose(result, result.nrow,
				     result.size() / result.nrow);
  for(i=1; i < LENGTH(args); i++)
    {
      v = bigintegerR::create_bignum(VECTOR_ELT(args,i));
      if(v.nrow == 0 )
	v.nrow = v.size();
      v = matrixz::bigint_transpose(v,v.nrow,v.size() / v.nrow);

      for(j=0; j< (int)v.size(); j++)
	result.push_back(v[j]);
      v.clear();
    }

  result = matrixz::bigint_transpose(result, result.nrow,
				     result.size() / result.nrow);
  return bigintegerR::create_SEXP(result);
}



namespace matrixz
{
  bigvec bigint_transpose ( bigvec & mat,int nr,int nc)
  {
    int i,j;

    /* cas: square matrix */
    bigvec matbis (nr * nc);
    matbis.nrow = nc;

    /* we compute transpose */
    for(i =0; i<nr; i++)
      for(j =0; j<nc; j++)
	matbis.set(j+i*nc,mat[i+j*nr]);

    return matbis;
  }

  /* return dimension in dima */
  int checkDims(int dima, int dimb)
  {
    if(dima > 0 && dimb > 0) {
      if (dimb != dima)
	error(_("Matrix dimensions do not match"));
    }
    else { /* either a or b is a matrix */
	if(dima == -1)
	  return(dimb);
    }
    return(dima);
  }
}
