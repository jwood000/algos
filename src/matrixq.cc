/************************************************************/
/*! \file matrixq.cc
 *  \brief C++ function to add matrix support
 *
 *  \version 1
 *
 *  \date Created: 19/02/06
 *  \date Last modified: Time-stamp: <2008-02-17 19:54:11 antoine>
 *
 *  \author A. Lucas
 *
 * \note
 *  as usually, matrix x[i,j] (n x p) is represented by a vector
 *              x[i + j x n]  (i=0..n-1 ; j=0..p-1)
 *
 *  \note Licence: GPL
 */

//#include <list>
#include <vector>

using namespace std;

#include "Rgmp.h"

#include "bigrational.h"
#include "bigrationalR.h"
//#include "bigintegerR.h"
#include "matrixq.h"


// function called by matrix.bigz()
SEXP as_matrixq (SEXP x, SEXP nrR, SEXP ncR, SEXP byrowR, SEXP den)
{
  /* get "bigz" vector, this makes all conversion int to bigz etc...*/
  bigvec_q mat = bigrationalR::create_bignum(x),
    denominator = bigrationalR::create_bignum(den);
  int
    nc= INTEGER(ncR)[0],
    nr= INTEGER(nrR)[0],
    byrow= INTEGER(byrowR)[0],
    lendat = mat.size();

  if(denominator.value.size()>0) // should be allways the case
      if(!denominator.value[0].isNA())
	{
	  for (unsigned int i = 0; i < mat.size(); i++)
	    if(mat.value[i].isNA() && (denominator.value[i%denominator.size()].sgn() != 0) )
	      mat.value[i].setDenValue (denominator.value[i%denominator.size()].getValueTemp());

	}

  /* A copy from R project to get correct dimension
   * all warnings...
   */
  if (nr == NA_INTEGER) /* This is < 0 */
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
	warning(_("data length [%d] is not a sub-multiple or multiple of the number of rows [%d] in matrix"), lendat, nr);
      else if (((lendat > nc) && (lendat / nc) * nc != lendat) ||
	       ((lendat < nc) && (nc / lendat) * lendat != nc))
	warning(_("data length [%d] is not a sub-multiple or multiple of the number of columns [%d] in matrix"), lendat, nc);
    }
    else if ((lendat > 1) && (nr * nc == 0)){
      warning(_("data length exceeds size of matrix"));
    }
  }

  /* update dimension parameters */
  if(nr == 1)
    nr = (int)ceil(lendat / (double) nc);
  if(nc == 1)
    nc = (int)ceil(lendat / (double)nr);

  /* when we extend "x"  */
  if(nc*nr > lendat)
    {
      mat.value.resize(nr*nc);
      for(int i = lendat; i < nr*nc; i++)
	mat.value[i] = mat.value[i % lendat];
    }
  mat.nrow = nr;
  if(byrow)
   {
      bigvec_q mat2 = matrixq::bigq_transpose (mat, nc,nr);
      mat2.nrow = nr;
      return( bigrationalR::create_SEXP (mat2));
   }

  return( bigrationalR::create_SEXP (mat));
}


// function called by t(m) when m is a bigrational
SEXP bigq_transposeR(SEXP x)
{
  SEXP dimAttr = Rf_getAttrib(x, Rf_mkString("nrow"));
  bigvec_q mat = bigrationalR::create_bignum(x);
  int nr, n = mat.size();

  if (dimAttr == R_NilValue) { // vector
    nr = n;
  } else if (TYPEOF(dimAttr) == INTSXP) {
    nr = INTEGER(dimAttr)[0];
  } else {
    error(_("argument must be a matrix of class \"bigq\""));
    nr = -1;// -Wall
  }

  int nc = (int) n / nr;
  bigvec_q mat_transp = matrixq::bigq_transpose(mat, nr,nc);
  mat_transp.nrow = nc; // FIXME - needed ?
  return( bigrationalR::create_SEXP( mat_transp));
}



/* \brief  matrix cross product
 *
 * returns  crossprod(a) := t(a) %*% a  [p x p]  or
 *         tcrossprod(a) := a %*% t(a)  [n x n]
 * \param a matrix (n x p)
 * \param trans if(trans), compute tcrossprod(), else crossprod()
 */
SEXP matrix_crossp_q (SEXP a, SEXP trans)
{
  bool tr = (bool)Rf_asLogical(trans);
  bigvec_q mat_a = bigrationalR::create_bignum(a);
  int
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
  bigvec_q res(m*m);
  res.nrow= m;

  mpq_t R_ij, tt;
  mpq_init(R_ij);
  mpq_init(tt);

  // here the computation:
  for(int i=0; i < m; i++)
    for(int j=0; j < m; j++) {
	mpq_set_ui(R_ij, 0,1);
	bool isna = false;
#define K_LOOP								\
	for(int k=0; k < p; k++) {					\
	  /* R_ij = \sum_{k=1}^p  a_{ik} b_{kj} */			\
	  if( !(A_I_K.isNA() || B_K_J.isNA())) {			\
	      mpq_mul(tt, A_I_K.getValueTemp(), B_K_J.getValueTemp());	\
	      mpq_add(R_ij, tt,R_ij);					\
	  }								\
	  else {							\
	    isna = true; break;						\
	  }								\
	}

	if(tr) {//------------- tcrossprod ---------------------------

#define A_I_K  mat_a[ i + k *a_nrow]
#define B_K_J  mat_a[ j + k *a_nrow]
	  K_LOOP
#undef A_I_K
#undef B_K_J

	} else {//------------- crossprod ---------------------------

#define A_I_K  mat_a[ k + i *a_nrow]
#define B_K_J  mat_a[ k + j *a_nrow]
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
        }

  mpq_clear(R_ij);
  mpq_clear(tt);

  return( bigrationalR::create_SEXP (res));
} // matrix_crossp_q()


/** \brief matrix multiplication
 *
 * returns matrix multiplication  T(a) %*% b  or  b %*% T(a)
 * \param a is of dimension n x p
 * \param b is of dimension p x m
 * \param op operation code: 0: %*%,  1: crossprod,  2: tcrossprod
 *     (same codes as in R's do_matprod() in src/main/array.c )
 */
SEXP matrix_mul_q (SEXP a, SEXP b, SEXP op)
{
  int o_ = Rf_asInteger(op); // INTEGER(op)[0]
  bigvec_q mat_a = bigrationalR::create_bignum(a),
           mat_b = bigrationalR::create_bignum(b);

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


  bigvec_q res(n*m);
  res.nrow=n;

  mpq_t tt, R_ij;
  mpq_init(R_ij);
  mpq_init(tt);

  // here the computation:
  for(int i=0; i < n; i++)
    for(int j=0; j < m; j++) {
	mpq_set_ui(R_ij, 0,1);
	bool isna = false;
#define K_LOOP								\
	for(int k=0; k < p; k++) {					\
	  /* R_ij = \sum_{k=1}^p  a_{ik} b_{kj} */			\
	  if( !(A_I_K.isNA() || B_K_J.isNA())) {			\
	      mpq_mul(tt, A_I_K.getValueTemp(), B_K_J.getValueTemp());	\
	      mpq_add(R_ij, tt,R_ij);					\
	  }								\
	  else {							\
	    isna = true; break;						\
	  }								\
	}

	if(o_ == 0) { //------------- %*% --------------------------------------

#define A_I_K  mat_a[ i + k *a_nrow]
#define B_K_J  mat_b[ k + j *b_nrow]
	  K_LOOP
#undef A_I_K
#undef B_K_J

	} else if(o_ == 1){//------------- crossprod ---------------------------

#define A_I_K  mat_a[ k + i *a_nrow]
#define B_K_J  mat_b[ k + j *b_nrow]
	  K_LOOP
#undef A_I_K
#undef B_K_J

	} else {//(o_ == 2) ------------- tcrossprod ---------------------------

#define A_I_K  mat_a[ i + k *a_nrow]
#define B_K_J  mat_b[ j + k *b_nrow]
	  K_LOOP
#undef A_I_K
#undef B_K_J

	}

	if(isna) {
	  res.value[i + j*n].setValue(0);
	  res.value[i + j*n].NA(true);
	}
	else
	  res.value[i + j*n].setValue(R_ij);
        }

  mpq_clear(R_ij);
  mpq_clear(tt);

  return( bigrationalR::create_SEXP (res));
} // matrix_mul_q()
#undef K_LOOP


SEXP bigrational_rbind(SEXP args)
{
  int i=0,j=0;
  bigvec_q result;
  bigvec_q v;

  result = bigrationalR::create_bignum(VECTOR_ELT(args,0));
  if(result.nrow ==0)
    result.nrow = result.size();

  result = matrixq::bigq_transpose(result, result.nrow,
				   result.size() / result.nrow);
  for(i=1; i< LENGTH(args); i++) {
    v = bigrationalR::create_bignum(VECTOR_ELT(args,i));
    if(v.nrow == 0 )
      v.nrow = v.size();
    v = matrixq::bigq_transpose(v,v.nrow,v.size() / v.nrow);

    for(j=0; j< (int)v.size(); j++)
      result.push_back(v[j]);
    v.clear();
  }

  result = matrixq::bigq_transpose(result, result.nrow,
				   result.size() / result.nrow);
  return bigrationalR::create_SEXP(result);
}


bigvec_q matrixq::bigq_transpose (const  bigvec_q & mat,int nr,int nc)
{
  int i,j;
  bigvec_q matbis (nr * nc);
  matbis.nrow = nc;
  /* we compute transpose */
  for(i=0; i < nr; i++)
    for(j=0; j < nc; j++)
      matbis.value[j+i*nc].setValue(mat.value[i+j*nr]);

  return(matbis);
}
