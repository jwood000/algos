/*! \file matrixq.h
 *  \brief header for rational matrix functions set
 *
 *  \date Created: 2005
 *  \date Last modified: Time-stamp: <2006-05-26 14:35:51 antoine>
 *
 *
 *  \note Licence: GPL
 */

#include "bigvec_q.h"

extern "C"
{

  /**
   * \brief build a matrix x with dimensions p&q byrow: 0 or 1
   */
  SEXP as_matrixq(SEXP x, SEXP p, SEXP q, SEXP byrow, SEXP mod);

  /**
   * \brief transpose a "matrix" n x p into the transpose p x n
   */
  SEXP bigq_transposeR(SEXP x);

  /** \brief  matrix cross product */
  SEXP matrix_crossp_q (SEXP a, SEXP trans);
  /** \brief  matrix multiplication */
  SEXP matrix_mul_q (SEXP a, SEXP b, SEXP op);

  /** \brief for function rbind
   */
  SEXP bigrational_rbind(SEXP args) ;
}

/**
 * \brief Function designed for bigrational matrix
 */
namespace matrixq
{
  /**
   * \brief C++ function use to transpose matrix
   */
  bigvec_q bigq_transpose (const  bigvec_q & mat,int nr,int nc);


}
