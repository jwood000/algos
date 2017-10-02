/*! \file bigve_cq.h
 *  \brief bigvec_q class definition
 *
 *  \version 1
 *
 *  \date Created: 2005
 *  \date Last modified: Time-stamp: <2006-06-16 20:19:01 antoine>
 *
 *
 *  \note Licence: GPL
 */

#ifndef BIGVEC_Q_HEADER_
#define BIGVEC_Q_HEADER_ 1

#include "bigrational.h"
#include "bigvec.h"


/** \brief class bigvec_q for vector of bigrational
 *
 * It is a class composed of a vector of bigrational
 * and an nrow parameter.
 */
class bigvec_q {
 public:

  /** \brief The real value */
  std::vector<bigrational> value;
  /** \brief optional parameter used with matrix -- set to -1 for non-matrix */
  int nrow ;

  /** \brief constructor
   */
  bigvec_q() :
    value(),
    nrow(-1)
    {      }

  /** \brief  strange constructors ! */
  bigvec_q(double * & a, double * b):
    value(a,b),nrow(-1) {};
  /** \brief strange constructors ! */
  bigvec_q(int * & a, int * b):
    value(a,b),nrow(-1) {};

  /**  \brief create empty vector of size i */
  bigvec_q(unsigned int i):
    value(i),
    nrow(-1){};

  /**  \brief copy constructor */
  bigvec_q(const bigvec_q & rhs);

  /** \brief create from bigintegers */
  bigvec_q(const bigvec & rhs);

  /**  \brief assignemt operator */
  bigvec_q & operator= (const bigvec_q& rhs);

  /**
   * \brief extract a value.
   * \note it returns a copy of the value.
   */
  bigrational  operator[] (unsigned int i) const;

  /**
   * \brief Set a value
   */
  void set(unsigned int i,const bigrational & val);

  /**
   * \brief set a value
   * \note it change bigrational na flag to false
   */
  void set(unsigned int i,const mpq_t & val);


  /**
   * \brief add a value
   */
  void push_back(const bigrational &i);

  /**
   * \brief return size
   */
  unsigned int size() const ;

  /**
   * \brief resize value
   */
  void resize(unsigned int n);

  /**
   * \brief clear all.
   */
  void clear();

  /**
   * \brief substract lambda[0] * line j to line i
   */
  void subLine(unsigned int i,unsigned int j,bigvec_q lambda);

  /**
   * \brief multiply line i by lambda[0]
   */
   void mulLine(unsigned int i, bigvec_q lambda);

   /** \brief print matrix to stdout
    *
    * use for debug purpose
    */
   void print();

};

/**  \brief comparison operator
 */
bool operator!= (const bigvec_q& rhs, const bigvec_q& lhs);


#endif
