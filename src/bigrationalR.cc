/************************************************************/
/*! \file bigrationalR.cc
 *  \brief C function to interface R and libgmp with big rational values
 *
 *  \version 1
 *
 *  \date Created: 12/12/04
 *  \date Last modified: $Id: bigrationalR.cc,v 1.22 2013-03-25 10:57:04 mmaechler Exp $
 *
 *  \author Antoine Lucas (adapted from biginteger class made by
 *                         Immanuel Scholz)
 *
 *  \note Licence: GPL
 */

#include "Rgmp.h"

#include "bigvec_q.h"

#include <vector>
#include <algorithm>
using namespace std;
#include "bigintegerR.h"
#include "bigrationalR.h"

#include "matrix.h"


namespace bigrationalR
{
  /** \brief create a vector of bigrationals, all without a denominator.
   */
  bigvec_q create_vector(SEXP param) {
    switch (TYPEOF(param)) {
    case NILSXP:
	return bigvec_q(); // = bigq(0)
    case RAWSXP:
      {
	// deserialise the vector. first int is the size.
	bigvec_q v;
	char* raw = (char*)RAW(param);
	int pos = sizeof(int); // position in raw[]. Starting after header.
	for (int i = 0; i < ((int*)raw)[0]; ++i) {
	  v.push_back(bigrational(bigrational((void*)&raw[pos])));
	  pos += v.value.back().raw_size(); // increment number of bytes read.
	}
	return v;
      }
    case REALSXP:
      {
	double* d = REAL(param);
	bigvec_q v(d,d+LENGTH(param));
	for (unsigned int j = 0; j < v.size(); ++j)
	  if ( R_FINITE ( d[j]) )
	    v.value[j].setValue(d[j]);
	  else
	    v.value[j].setValue();
	return v;
      }
    case INTSXP:
    case LGLSXP:
      {
	int* i = INTEGER(param);
	bigvec_q v(i,i+LENGTH(param));
	for (unsigned int j = 0; j < v.size(); ++j)
	  if (i[j] == NA_INTEGER)
	    v.value[j].setValue();
	return v;
      }
    case STRSXP:
      {
	bigvec_q v;
	v.value.reserve(LENGTH(param));
	for (int i = 0; i < LENGTH(param); ++i) {
	  if (STRING_ELT(param,i) == NA_STRING)
	    v.value.push_back(bigrational());
	  else
	    v.value.push_back(bigrational(std::string(CHAR(STRING_ELT(param,i)))));
	}
	return v;
      }
    default:
      {
	// no longer: can be fatal later! /* vector of size 0 */ return bigvec_q();
	error(_("only logical, numeric or character (atomic) vectors can be coerced to 'bigq'"));
      }
    }
  }

  bigvec_q create_bignum(SEXP param)
  {
    /*
    // First check if it is a bigz class
    // This should occur very rarely except if force
    // a bigq operator like add.bigq(z1,z2)...

    // this should be decommented when
    // biqg will not be exported in mpz format
    SEXP className;
    PROTECT(className = Rf_allocVector(STRSXP,1));
    SET_STRING_ELT(className, 0, Rf_mkChar("class"));
    SEXP classAttr = Rf_getAttrib(param, className);
    UNPROTECT(1);
    if (TYPEOF(classAttr) == STRSXP)
      if( CHAR(STRING_ELT(param,0)) == "bigz")
	return(bigvec_q(bigintegerR::create_bignum(param)) );
    */

    bigvec_q v = bigrationalR::create_vector(param);
    SEXP denAttr = Rf_getAttrib(param, Rf_mkString("denominator"));
    SEXP dimAttr = Rf_getAttrib(param, Rf_mkString("nrow"));
    if (TYPEOF(dimAttr) == INTSXP)
      v.nrow = INTEGER(dimAttr)[0];
    else {
	// catch to get std matrix dimensions value
	dimAttr = Rf_getAttrib(param, Rf_mkString("dim"));
	v.nrow = (TYPEOF(dimAttr) == INTSXP) ? INTEGER(dimAttr)[0] : -1;
      }
    if (TYPEOF(denAttr) != NILSXP)
      {
	bigvec_q attrib = bigrationalR::create_vector(denAttr);
	if (!attrib.value.empty()) // sanity check
	  for (unsigned int i = 0; i < v.size(); ++i)
	    {
	      if( attrib[i%attrib.size()].sgn() != 0)
		v.value[i].setDenValue (attrib.value[i%attrib.size()].getValueTemp()) ;
	    }
      }
    return v;
  }

  SEXP create_SEXP(const bigvec_q & v)
  {

    SEXP ans, R_denom;
    int sizenum = sizeof(int), // starting with vector-size-header
	sizedenum = sizenum;
    unsigned int i;
    mpz_t  num, den;

    mpz_init(num);
    mpz_init(den);
    mpz_t_sentry val_n(num);
    mpz_t_sentry val_d(den);

    int numb = 8*sizeof(int);

    //return sizeof(int) * (2 + (mpz_sizeinbase(value,2)+numb-1) / numb);

    for (i = 0; i < v.size(); ++i)
      {
	if(v.value[i].isNA())
	  {
	    sizenum += sizeof(int);
	    sizedenum += sizeof(int);
	  }
	else
	  {
	    // *num = mpq_numref(v.value[i].getValueTemp());
	    // *den = mpq_denref(v.value[i].getValueTemp());
	    mpq_get_num(num,v.value[i].getValueTemp());
	    mpq_get_den(den,v.value[i].getValueTemp());

	    sizenum += sizeof(int) * (2 + (mpz_sizeinbase(num,2)+numb-1) / numb); // adding each bigint's needed size
	    sizedenum += sizeof(int) * (2 + (mpz_sizeinbase(den,2)+numb-1) / numb);
	  }
      }
    PROTECT(ans = Rf_allocVector(RAWSXP, sizenum));
    PROTECT(R_denom = Rf_allocVector(RAWSXP, sizedenum));
    char* r = (char*)RAW(ans);
    char* rdenom = (char*)RAW(R_denom);
    ((int*)r)[0] =((int*)rdenom)[0] = v.size(); // first int is vector-size-header
    int posnum = sizeof(int); // current position in r[] (starting after vector-size-header)
    int posdenum = sizeof(int); // current position in r[] (starting after vector-size-header)
    for (i = 0; i < v.size(); ++i)
      {
	mpq_get_num(num,v.value[i].getValueTemp());
	mpq_get_den(den,v.value[i].getValueTemp());
	posnum += as_raw(&r[posnum],num,v.value[i].isNA());
	posdenum += as_raw(&rdenom[posdenum],den,v.value[i].isNA());
      }

    // set the class attribute to "bigrational"
    Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigq"));
    Rf_setAttrib(ans, Rf_mkString("denominator"), R_denom);

    // set the dim attribute to "bigq"
    if(v.nrow >= 0)
      Rf_setAttrib(ans, Rf_mkString("nrow"), Rf_ScalarInteger((int) v.nrow));

    UNPROTECT(2);
    return ans;
  }


  /**
   * \brief Main function of doing a binary operation on bigrationals.
   * It calls a function argument for doing the correct thing.
   * This could also be written as a class functor (template)
   * to save one function call, but then code bloat will happen.
   */
  SEXP bigrational_binary_operation(SEXP a, SEXP b, bigrational_binary_fn f)
  {
    bigvec_q va = bigrationalR::create_bignum(a);
    bigvec_q vb = bigrationalR::create_bignum(b), result;
    // MM: I think, 0 length vector operations should work too!
    // if (va.value.empty() || vb.value.empty())
    //   error(_("argument must not be an empty list"));
    int size = (va.value.empty() || vb.value.empty()) ? 0 : max(va.size(), vb.size());
    result.value.reserve(size);
    for (int i = 0; i < size; ++i)
      result.push_back(f(va.value[i%va.size()], vb.value[i%vb.size()]));
    result.nrow = matrixz::checkDims(va.nrow,vb.nrow) ;
    return bigrationalR::create_SEXP(result);
  }

  SEXP bigrational_bigz_binary_operation(SEXP a, SEXP b, bigrational_bigz_binary_fn f)
  {
    bigvec_q va = bigrationalR::create_bignum(a), result;
    bigvec vb = bigintegerR::create_bignum(b);
    int size = (va.value.empty() || vb.value.empty()) ? 0 : max(va.size(), vb.size());
    result.value.reserve(size);
    for (int i = 0; i < size; ++i)
      result.push_back(f(va.value[i%va.size()], vb.value[i%vb.size()]));
    result.nrow = matrixz::checkDims(va.nrow,vb.nrow) ;
    return bigrationalR::create_SEXP(result);
  }

  SEXP bigrational_logical_binary_operation(SEXP a, SEXP b, bigrational_logical_binary_fn f)
  {
    bigvec_q va = bigrationalR::create_bignum(a);
    bigvec_q vb = bigrationalR::create_bignum(b), result;
    // MM: I think, 0 length vector operations should work too!
    // if (va.value.empty() || vb.value.empty())
    //   error(_("argument must not be an empty list"));
    int size = (va.value.empty() || vb.value.empty()) ? 0 : max(va.size(), vb.size());
    SEXP ans = PROTECT(Rf_allocVector(LGLSXP, size));
    for (int i = 0; i < size; ++i) {
      bigrational am = va.value[i % va.size()];
      bigrational bm = vb.value[i % vb.size()];
      if (am.isNA() || bm.isNA())
	LOGICAL(ans)[i] = NA_LOGICAL;
      else
	LOGICAL(ans)[i] = f(va[i%va.size()], vb[i%vb.size()]) ? 1 : 0;
    }

    int nrow = matrixz::checkDims(va.nrow,vb.nrow) ;

    // Add dimension parameter when available
    if(nrow >= 0)
      {
	SEXP dimVal;
	PROTECT(dimVal = Rf_allocVector(INTSXP, 2));
	INTEGER(dimVal)[0] = (int) nrow;
	INTEGER(dimVal)[1] = (int) size / nrow;
	Rf_setAttrib(ans, Rf_mkString("dim"), dimVal);
	UNPROTECT(1);
      }

    UNPROTECT(1);
    return ans;
  }

  bool lt(const bigrational& lhs, const bigrational& rhs)
  {
    return mpq_cmp(lhs.getValueTemp(), rhs.getValueTemp()) < 0;
  }
  bool gt(const bigrational& lhs, const bigrational& rhs)
  {
    return mpq_cmp(lhs.getValueTemp(), rhs.getValueTemp()) > 0;
  }
  bool lte(const bigrational& lhs, const bigrational& rhs)
  {
    return mpq_cmp(lhs.getValueTemp(), rhs.getValueTemp()) <= 0;
  }
  bool gte(const bigrational& lhs, const bigrational& rhs)
  {
    return mpq_cmp(lhs.getValueTemp(), rhs.getValueTemp()) >= 0;
  }
  bool eq(const bigrational& lhs, const bigrational& rhs)
  {
    return mpq_cmp(lhs.getValueTemp(), rhs.getValueTemp()) == 0;
  }
  bool neq(const bigrational& lhs, const bigrational& rhs)
  {
    return mpq_cmp(lhs.getValueTemp(), rhs.getValueTemp()) != 0;
  }


  // Create a bigrational from a binary combination of two other bigrationals
  bigrational create_bigrational(const bigrational& lhs, const bigrational& rhs,
				 gmpq_binary f,  bool zeroRhsAllowed ) {
    if (lhs.isNA() || rhs.isNA())
      return bigrational();

    if (!zeroRhsAllowed && ( mpq_sgn(rhs.getValueTemp()) == 0))
      error(_("division by zero"));

    mpq_t val;
    mpq_init(val);
    mpq_t_sentry val_s(val);

    f(val,lhs.getValueTemp(),rhs.getValueTemp());

    /* Simplify numerator and denominator */
    mpq_canonicalize(val);

    return bigrational(val);
  }

  // Create a bigrational from a binary combination of  (bigrational , biginteger)
  bigrational create_bigrational_z(const bigrational& lhs, const biginteger& rhs,
				   gmp_qz_binary f, bool zeroRhsAllowed)
  {
    if (lhs.isNA() || rhs.isNA())
      return bigrational();
    if (!zeroRhsAllowed && ( mpz_sgn(rhs.getValueTemp()) == 0))
      error(_("division by zero"));

    mpq_t val; mpq_init(val); mpq_t_sentry val_s(val);

    f(val,lhs.getValueTemp(),rhs.getValueTemp());
    /* Simplify numerator and denominator and return */
    mpq_canonicalize(val);
    return bigrational(val);
  }

  //  x ^ y  (for biginteger y):
  void mpqz_pow(mpq_t result, const mpq_t x, const mpz_t y)
  {
    if(!mpz_fits_slong_p(y))
      error(_("exponent 'y' too large in 'x^y'"));

    mpz_t num, den; mpz_init(num); mpz_init(den);
    mpz_t_sentry val_n(num); mpz_t_sentry val_d(den);
    int yi = mpz_get_si(y);
    bool neg =(yi < 0);
    // *num = mpq_numref(x);
    // *den = mpq_denref(y);
    mpq_get_num(num, x);
    mpq_get_den(den, x);
    if(neg) {
	if(mpz_sgn(num) == 0)
	    error(_("0 ^ <negative> is a division by zero"));
	yi = -yi;
    }
    mpz_pow_ui(num, num, yi); // num := num ^ |y|
    mpz_pow_ui(den, den, yi); // den := den ^ |y|

    if(neg) { // Q^{-n} = 1 / Q^|n| --> result := den / num = as.bigq(den, num) :
	mpz_set(mpq_numref(result), den);
	mpz_set(mpq_denref(result), num);
    } else {
	// result := num / den = as.bigq(num, den) :
	mpz_set(mpq_numref(result), num);
	mpz_set(mpq_denref(result), den);
    }
    /* Simplify numerator and denominator */
    mpq_canonicalize(result);
  }


}
// End of namespace bigrationalR ----------------------------------------------



SEXP bigrational_add (SEXP a, SEXP b) {return bigrationalR::bigrational_binary_operation(a,b,operator+);}
SEXP bigrational_sub (SEXP a, SEXP b) {return bigrationalR::bigrational_binary_operation(a,b,operator-);}
SEXP bigrational_mul (SEXP a, SEXP b) {return bigrationalR::bigrational_binary_operation(a,b,operator*);}
SEXP bigrational_div (SEXP a, SEXP b) {return bigrationalR::bigrational_binary_operation(a,b,operator/);}
SEXP bigrational_pow (SEXP a, SEXP b) {
    return bigrationalR::bigrational_bigz_binary_operation(a,b,operator^); //-> mpqz_pow() above
}
SEXP bigrational_as (SEXP n, SEXP d) {return bigrationalR::bigrational_binary_operation(n,d,set_denominator);} //-> mpq_div()

SEXP bigrational_lt (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::lt);}
SEXP bigrational_gt (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::gt);}
SEXP bigrational_lte (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::lte);}
SEXP bigrational_gte (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::gte);}
SEXP bigrational_eq (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::eq);}
SEXP bigrational_neq (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::neq);}



SEXP bigrational_as_character(SEXP a, SEXP b)
{
  bigvec_q v = bigrationalR::create_bignum(a);

  int base = Rf_asInteger(b);
  SEXP ans = PROTECT(Rf_allocVector(STRSXP, v.size()));
  for (unsigned int i = 0; i < v.size(); ++i)
    SET_STRING_ELT(ans, i, Rf_mkChar(v.value[i].str(base).c_str()));

  // matrix part
  if(v.nrow >= 0)
    {
      SEXP nRow = PROTECT(Rf_allocVector(INTSXP, 2));
      INTEGER(nRow)[0] = v.nrow;
      INTEGER(nRow)[1] = v.value.size() / v.nrow;
      Rf_setAttrib(ans, Rf_mkString("dim"), nRow);
      UNPROTECT(1);
    }

  UNPROTECT(1);
  return ans;
}

SEXP bigrational_as_numeric(SEXP a)
{
  bigvec_q v = bigrationalR::create_bignum(a);
  SEXP ans = PROTECT(Rf_allocVector(REALSXP,v.size()));
  double *r = REAL(ans);
  for (unsigned int i = 0; i < v.size(); ++i)
    r[i] = v.value[i].isNA() ? NA_REAL : v.value[i].as_double();
  UNPROTECT(1);
  return ans;
}

SEXP bigrational_get_at(SEXP a, SEXP b)
{
  bigvec_q va = bigrationalR::create_bignum(a);
  vector<int> vb = bigintegerR::create_int(b);
  bigvec_q result;
  if (TYPEOF(b) == LGLSXP) {
    for (unsigned int i = 0; i < va.size(); ++i)
      if (vb[i%vb.size()])
	result.push_back(va.value[i]);
  } else {
    std::remove(vb.begin(), vb.end(), 0); // remove all zeroes
    if (vb.empty())
      return bigrationalR::create_SEXP(bigvec_q());
    if (vb[0] < 0) {
      for (vector<int>::iterator it = vb.begin(); it != vb.end(); ++it)
	if (*it > 0)
	  error(_("only 0's may mix with negative subscripts"));
	else if (-(*it)-1 >= (int)va.size())
	  error(_("subscript out of bounds"));
      // TODO: This is optimized for large va.size and small vb.size.
      // Maybe add a condition to use a different approach for large vb's
      result.value.reserve(va.size()-vb.size());
      for (int i = 0; i < (int)va.size(); ++i)
	if (find(vb.begin(), vb.end(), -i-1) == vb.end())
	  result.push_back(va.value[i]);
    } else {
      result.value.reserve(vb.size());
      for (vector<int>::iterator it = vb.begin(); it != vb.end(); ++it) {
	if (*it < 0)
	  error(_("only 0's may mix with negative subscripts"));
	if (*it <= (int)va.size())
	  result.push_back(va.value[(*it)-1]);
	else
	  result.push_back(bigrational()); // NA for out of range's
      }
    }
  }
  return bigrationalR::create_SEXP(result);
}

SEXP bigrational_set_at(SEXP src, SEXP idx, SEXP value)
{
  vector<int>::iterator it;
  int i;
  bigvec_q result = bigrationalR::create_bignum(src);
  bigvec_q vvalue = bigrationalR::create_bignum(value);
  vector<int> vidx = bigintegerR::create_int(idx);
  int pos = 0;

  if(vvalue.size() == 0) {
      if(result.size() == 0)
	  return bigrationalR::create_SEXP(result);
      else
	  error(_("replacement has length zero"));
  }
  if (TYPEOF(idx) == LGLSXP) {
    for (i = 0; i < (int)result.size(); ++i)
      if (vidx[i%vidx.size()])
	result.value[i] = vvalue.value[pos++%vvalue.size()];
  } else {
    std::remove(vidx.begin(), vidx.end(), 0); // remove all zeroes
    if (vidx.empty())
      return bigrationalR::create_SEXP(result);
    if (vidx[0] < 0) {
      for ( it = vidx.begin(); it != vidx.end(); ++it)
	if (*it > 0)
	  error(_("only 0's may mix with negative subscripts"));
	else if (-(*it)-1 >= (int)result.size())
	  error(_("subscript out of bounds"));
      pos = 0;
      for ( i = 0; i < (int)result.size(); ++i)
	if (find(vidx.begin(), vidx.end(), -i-1) == vidx.end())
	  result.value[i] = vvalue.value[pos++%vvalue.size()];
    } else {
      // finding maximum to resize vector if needed
      int maximum = INT_MIN;
      for (it = vidx.begin(); it != vidx.end(); ++it)
	maximum = max(maximum, *it);
      if (maximum > (int)result.size())
	result.value.resize(maximum);
      pos = 0;
      for (it = vidx.begin(); it != vidx.end(); ++it) {
	if (*it < 0)
	  error(_("only 0's may mix with negative subscripts"));
	result.value[(*it)-1] = vvalue[pos++%vvalue.size()];
      }
    }
  }
  return bigrationalR::create_SEXP(result);
}

SEXP bigrational_length(SEXP a)
{
  return Rf_ScalarInteger(bigrationalR::create_bignum(a).size());
}

SEXP bigrational_den(SEXP a)
{
  mpz_t z_tmp;
  mpz_init(z_tmp);
  bigvec_q v =bigrationalR::create_bignum(a);
  bigvec result;
  result.value.resize(v.size());

  for (unsigned int i = 0; i < v.size(); ++i) {
    mpq_get_den(z_tmp,v.value[i].getValueTemp());
    result.value[i].setValue(z_tmp);
  }
  mpz_clear(z_tmp);
  return bigintegerR::create_SEXP(result);
}

SEXP bigrational_num(SEXP a)
{
  mpz_t z_tmp;
  mpz_init(z_tmp);
  bigvec_q v =bigrationalR::create_bignum(a);
  bigvec result;
  result.resize(v.size());

  for (unsigned int i = 0; i < v.size(); ++i) {
    if(!v.value[i].isNA()) {
      mpq_get_num(z_tmp,v.value[i].getValueTemp());
      result.value[i].setValue(z_tmp);
    } // else: uninitialized, i.e., NA
  }
  mpz_clear(z_tmp);
  return bigintegerR::create_SEXP(result);
}

SEXP bigrational_setlength(SEXP vec, SEXP value)
{
  int len = 0;
  switch (TYPEOF(value)) {
  case INTSXP:
  case LGLSXP:
    if (LENGTH(value) != 1)
      error(_("invalid second argument"));
    len = *INTEGER(value);
    if (len < 0)
      error(_("vector size cannot be negative"));
    else if (len == NA_INTEGER)
      error(_("vector size cannot be NA"));
    break;
  case REALSXP:
    if (LENGTH(value) != 1)
      error(_("invalid second argument"));
    len = (int)*REAL(value);
    if (len < 0)
      error(_("vector size cannot be negative"));
    else if (! (R_FINITE (len ) ) )
      error(_("vector size cannot be NA, NaN, or Inf"));
    break;
  case STRSXP:
    // dunno why R spits out this strange error on "Length(foo) <- -1"
    // but I always follow the holy standard ;-)
    error(_("negative length vectors are not allowed"));
  default:
    error(_("invalid second argument"));
  }
  bigvec_q v =bigrationalR::create_bignum(vec);
  v.value.resize(len);
  return bigrationalR::create_SEXP(v);
}

SEXP bigrational_is_na(SEXP a)
{
  bigvec_q v = bigrationalR::create_bignum(a);
  SEXP ans = PROTECT(Rf_allocVector(LGLSXP, v.size()));
  int *a_ = LOGICAL(ans);
  for (unsigned int i = 0; i < v.size(); ++i)
    a_[i] = v.value[i].isNA();
  UNPROTECT(1);
  return ans;
}


SEXP bigrational_is_int(SEXP a)
{
  bigvec_q v = bigrationalR::create_bignum(a);
  SEXP ans = PROTECT(Rf_allocVector(LGLSXP, v.size()));
  int *a_ = LOGICAL(ans);
  mpz_t z_tmp;
  mpz_init(z_tmp);

  for (unsigned int i = 0; i < v.size(); ++i) {
    mpq_get_den(z_tmp,v.value[i].getValueTemp());
    a_[i] = mpz_cmp_ui (z_tmp, 1) == 0; // <==> numerator == 1
  }
  mpz_clear(z_tmp);
  UNPROTECT(1);
  return ans;
}

SEXP bigrational_c(SEXP args)
{
  //  if(TYPEOF( args ) != LISTSXP)
  //  error(_("should be a list"));
  bigvec_q result;
  for(int i = 0; i < Length(args); i++) {
      bigvec_q v = bigrationalR::create_bignum(VECTOR_ELT(args,i));
      for(unsigned int j=0; j < v.size(); j++)
	result.push_back(v.value[j]);
      v.value.clear();
    }

  return bigrationalR::create_SEXP(result);
}


SEXP bigrational_cbind(SEXP args)
{
  bigvec_q result = bigrationalR::create_bignum(VECTOR_ELT(args,0));
  if(result.nrow <= 0)
    result.nrow = result.size();

  for(int i= 1; i < Length(args); i++)
    {
      bigvec_q v = bigrationalR::create_bignum(VECTOR_ELT(args,i));
      for(int j=0; j < (int)v.size() ; j++)
	result.push_back(v[j]);
      v.clear();
    }

  return bigrationalR::create_SEXP(result);
}


SEXP bigrational_rep(SEXP x, SEXP times)
{
  bigvec_q v = bigrationalR::create_bignum(x), result;
  unsigned int i,j, rep = (unsigned int )INTEGER(AS_INTEGER(times))[0];

  result.value.reserve(v.size()*rep);
  for(i = 0 ; i< rep ; i++)
    for(j = 0 ; j < v.size() ; j++)
      result.push_back(v.value[j]);

  return bigrationalR::create_SEXP(result);
}



// Return max
SEXP bigrational_max(SEXP a, SEXP narm)
{
  bigvec_q va = bigrationalR::create_bignum(a), result;
  if(! va.size())
    return bigrationalR::create_SEXP(result);

  unsigned int maximum = 0;
  int na_remove = Rf_asInteger(narm);

  for(unsigned int i = 1 ; i < va.size(); ++i)
    {
      if(va.value[i].isNA() && ! na_remove)
	return(bigrationalR::create_SEXP(result));
      else
	if(!(va.value[i] <  va.value[maximum] ))
	  maximum = i; // if va.value[maximum = 0] is NA => false for the "<" => maximum changed = good
    }

  result.push_back(va.value[maximum]);

  return bigrationalR::create_SEXP(result);
}


// Return min

SEXP bigrational_min(SEXP a, SEXP narm)
{
  bigvec_q result, va = bigrationalR::create_bignum(a);
  if (! va.size())
    return bigrationalR::create_SEXP(result);

  unsigned int minimum = 0;
  int na_remove = Rf_asInteger(narm);

  for(unsigned int i = 1 ; i < va.size(); ++i)
    {
      if(va.value[i].isNA() && !na_remove)
	return(bigrationalR::create_SEXP(result));
      else
	if(!(va.value[i] >  va.value[minimum] ))
	  minimum = i; // if va.value[maximum = 0] is NA => false for the "<" => maximum changed = good
    }

  result.push_back(va.value[minimum]);

  return bigrationalR::create_SEXP(result);
}

// Return cumsum
SEXP bigrational_cumsum(SEXP a)
{
  bigvec_q result, va = bigrationalR::create_bignum(a);
  result.value.resize(va.value.size());

  mpq_t val;
  mpq_init(val);
  mpq_t_sentry val_s(val);

  for(unsigned int i = 0 ; i < va.size(); ++i) {
      if(va.value[i].isNA() ) {
	  break; // all last values are NA.
      }
      mpq_add(val,val,va.value[i].getValueTemp());

      result.value[i].setValue(val);
  }
  return(bigrationalR::create_SEXP(result));
}


// Return sum
SEXP bigrational_sum(SEXP a)
{
  bigvec_q result, va = bigrationalR::create_bignum(a);
  result.value.resize(1);

  mpq_t val;
  mpq_init(val);
  mpq_t_sentry val_s(val);

  for(unsigned int i = 0 ; i < va.size(); ++i) {
    if(va.value[i].isNA()) {
      break; // all last values are NA.
    }
    mpq_add(val,val,va.value[i].getValueTemp());
  }
  result.value[0].setValue(val);
  return(bigrationalR::create_SEXP(result));
}


// Return prod

SEXP bigrational_prod(SEXP a)
{
  bigvec_q result, va = bigrationalR::create_bignum(a);
  result.value.resize(1);

  mpq_t val;
  mpq_init(val);
  mpq_set_ui(val,1,1);
  mpq_t_sentry val_s(val);

  for(unsigned int i = 0 ; i < va.size(); ++i) {
    if(va.value[i].isNA() ) {
      return (bigrationalR::create_SEXP(result));
    }
    mpq_mul(val,val,va.value[i].getValueTemp());
  }
  result.value[0].setValue(val);
  return(bigrationalR::create_SEXP(result));
}

/* ================================ UNUSED ============================*/
// return x ^ y   x: "bigq"  y: "bigz" (or INTSXP, REALSXP -- TODO ?)
SEXP bigrational_R_pow(SEXP x, SEXP y)
{
  bigvec_q result, vx = bigrationalR::create_bignum(x);
  bigvec vy = bigintegerR::create_bignum(y);
  int size = (vx.value.empty() || vy.value.empty()) ? 0 : max(vx.size(), vy.size());

  mpq_t val; mpq_init(val); mpq_t_sentry val_s(val);
  mpz_t num, den; mpz_init(num); mpz_init(den);
  mpz_t_sentry val_n(num); mpz_t_sentry val_d(den);

  result.value.resize(size);

  for (int i = 0 ; i < size; i++) {
    int i_x = i % vx.value.size(),
        i_y = i % vy.value.size();

    //fails: val.NA(false);
    if(vx.value[i_x].isNA() ||
       vy.value[i_y].isNA())
	break;
    // else -- res[i] :=  x[i] ^ y[i]  =  (num ^ y_i) / (den ^ y_i) ----

    if(mpz_sgn(vy.value[i_y].getValueTemp()) < 0)
	error(_("Negative powers not yet implemented [i = %d]"), i_y +1);
    if (!mpz_fits_ulong_p(vy.value[i_y].getValueTemp()))
	error(_("exponent too large for pow  [i = %d]"), i_y +1);
    int y_i = mpz_get_ui(vy.value[i_y].getValueTemp());
    // *num = mpq_numref(vx.value[i].getValueTemp());
    // *den = mpq_denref(vx.value[i].getValueTemp());
    mpq_get_num(num, vx.value[i_x].getValueTemp());
    mpq_get_den(den, vx.value[i_x].getValueTemp());
    mpz_pow_ui(num, num, y_i); // num := num ^ y_i
    mpz_pow_ui(den, den, y_i); // den := den ^ y_i

    // val := as.bigq(num, den) :
    mpz_set(mpq_numref(val), num);
    mpz_set(mpq_denref(val), den);
    /* Simplify numerator and denominator */
    mpq_canonicalize(val);

    result.value[i].setValue(val);
  }

  return bigrationalR::create_SEXP(result);
} // ..._pow()

