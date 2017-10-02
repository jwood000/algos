#include "extract_matrix.h"

#include "bigrationalR.h"

// for something like x = A[indi, indj], but also simply  A[ind]
SEXP matrix_get_at_q(SEXP A,SEXP INDI, SEXP INDJ)
{
  bigvec_q mat = bigrationalR::create_bignum(A);

  return(bigrationalR::create_SEXP(extract_gmp_R::get_at( mat,INDI,INDJ)));
}

// for something like x = A[indi, indj], but also simply  A[ind]
SEXP matrix_get_at_z(SEXP A,SEXP INDI, SEXP INDJ)
{
  bigvec mat = bigintegerR::create_bignum(A);
  bigvec mat2 = extract_gmp_R::get_at( mat,INDI,INDJ);

  // now modulus !
  // cell based modulus
  if(mat.modulus.size() == mat.value.size())
    {
      for(unsigned int i = 0; i< mat.size(); ++i)
	mat.value[i] = mat.modulus[i];

      mat = extract_gmp_R::get_at( mat,INDI,INDJ);

      mat2.modulus.resize(mat.size());
      for(unsigned int i = 0; i< mat.size(); ++i)
	mat2.modulus[i] = mat.value[i];
    }
  // row base modulus
  else if((int)mat.modulus.size() == mat.nrow)
    {
      for(unsigned int i = 0; i< mat.size(); ++i)
	mat.value[i] = mat.modulus[i];

      mat.modulus.clear();

      mat = bigintegerR::biginteger_get_at_C(mat,INDI);

      mat2.modulus.resize(mat.size());
      for(unsigned int i = 0; i< mat.size(); ++i)
	mat2.modulus[i] = mat.value[i];

    }
  //global modulus
  else if(mat.modulus.size() == 1)
    {
      mat2.modulus.resize(1);
      mat2.modulus[0] = mat.modulus[0];
    }

  return(bigintegerR::create_SEXP(mat2) );
}



// for something like A[indi, indj] <- val
SEXP matrix_set_at_z(SEXP A, SEXP VAL, SEXP INDI, SEXP INDJ)
{
  bigvec mat = bigintegerR::create_bignum(A);

  if(TYPEOF(INDI) != LGLSXP ) {
      if(!length(INDI)) return(A);
      std::vector<int> vidx = bigintegerR::create_int(INDI);
      for(std::vector<int>::const_iterator it = vidx.begin();
	  it != vidx.end();
	  ++it)
	if(*it >= static_cast<int>(mat.size())) // in this case: we extend the vector
	  return( biginteger_set_at(A,INDI,VAL) );
    }
  bigvec val = bigintegerR::create_bignum(VAL);
  extract_gmp_R::set_at( mat,val,INDI,INDJ);
  return(bigintegerR::create_SEXP(mat));

}

// for something like A[indi, indj] <- val
SEXP matrix_set_at_q(SEXP A,SEXP VAL ,SEXP INDI, SEXP INDJ)
{
  bigvec_q mat = bigrationalR::create_bignum(A);

  if(TYPEOF(INDI) != LGLSXP ) {
      if(!length(INDI)) return(A);
      std::vector<int> vidx = bigintegerR::create_int(INDI);
      for(std::vector<int>::const_iterator it = vidx.begin();
	  it != vidx.end();
	  ++it)
	if(*it >= static_cast<int>(mat.size())) // in this case: we extend the vector
	  return( bigrational_set_at(A,INDI,VAL) );
    }

  bigvec_q val = bigrationalR::create_bignum(VAL);

  extract_gmp_R::set_at( mat,val,INDI,INDJ);

  return(bigrationalR::create_SEXP(mat));

}


//
// return a vector of n boolean corresponding to values that must be affected.
//
std::vector<bool> extract_gmp_R::indice_set_at (unsigned int n , SEXP & IND)
{


  std::vector<int> vidx = bigintegerR::create_int(IND);
  std::vector<bool> result (n,false);


  if(TYPEOF(IND) != NILSXP)
    //LOCICAL
    if (TYPEOF(IND) == LGLSXP)
      {
	for(unsigned int i = 0; i< n; ++i)
	  result[i] = static_cast<bool>( vidx[i % vidx.size() ] );
      }
    else
      //INTEGERS
      {
	//negatives integers: all except indices will be modified
	if (vidx[0] < 0)
	  {
	    for (std::vector<bool>::iterator it = result.begin(); it != result.end(); ++it)
	      *it = true;
	    for (std::vector<int>::const_iterator jt = vidx.begin(); jt != vidx.end(); ++jt)
	      {
		if(*jt > 0)
		  error(_("only 0's may mix with negative subscripts"));
		if( (*jt != 0) && (*jt >= - static_cast<int>(n)))
		  result[-(*jt)-1] = false;
	      }
	  }
	else
	  //INTEGERS (and positive)
	  for (std::vector<int>::const_iterator jt = vidx.begin(); jt != vidx.end(); ++jt)
	    {
	      if(*jt < 0)
		error(_("only 0's may mix with negative subscripts"));

	      if((*jt != 0) && (*jt <= static_cast<int>(n)))
		result[*jt-1] = true;
	    }
      }
  else
    // NILSXP: return true
    for (std::vector<bool>::iterator it = result.begin(); it != result.end(); ++it)
      *it = true;


  return(result);

}//end of indice_set_at

