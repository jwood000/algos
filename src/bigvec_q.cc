#include "bigvec_q.h"



bigvec_q::bigvec_q(const bigvec_q & rhs):
  value(rhs.value.size()),
  nrow(0)
{
  *this = rhs;
}


bigvec_q::bigvec_q(const bigvec & rhs):
  value(rhs.value.size()),
  nrow(rhs.nrow)
{
  for(unsigned int i=0; i< rhs.size(); ++i)
    {
      value[i].setValue(rhs.value[i]);
    }
}

bigvec_q & bigvec_q::operator= (const bigvec_q & rhs)

{
  if(this != &rhs)
    {
      nrow = rhs.nrow;
      value.resize(rhs.value.size());
      std::vector<bigrational>::iterator it = value.begin();
      std::vector<bigrational>::const_iterator jt = rhs.value.begin();
      while(it != value.end())
	{
	  *it = *jt;
	  ++it;
	  ++jt;
	}
    }
  return(*this);

}


bigrational  bigvec_q::operator[] (unsigned int i) const
{
  return(value[i]);
}

void bigvec_q::set(unsigned int i,const bigrational & val)
{
  //DEBUG !!
  if(i>=value.size())
    {
      Rprintf("t nul a bigvec_q_set\n");
      return;
    }
  value[i] = val;
}
void bigvec_q::set(unsigned int i,const mpq_t & val)
{

  if(i>=value.size())
    {
      Rprintf("t nul a bigvec_q_set_mpq __LINE__ \n");
      return;
    }

  value[i].setValue(val);

}

void bigvec_q::push_back(const bigrational &i)
{
  value.push_back(i);
}

unsigned int bigvec_q::size() const
{
  return(value.size());
}

void bigvec_q::resize(unsigned int n)
{
  value.resize(n);
}

void bigvec_q::clear()
{
  value.clear();
  nrow=0;
}



//
// \brief substract lambda[0] * line j to line i
//
void bigvec_q::subLine(unsigned int i,unsigned int j,bigvec_q lambda)
{
  if(nrow <= 0)
    error(_("Need matrix with at least one row to do this operation"));

  unsigned int k, n = value.size() /  nrow;
  for(k=0; k < n; ++k)
      value[i + k*nrow] =  value[i + k*nrow] - ( value[j + k*nrow] * lambda.value[0] ) ;
}


/*
 * \brief multiply line i by lambda[0]
 */
void bigvec_q::mulLine(unsigned int i, bigvec_q lambda)
{
  if(nrow <= 0)
    error(_("Need matrix with at least one row to do this operation"));

  unsigned int k, n = value.size() / nrow;
  for(k=0; k < n; ++k)
      value[i + k*nrow] =  value[i + k*nrow]  * lambda.value[0];
}

// never used
void bigvec_q::print()
{
  if(nrow>0) {
    for(int i=0; i < nrow; ++i)
      {
	for(unsigned int j=0; j < (value.size() / nrow); ++j)
	  Rprintf("%s\t", value[i+j* nrow].str(10).c_str() );
	Rprintf("\n");
      }
  }
  else {
    for(unsigned int i=0; i < value.size(); ++i)
      Rprintf("%s\t", value[i].str(10).c_str() );
    Rprintf("\n");
  }
}
