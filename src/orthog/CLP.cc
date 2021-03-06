//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CLP.cc
 *  @brief  CLP
 *  @author Jeremy Roberts
 *  @date   Jan 8, 2013
 */
//---------------------------------------------------------------------------//

#include "CLP.hh"
#ifdef DETRAN_ENABLE_BOOST
#include <boost/math/special_functions/legendre.hpp>
#endif

namespace detran_orthog
{

//---------------------------------------------------------------------------//
CLP::CLP(const size_t   order,
         const vec_dbl &x,
         const vec_dbl &qw,
         const double   x_0,
         const double   x_1)
  : ContinuousOrthogonalBasis(order, x, qw)
{
#ifndef DETRAN_ENABLE_BOOST
  THROW("CLP needs boost to be enabled.");
#else
  Require(x_1 > x_0);

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_order + 1, d_size, 0.0);

  // Allocate the normalization array
  d_a = Vector::Create(d_order + 1, 0.0);

  // The weights are just the qw's.
  double L = x_1 - x_0;
  for (size_t i = 0; i < d_w->size(); ++i)
  {
    d_x[i] =  2.0*(d_x[i] - x_0)/L - 1.0;
    (*d_w)[i] = d_qw[i];
  }

  // Build the basis
  for (size_t l = 0; l <= d_order; ++l)
  {
    for (size_t i = 0; i < d_size; ++i)
    {
      (*d_basis)(l, i) = boost::math::legendre_p(l, d_x[i]);
    }
    // Inverse of normalization coefficient.
    (*d_a)[l] = (2.0 * l + 1.0) / 2.0;
  }
  d_orthonormal = true;
  compute_a();
#endif
}

} // end namespace detran_orthog

//---------------------------------------------------------------------------//
//              end of file CLP.cc
//---------------------------------------------------------------------------//
