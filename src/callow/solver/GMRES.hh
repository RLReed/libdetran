//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GMRES.hh
 *  @brief GMRES class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_GMRES_HH_
#define callow_GMRES_HH_

#include "LinearSolver.hh"

namespace callow
{

class MatrixDense;

/**
 *  @class GMRES
 *  @brief Uses preconditioned GMRES(m) iteration to solve a system
 *
 *  GMRES seeks to find the best solution \f$ x \f$
 *  to the linear system
 *  @f[
 *     \mathbf{A}x = b
 *  @f]
 *  such that \f$ x \in \mathcal{K}_n \f$, where the
 *  Krylov subspace is defined
 *  @f[
 *     \mathcal{K}_n \equiv
 *       [b, \mathbf{A}b, \mathbf{A}^2 b, \ldots, \mathbf{A}^{n-1} b] \, .
 *  @f]
 *  Specifically, GMRES finds \f$ x_n \f$ that satisfies
 *  @f[
 *     \min_{x_n \in \mathcal{K}_n} ||b-\mathbf{A}x_n||_2 \, ,
 *  @f]
 *  i.e. it finds the least squares fit within the current Krylov
 *  subspace at every iteration, stopping when that residual
 *  is small enough.
 *
 *  In this implementation, we employ GMRES(m) as described in
 *  Kelley's red book.  A key feature is its use of Givens
 *  rotation for incremental conversion of the upper Hessenberg
 *  matrix \f$ H \f$ to an upper triangle matrix \f$ R \f$.
 *
 */
class GMRES: public LinearSolver
{

public:

  typedef LinearSolver Base;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  GMRES(SP_db db);

  virtual ~GMRES();

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// maximum size of krylov subspace
  int d_restart;

  /// reorthogonalize flag [0 = none, 1 = formula, 2 = always]
  int d_reorthog;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /**
   *  @param b  right hand side
   *  @param x  unknown vector
   */
  void solve_impl(const Vector &b, Vector &x);

  /// apply givens rotation to H
  void apply_givens(MatrixDense &H, Vector &c, Vector &s, const int k);

  void compute_y(MatrixDense &H, Vector &y, const Vector &g, const int k);

};

} // end namespace callow

// Inline member definitions
#include "GMRES.i.hh"

#endif /* callow_GMRES_HH_ */
