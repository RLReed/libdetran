//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGTransportSolver.hh
 *  @author robertsj
 *  @date   Jun 19, 2012
 *  @brief  MGTransportSolver inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGTRANSPORTSOLVER_I_HH_
#define detran_MGTRANSPORTSOLVER_I_HH_

#include "detran_config.hh"
#include "WGSolverSI.hh"
#ifdef DETRAN_ENABLE_PETSC
#include "WGSolverGMRES.hh"
#endif

namespace detran
{

template <class D>
MGTransportSolver<D>::MGTransportSolver(SP_state                  state,
                                        SP_material               material,
                                        SP_boundary               boundary,
                                        const vec_externalsource &q_e,
                                        SP_fissionsource          q_f,
                                        bool                      multiply)
  : Base(state, material, boundary, q_e, q_f, multiply)
{

  // Get the quadrature from the state
  d_quadrature = d_state;
  Ensure(d_quadrature);

  // Get the inner solver type and create.
  std::string wg_solver = "SI";
  if (d_input->check("inner_solver"))
  {
    wg_solver = d_input->template get<std::string>("inner_solver");
  }
  if (wg_solver == "SI")
  {
    d_wg_solver = new WGSolverSI<D>(d_state, d_material, d_quadrature,
                                    d_boundary, d_externalsources,
                                    d_fissionsource);
  }
  else if (wg_solver == "GMRES")
  {
#ifdef DETRAN_ENABLE_PETSC
    d_wg_solver = new WGSolverGMRES<D>(d_state, d_material, d_quadrature,
                                       d_boundary, d_externalsources,
                                       d_fissionsource);
#else
    THROW("InnerGMRES is not available because PETSc is not enabled.");
#endif
  }
  else
  {
    THROW("Unsupported inner solver type selected: " + wg_solver);
  }

}

} // end namespace detran

#endif /* detran_MGTRANSPORTSOLVER_I_HH_ */
