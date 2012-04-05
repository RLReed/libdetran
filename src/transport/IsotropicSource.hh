//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   IsotropicSource.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  IsotropicSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef ISOTROPICSOURCE_HH_
#define ISOTROPICSOURCE_HH_

#include "ExternalSource.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class IsotropicSource
 * \brief Isotropic volume source.
 */
//===========================================================================//

class IsotropicSource : public ExternalSource
{

public:

  // Source types
  typedef detran_utils::SP<IsotropicSource>  SP_source;
  typedef detran_utils::vec2_dbl             spectra_type;

  IsotropicSource(SP_mesh mesh,
                  SP_quadrature quadrature,
                  int number_groups)
    :  ExternalSource(mesh, quadrature, number_groups)
  { /* ... */ }

  virtual double source(int cell, int group)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    return d_spectra[d_source_map[cell]][group];
  }

  virtual double source(int cell, int group, int angle)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    Require(angle >=0);
    Require(angle < d_number_angles);
    return d_spectra[d_source_map[cell]][group] * d_quadrature->angular_norm();
  }

  void set_source(spectra_type &spectra, Mesh::vec_int &map)
  {
    Require(spectra.size()>0);
    Require(spectra[0].size() == d_number_groups);
    d_source_spectra = spectra;
    d_mesh->add_coarse_mesh_map("ISOTROPICSOURCE", map);
    d_source_map = d_mesh->mesh_map("ISOTROPICSOURCE");
  }

private:

  /// Source spectra
  Mesh::vec_dbl d_source_spectra;

  /// Fine mesh source map
  Mesh::vec_int d_source_map;

};

} // end namespace detran

#endif /* ISOTROPICSOURCE_HH_ */

//---------------------------------------------------------------------------//
//              end of IsotropicSource.hh
//---------------------------------------------------------------------------//