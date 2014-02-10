//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Mesh.hh
 *  @brief Mesh class definition.
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_CARTESIANMESH_HH_
#define detran_geometry_CARTESIANMESH_HH_

#include "geometry/geometry_export.hh"
#include "geometry/Mesh.hh"
#include "geometry/TrackDB.hh"
#include "geometry/Point.hh"
#include "utilities/Definitions.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "utilities/Warning.hh"
#include <cmath>
#include <map>
#ifdef DETRAN_ENABLE_BOOST
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#endif

namespace detran_geometry
{

//----------------------------------------------------------------------------//
/**
 *  @class Mesh
 *  @brief Abstract Cartesian mesh class.
 *
 *  Note, the constructors are protected to forbid direct instantiation of
 *  the Mesh class.  Rather, use the dimension-specific subclasses.  We
 *  could use a pure virtual destructor as an alternative.
 */
//----------------------------------------------------------------------------//

class GEOMETRY_EXPORT CartesianMesh: public Mesh
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<CartesianMesh>   SP_cartesianmesh;

  //--------------------------------------------------------------------------//
  // CONSTRUCTORS & DESTRUCTORS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *  @param    dim         Spatial dimension
   *  @param    xfm         Fine meshes per coarse mesh in x dimension.
   *  @param    yfm         Fine meshes per coarse mesh in y dimension.
   *  @param    zfm         Fine meshes per coarse mesh in z dimension.
   *  @param    xcme        Coarse mesh edges x dimension.
   *  @param    ycme        Coarse mesh edges y dimension.
   *  @param    zcme        Coarse mesh edges z dimension.
   */
  CartesianMesh(size_t dim,
                vec_int xfm,  vec_int yfm,  vec_int zfm,
                vec_dbl xcme, vec_dbl ycme, vec_dbl zcme,
                vec_int mat_map);

  /**
   *  @brief Constructor.
   *  @param    dim         Spatial dimension
   *  @param    xfme        Fine mesh edges x dimension.
   *  @param    yfme        Fine mesh edges y dimension.
   *  @param    zfme        Fine mesh edges z dimension.
   */
  CartesianMesh(size_t dim,
                vec_dbl xfme, vec_dbl yfme, vec_dbl zfme,
                vec_int mat_map);

  /// Virtual destructor
  virtual ~CartesianMesh(){}

  //--------------------------------------------------------------------------//
  // Setters
  //--------------------------------------------------------------------------//

  /**
   * @brief  Add map of coarse mesh integer properties.
   *
   * This is an easy way to set mesh properties for meshes based on
   * the coarse mesh regions used to create the mesh.  If the key exists,
   * this function overwrites the map.
   *
   * @param  map_key   String description of map.
   * @param  mesh_map  Logically multi-dimensional map as 1-d vector.
   */
  void add_coarse_mesh_map(const std::string &map_key,
                           const vec_int     &mesh_map);

  //------------------------------------------------------------------------//
  // Getters
  //------------------------------------------------------------------------//

  /// Get number of cells in specified dimension
  //@{
  using Mesh::number_cells;
  size_t number_cells(const size_t dim) const;
  size_t number_cells_x() const;
  size_t number_cells_y() const;
  size_t number_cells_z() const;
  //@}

  /// Get cell width in specified dimension
  //@{
  double width(const size_t dim, const size_t ijk) const;
  double dx(const size_t i) const;
  double dy(const size_t j) const;
  double dz(const size_t k) const;
  //@}

  double volume(const size_t cell) const;

  /// Get domain width along x axis
  double total_width_x() const;
  /// Get domain width along y axis
  double total_width_y() const;
  /// Get domain width along z axis
  double total_width_z() const;

  /**
   * @brief   Returns the cardinal index for i, j, and k
   * @param   i  Index along x axis.
   * @param   j  Index along y axis.
   * @param   k  Index along z axis.
   * @return     Cardinal index.
   */
  size_t index(const size_t i, const size_t j = 0, const size_t k = 0) const;

  /// Return the x, y, or z index given the cardinal index
  //@{
  size_t cell_to_i(const size_t cell) const;
  size_t cell_to_j(const size_t cell) const;
  size_t cell_to_k(const size_t cell) const;
  size_t cell_to_ijk(const size_t dim, const size_t cell) const;
  //@}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MESHES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  /// Find the cell containing a point
  int find_cell(const Point &p) const;

  /// Display some key features
  void display() const;

protected:

  /// x fine meshes in each x coarse mesh
  vec_int d_xfm;
  /// y fine meshes in each y coarse mesh
  vec_int d_yfm;
  /// z fine meshes in each y coarse mesh
  vec_int d_zfm;
  /// x coarse mesh edges
  vec_dbl d_xcme;
  /// y coarse mesh edges
  vec_dbl d_ycme;
  /// z coarse mesh edges
  vec_dbl d_zcme;
  /// x widths
  vec_dbl d_dx;
  /// y widths
  vec_dbl d_dy;
  /// z widths
  vec_dbl d_dz;
  /// total x width
  double d_total_width_x;
  /// total y width
  double d_total_width_y;
  /// total z width
  double d_total_width_z;
  /// Number of cells in x direction
  size_t d_number_cells_x;
  /// Number of cells in y direction
  size_t d_number_cells_y;
  /// Number of cells in y direction
  size_t d_number_cells_z;

private:

  /// Common construction tasks.
  void setup();

#ifdef DETRAN_ENABLE_BOOST

  /// Default constructor.  This is required for serialization and
  /// is of no use otherwise.
  CartesianMesh(){}

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & d_dimension;
    ar & d_xfm;
    ar & d_yfm;
    ar & d_zfm;
    ar & d_xcme;
    ar & d_ycme;
    ar & d_zcme;
    ar & d_dx;
    ar & d_dy;
    ar & d_dz;
    ar & d_total_width_x;
    ar & d_total_width_y;
    ar & d_total_width_z;
    ar & d_number_cells;
    ar & d_number_cells_x;
    ar & d_number_cells_y;
    ar & d_number_cells_z;
    ar & d_mesh_map;
  }

#endif

};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<CartesianMesh>)

//----------------------------------------------------------------------------//
struct GEOMETRY_EXPORT Mesh1D
{
  typedef detran_utilities::vec_dbl vec_dbl;
  typedef detran_utilities::vec_int vec_int;
  static CartesianMesh::SP_cartesianmesh
  Create(vec_int xfm, vec_dbl xcm, vec_int mat_map)
  {
    vec_int yfm(1, 1);
    vec_dbl ycm(2, 0.0);
    ycm[1] = 1.0;
    CartesianMesh::SP_cartesianmesh
    p(new CartesianMesh(1, xfm, yfm, yfm, xcm, ycm, ycm, mat_map));
    return p;
  }
};

//----------------------------------------------------------------------------//
struct GEOMETRY_EXPORT Mesh2D
{
  typedef detran_utilities::vec_dbl vec_dbl;
  typedef detran_utilities::vec_int vec_int;
  static CartesianMesh::SP_cartesianmesh
  Create(vec_int xfm, vec_int yfm, vec_dbl xcm, vec_dbl ycm, vec_int mat_map)
  {
    vec_int zfm(1, 1);
    vec_dbl zcm(2, 0.0);
    zcm[1] = 1.0;
    CartesianMesh::SP_cartesianmesh
    p(new CartesianMesh(2, xfm, yfm, zfm, xcm, ycm, zcm, mat_map));
    return p;
  }
};

//----------------------------------------------------------------------------//
struct GEOMETRY_EXPORT Mesh3D
{
  typedef detran_utilities::vec_dbl vec_dbl;
  typedef detran_utilities::vec_int vec_int;
  static CartesianMesh::SP_cartesianmesh
  Create(vec_int xfm, vec_int yfm, vec_int zfm,
         vec_dbl xcm, vec_dbl ycm, vec_dbl zcm,
         vec_int mat_map)
  {
    CartesianMesh::SP_cartesianmesh
    p(new CartesianMesh(3, xfm, yfm, zfm, xcm, ycm, zcm, mat_map));
    return p;
  }
};

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#include "CartesianMesh.i.hh"

#endif /* detran_geometry_CARTESIANMESH_HH_ */

//----------------------------------------------------------------------------//
//              end of Mesh.hh
//----------------------------------------------------------------------------//
