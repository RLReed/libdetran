//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InputHDF5.hh
 * \brief  InputHDF5 
 * \author Jeremy Roberts
 * \date   Jul 29, 2012
 */
//---------------------------------------------------------------------------//

#ifndef INPUTHDF5_HH_
#define INPUTHDF5_HH_

// Utilities
#include "DBC.hh"
#include "InputDB.hh"

// System
#include <string>
#include "hdf5.h"

namespace detran
{

/*!
 *  \class IO_HDF5
 *  \brief Maps an HDF5 database to an InputDB and vice versa.
 *
 *  A problem can be specified completely using the input,
 *  a material definition, and the mesh definition, all of
 *  which can actually live in the InputDB (as done in
 *  \ref StupidParser).  By mapping an HDF5 file to an
 *  input object, we can eliminate the fickle text processing,
 *  using Python or another interface to construct the HDF5
 *  file.  Then, the executable detran can be used, which is
 *  really handy for profiling and other analyses.
 *
 *  For now, we'll use a single group, /input.  All
 *  entries will go in that group.  Later, it might be useful
 *  to add ones for the mesh and material specification.
 */
class IO_HDF5: public Object
{

public:

  typedef InputDB::SP_input SP_input;

  /*!
   *  \brief Constructor
   *  \param filename HDF5 filename
   */
  IO_HDF5(std::string filename);

  /*!
   *  \brief Write the input database into an HDF5 file.
   *  \param input    Input database to be written
   */
  void write(SP_input input);

  /*!
   *  \brief Close any open HDF5 file
   */
  void close();

  /*!
   *  \brief Fill an input database from an HDF5 file.
   *  \param input    Input database to be filled
   *  \param filename HDF5 filename
   */
  void read(SP_input input, std::string filename);

private:

  /// \name Private Data
  /// \{

  /// HDF5 file id
  hid_t   d_file_id;

  /// HDF5 filename
  std::string d_filename;


  /// \}

  /// \name Implementation
  /// \{


  /// \}

};

} // end namespace detran

#endif // INPUTHDF5_HH_ 

//---------------------------------------------------------------------------//
//              end of file InputHDF5.hh
//---------------------------------------------------------------------------//