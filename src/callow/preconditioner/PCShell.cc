//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PCShell.cc
 *  @brief PCShell member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "PCShell.hh"

namespace callow
{

//----------------------------------------------------------------------------//
PCShell::PCShell(std::string name, void* context)
  : Preconditioner(name)
  , d_context(context)
{
  /* ... */
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of PCShell.cc
//----------------------------------------------------------------------------//
