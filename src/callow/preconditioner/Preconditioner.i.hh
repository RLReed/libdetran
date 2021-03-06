//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Preconditioner.i.hh
 *  @brief  Preconditioner inline member definitions
 *  @author Jeremy Roberts
 *  @date   Sep 20, 2012
 */
//---------------------------------------------------------------------------//

#ifndef callow_PRECONDITIONER_I_HH_
#define callow_PRECONDITIONER_I_HH_

namespace callow
{

#ifdef DETRAN_ENABLE_PETSC

//---------------------------------------------------------------------------//
inline void Preconditioner::set_petsc_pc(PC pc)
{
  // save the pc object
  d_petsc_pc = pc;

  PetscErrorCode ierr;
  // set the shell preconditioner
  ierr = PCSetType(pc, PCSHELL);
  // set the PC context
  ierr = PCShellSetContext(d_petsc_pc, this);
  Insist(!ierr, "Error setting shell preconditioner context.");
  // set the PC operator
  ierr = PCShellSetApply(d_petsc_pc, pc_apply_wrapper);
  Insist(!ierr, "Error setting shell preconditioner operator.");
  // set the PC name for good measure
  ierr = PCShellSetName(d_petsc_pc, d_name.c_str());
  Insist(!ierr, "Error within-group preconditioner name.");

  Ensure(!ierr);
}

//---------------------------------------------------------------------------//
inline PetscErrorCode pc_apply_wrapper(PC pc, Vec b, Vec x)
{
  // get the context and cast
  PetscErrorCode ierr;
  void *context;
  ierr = PCShellGetContext(pc, &context); CHKERRQ(ierr);
  Preconditioner *foo = (Preconditioner *) context;
  // wrap the petsc vectors
  Vector B(b);
  Vector X(x);
//  B.display();
//  X.display();
  // call the actual apply operator.
  foo->apply(B, X);
  return ierr;
}

#endif

} // end namespace detran

#endif // callow_PRECONDITIONER_I_HH_

//---------------------------------------------------------------------------//
//              end of file Preconditioner.i.hh
//---------------------------------------------------------------------------//
