!< Define the abstract type [[integrator_multistage_explicit_object]] of FOODIE ODE integrators.

module foodie_integrator_multistage_explicit_object
!< Define the abstract type [[integrator_multistage_explicit_object]] of FOODIE ODE integrators.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foodie_integrand_object, only : integrand_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_multistage_explicit_object

type, extends(integrator_object), abstract :: integrator_multistage_explicit_object
   !< Abstract type of FOODIE ODE integrators of the multistage family.
   integer(I_P) :: stages=0 !< Number of stages.
   contains
      ! deferred methods
      procedure(integrate_interface),      pass(self), deferred :: integrate      !< Integrate integrand field.
      procedure(integrate_fast_interface), pass(self), deferred :: integrate_fast !< Integrate integrand field, fast mode.
      ! implemented deferred methods of parent
      procedure, pass(self) :: is_multistage !< Return .true. for multistage integrator.
      procedure, pass(self) :: is_multistep  !< Return .true. for multistep integrator.
      procedure, pass(self) :: stages_number !< Return number of stages used.
      procedure, pass(self) :: steps_number  !< Return number of steps used.
      ! public methods
      procedure, pass(self) :: destroy_multistage !< Destroy the integrator.
endtype integrator_multistage_explicit_object

abstract interface
  !< Abstract interfaces of deferred methods of [[integrator_object]].
  subroutine integrate_interface(self, U, stage, Dt, t, new_Dt)
  !< Integrate integrand field.
  import :: integrand_object, integrator_multistage_explicit_object, R_P
  class(integrator_multistage_explicit_object), intent(in)    :: self      !< Integrator.
  class(integrand_object),                      intent(inout) :: U         !< Integrand.
  class(integrand_object),                      intent(inout) :: stage(1:) !< Runge-Kutta stages.
  real(R_P),                                    intent(in)    :: Dt        !< Time step.
  real(R_P),                                    intent(in)    :: t         !< Time.
  real(R_P), optional,                          intent(out)   :: new_Dt    !< New adapted time step.
  endsubroutine integrate_interface

  subroutine integrate_fast_interface(self, U, stage, buffer, Dt, t, new_Dt)
  !< Integrate integrand field, fast mode.
  import :: integrand_object, integrator_multistage_explicit_object, R_P
  class(integrator_multistage_explicit_object), intent(in)    :: self      !< Integrator.
  class(integrand_object),                      intent(inout) :: U         !< Field to be integrated.
  class(integrand_object),                      intent(inout) :: stage(1:) !< Runge-Kutta stages.
  class(integrand_object),                      intent(inout) :: buffer    !< Temporary buffer for doing fast operation.
  real(R_P),                                    intent(in)    :: Dt        !< Time step.
  real(R_P),                                    intent(in)    :: t         !< Time.
  real(R_P), optional,                          intent(out)   :: new_Dt    !< New adapted time step.
  endsubroutine integrate_fast_interface
endinterface

contains
   ! deferred methods
   elemental function is_multistage(self)
   !< Return .true. for multistage integrator.
   class(integrator_multistage_explicit_object), intent(in) :: self          !< Integrator.
   logical                                                  :: is_multistage !< Inquire result.

   is_multistage = .true.
   endfunction is_multistage

   elemental function is_multistep(self)
   !< Return .true. for multistage integrator.
   class(integrator_multistage_explicit_object), intent(in) :: self         !< Integrator.
   logical                                                  :: is_multistep !< Inquire result.

   is_multistep = .false.
   endfunction is_multistep

   elemental function stages_number(self)
   !< Return number of stages used.
   class(integrator_multistage_explicit_object), intent(in) :: self          !< Integrator.
   integer(I_P)                                             :: stages_number !< Number of stages used.

   stages_number = self%stages
   endfunction stages_number

   elemental function steps_number(self)
   !< Return number of steps used.
   class(integrator_multistage_explicit_object), intent(in) :: self         !< Integrator.
   integer(I_P)                                             :: steps_number !< Number of steps used.

   steps_number = 0
   endfunction steps_number

   ! public methods
   elemental subroutine destroy_multistage(self)
   !< Destroy the integrator.
   class(integrator_multistage_explicit_object), intent(inout) :: self !< Integrator.

   call self%destroy_abstract
   self%stages = 0
   endsubroutine destroy_multistage
endmodule foodie_integrator_multistage_explicit_object
