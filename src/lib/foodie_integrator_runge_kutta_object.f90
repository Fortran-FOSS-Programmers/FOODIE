!< Define the abstract type [[integrator_runge_kutta_object]] of FOODIE ODE integrators.

module foodie_integrator_runge_kutta_object
!< Define the abstract type [[integrator_runge_kutta_object]] of FOODIE ODE integrators.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foodie_integrand_object, only : integrand_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_runge_kutta_object

logical, parameter :: is_multistage_=.true. !< Flag to check if integrator is multistage.
logical, parameter :: is_multistep_=.false. !< Flag to check if integrator is multistep.

type, extends(integrator_object), abstract :: integrator_runge_kutta_object
   !< Abstract type of FOODIE ODE integrators of the family Runge-Kutta.
   contains
      ! deferred methods
      procedure(integrate_interface),      pass(self), deferred :: integrate      !< Integrate integrand field.
      procedure(integrate_fast_interface), pass(self), deferred :: integrate_fast !< Integrate integrand field, fast mode.
      ! deferred methods of parent
      procedure, pass(self) :: is_multistage !< Return .true. for multistage integrator.
      procedure, pass(self) :: is_multistep  !< Return .true. for multistep integrator.
endtype integrator_runge_kutta_object

abstract interface
  !< Abstract interfaces of deferred methods of [[integrator_object]].
  subroutine integrate_interface(self, U, stage, Dt, t)
  !< Integrate integrand field.
  import :: integrand_object, integrator_runge_kutta_object, R_P
  class(integrator_runge_kutta_object), intent(in)    :: self      !< Integrator.
  class(integrand_object),              intent(inout) :: U         !< Integrand.
  class(integrand_object),              intent(inout) :: stage(1:) !< Runge-Kutta stages.
  real(R_P),                            intent(in)    :: Dt        !< Time step.
  real(R_P),                            intent(in)    :: t         !< Time.
  endsubroutine integrate_interface

  subroutine integrate_fast_interface(self, U, stage, buffer, Dt, t)
  !< Integrate integrand field, fast mode.
  import :: integrand_object, integrator_runge_kutta_object, R_P
  class(integrator_runge_kutta_object), intent(in)    :: self      !< Integrator.
  class(integrand_object),              intent(inout) :: U         !< Field to be integrated.
  class(integrand_object),              intent(inout) :: stage(1:) !< Runge-Kutta stages.
  class(integrand_object),              intent(inout) :: buffer    !< Temporary buffer for doing fast operation.
  real(R_P),                            intent(in)    :: Dt        !< Time step.
  real(R_P),                            intent(in)    :: t         !< Time.
  endsubroutine integrate_fast_interface
endinterface

contains
  ! public methods
  elemental function is_multistage(self)
  !< Return .true. for multistage integrator.
  class(integrator_runge_kutta_object), intent(in) :: self          !< Integrator.
  logical                                          :: is_multistage !< Inquire result.

  is_multistage = is_multistage_
  endfunction is_multistage

  elemental function is_multistep(self)
  !< Return .true. for multistage integrator.
  class(integrator_runge_kutta_object), intent(in) :: self         !< Integrator.
  logical                                          :: is_multistep !< Inquire result.

  is_multistep = is_multistep_
  endfunction is_multistep
endmodule foodie_integrator_runge_kutta_object
