!< Define Euler 1D field that is a concrete extension of the abstract integrand type.
module type_euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define Euler 1D field that is a concrete extension of the abstract integrand type.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P
use foodie, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(integrand) :: euler_1D
  !< Euler 1D PDEs system field.
  !<
  !< It is a FOODiE integrand class.
  private
  ! contains
    ! ! public methods
    ! ! auxiliary methods
    ! procedure, pass(self), public :: init             !< Init field.
    ! procedure, pass(self), public :: output           !< Extract Burgers field.
    ! procedure, pass(self), public :: dt => compute_dt !< Compute the current time step, by means of CFL condition.
    ! ! type_integrand deferred methods
    ! procedure, pass(self), public :: t => dBurgers_dt                                        !< Time derivate, residuals function.
    ! procedure, pass(self), public :: update_previous_steps                                    !< Update previous time steps.
    ! procedure, pass(lhs),  public :: integrand_multiply_integrand => burgers_multiply_burgers !< Burgers * burgers operator.
    ! procedure, pass(lhs),  public :: integrand_multiply_real => burgers_multiply_real         !< Burgers * real operator.
    ! procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_burgers         !< Real * Burgers operator.
    ! procedure, pass(lhs),  public :: add => add_burgers                                       !< Burgers + Burgers oprator.
    ! procedure, pass(lhs),  public :: sub => sub_burgers                                       !< Burgers - Burgers oprator.
    ! procedure, pass(lhs),  public :: assign_integrand => burgers_assign_burgers               !< Burgers = Burgers.
    ! procedure, pass(lhs),  public :: assign_real => burgers_assign_real                       !< Burgers = real.
    ! ! operators overloading
    ! generic, public :: operator(-) => sub
    ! ! private methods
    ! procedure, pass(self), private :: x  => dBurgers_dx   !< 1st derivative.
    ! procedure, pass(self), private :: xx => d2Burgers_dx2 !< 2nd derivative.
endtype euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------
! contains
  ! public methods
  ! private methods
endmodule type_euler_1D
