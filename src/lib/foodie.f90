!< FOODIE, Fortran Object oriented Ordinary Differential Equations integration library.

module foodie
!< FOODIE, Fortran Object oriented Ordinary Differential Equations integration library.
!<
!< FOODIE is a KISS library for solving systems of Ordinary Differential Equation (ODE) into the Initial Values Problems (IVP)
!< contest. The mathematical formulation of the problem is:
!<
!< $$U_t = R(t,U)$$
!< $$U_0 = F(0)$$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function and *F* is the (vectorial) initial conditions function.
!<
!< FOODIE is aimed to the time-like integration of the above system of ODE. To this aim, different numerical schemes are provided:
!< see [FOODIE home page](https://github.com/Fortran-FOSS-Programmers/FOODIE) for more details about available integrators.
!<
!<### Usage
!<
!< FOODIE schemes must be applied to only subclass extensions of the abstract class *integrand*.
!<
!< To use FOODIE you must:
!<
!<#### extend integrand abstract class provided by FOODIE implementing your concrete integrand field
!<
!< For example for the Lorenz' ODE system
!<
!<```fortran
!< type, extends(integrand_object) :: lorenz
!<   !< Lorenz equations field.
!<   !<
!<   !< It is a FOODIE integrand class.
!<   private
!<   real(R_P), dimension(:), allocatable :: state        !< Solution vector.
!<   real(R_P)                            :: sigma=0._R_P !< Lorenz \(\sigma\).
!<   real(R_P)                            :: rho=0._R_P   !< Lorenz \(\rho\).
!<   real(R_P)                            :: beta=0._R_P  !< Lorenz \(\beta\).
!<   contains
!<     procedure, pass(self), public :: t => dLorenz_dt                                 !< Time derivate, resiuduals function.
!<     procedure, pass(lhs),  public :: integrand_multiply_real => lorenz_multiply_real !< lorenz * real operator.
!<     procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_lorenz !< Real * Lorenz operator.
!<     procedure, pass(lhs),  public :: add => add_lorenz                               !< Lorenz + Lorenz oprator.
!<     procedure, pass(lhs),  public :: assign_integrand => lorenz_assign_lorenz        !< Lorenz = Lorenz.
!<     procedure, pass(lhs),  public :: assign_real => lorenz_assign_real               !< Lorenz = real.
!<     ...
!< endtype lorenz
!<```
!<
!<#### use one of the provided FOODIE integrator
!<
!< For example using the forward explicit Euler scheme to the above Lorenz' ODE system
!<
!<```fortran
!< use foodie, only : euler_explicit_integrator
!< use type_lorenz, only : lorenz
!< type(euler_explicit_integrator) :: euler_integrator
!< type(lorenz)                    :: attractor
!< real                            :: dt=0.01
!< do step = 1, num_steps
!<   call euler_integrator%integrate(field=attractor, dt=dt)
!< enddo
!<```

use foodie_adt_integrand, only : integrand
use foodie_error_codes, only : error_missing_steps_number
use foodie_integrator_object, only : integrator_object
use foodie_kinds, only : I_P, R_P
use foodie_integrator_adams_bashforth, only : integrator_adams_bashforth
use foodie_integrator_adams_bashforth_moulton, only : integrator_adams_bashforth_moulton
use foodie_integrator_adams_moulton, only : integrator_adams_moulton
use foodie_integrator_backward_differentiation_formula, only : back_df_integrator
use foodie_integrator_emd_runge_kutta, only : emd_runge_kutta_integrator
use foodie_integrator_euler_explicit, only : euler_explicit_integrator
use foodie_integrator_leapfrog, only : leapfrog_integrator
use foodie_integrator_low_storage_runge_kutta, only : ls_runge_kutta_integrator
use foodie_integrator_tvd_runge_kutta, only : tvd_runge_kutta_integrator

implicit none
private
public :: integrand
public :: integrator_object
public :: integrator_adams_bashforth
public :: integrator_adams_bashforth_moulton
public :: integrator_adams_moulton
public :: back_df_integrator
public :: emd_runge_kutta_integrator
public :: euler_explicit_integrator
public :: leapfrog_integrator
public :: ls_runge_kutta_integrator
public :: tvd_runge_kutta_integrator

contains
  function foodie_integrator(scheme, steps) result(integrator)
  !< Return a concrete instance of [[integrator]] given a scheme selection.
  !<
  !< This is the FOODIE integrators factory.
  !<
  !< @note If an error occurs the error status of [[integrator]] is updated.
  character(*), intent(in)              :: scheme     !< Selected integrator given.
  integer(I_P), intent(in), optional    :: steps      !< Number of time steps used in multi-step schemes.
  class(integrator_object), allocatable :: integrator !< The FOODIE integrator.

  select case(trim(adjustl(scheme)))
  case('adams_bashforth')
    allocate(integrator_adams_bashforth :: integrator)
    if (present(steps)) then
      select type(integrator)
      type is(integrator_adams_bashforth)
        call integrator%init(steps=steps)
      endselect
    else
      call integrator%trigger_error(error=ERROR_MISSING_STEPS_NUMBER,                                  &
                                    error_message='missing steps number for initializing integrator!', &
                                    is_severe=.true.)
    endif
  case('adams_bashforth_moulton')
    allocate(integrator_adams_bashforth_moulton :: integrator)
    if (present(steps)) then
      select type(integrator)
      type is(integrator_adams_bashforth_moulton)
        call integrator%init(steps=steps)
      endselect
    else
      call integrator%trigger_error(error=ERROR_MISSING_STEPS_NUMBER,                                  &
                                    error_message='missing steps number for initializing integrator!', &
                                    is_severe=.true.)
    endif
  case('adams_moulton')
    allocate(integrator_adams_moulton :: integrator)
    if (present(steps)) then
      select type(integrator)
      type is(integrator_adams_moulton)
        call integrator%init(steps=steps)
      endselect
    else
      call integrator%trigger_error(error=ERROR_MISSING_STEPS_NUMBER,                                  &
                                    error_message='missing steps number for initializing integrator!', &
                                    is_severe=.true.)
    endif
  endselect
  endfunction foodie_integrator
endmodule foodie
