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
!< use foodie, only : integrator_euler_explicit
!< use type_lorenz, only : lorenz
!< type(integrator_euler_explicit) :: euler_integrator
!< type(lorenz)                    :: attractor
!< real                            :: dt=0.01
!< do step = 1, num_steps
!<   call euler_integrator%integrate(field=attractor, dt=dt)
!< enddo
!<```

use foodie_adt_integrand, only : integrand
use foodie_error_codes, only : ERROR_MISSING_STAGES_NUMBER, ERROR_MISSING_STEPS_NUMBER
use foodie_integrator_object, only : integrator_object
use foodie_kinds, only : I_P, R_P
use foodie_integrator_adams_bashforth, only : integrator_adams_bashforth
use foodie_integrator_adams_bashforth_moulton, only : integrator_adams_bashforth_moulton
use foodie_integrator_adams_moulton, only : integrator_adams_moulton
use foodie_integrator_backward_differentiation_formula, only : integrator_back_df
use foodie_integrator_euler_explicit, only : integrator_euler_explicit
use foodie_integrator_leapfrog, only : integrator_leapfrog
use foodie_integrator_lmm_ssp, only : integrator_lmm_ssp
use foodie_integrator_lmm_ssp_vss, only : integrator_lmm_ssp_vss
use foodie_integrator_runge_kutta_emd, only : integrator_runge_kutta_emd
use foodie_integrator_runge_kutta_low_storage, only : integrator_runge_kutta_ls
use foodie_integrator_runge_kutta_tvd, only : integrator_runge_kutta_tvd

implicit none
private
public :: integrand
public :: integrator_object
public :: integrator_adams_bashforth
public :: integrator_adams_bashforth_moulton
public :: integrator_adams_moulton
public :: integrator_back_df
public :: integrator_euler_explicit
public :: integrator_leapfrog
public :: integrator_lmm_ssp
public :: integrator_lmm_ssp_vss
public :: integrator_runge_kutta_emd
public :: integrator_runge_kutta_ls
public :: foodie_integrator
public :: integrator_runge_kutta_tvd

contains
  function foodie_integrator(scheme, steps, stages, order, tolerance, nu, alpha) result(integrator)
  !< Return a concrete instance of [[integrator_object]] given a scheme selection.
  !<
  !< This is the FOODIE integrators factory.
  !<
  !< @note If an error occurs the error status of [[integrator_object]] is updated.
  character(*), intent(in)              :: scheme     !< Selected integrator given.
  integer(I_P), intent(in), optional    :: steps      !< Number of time steps used in multi-step schemes.
  integer(I_P), intent(in), optional    :: stages     !< Number of Runge-Kutta stages used in multi-stage schemes.
  integer(I_P), intent(in), optional    :: order      !< Order of accuracy.
  real(R_P),    intent(in), optional    :: tolerance  !< Tolerance on the local truncation error.
  real(R_P),    intent(in), optional    :: nu         !< Williams-Robert-Asselin filter coefficient.
  real(R_P),    intent(in), optional    :: alpha      !< Robert-Asselin filter coefficient.
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
  case('back_df')
    allocate(integrator_back_df :: integrator)
    if (present(steps)) then
      select type(integrator)
      type is(integrator_back_df)
        call integrator%init(steps=steps)
      endselect
    else
      call integrator%trigger_error(error=ERROR_MISSING_STEPS_NUMBER,                                  &
                                    error_message='missing steps number for initializing integrator!', &
                                    is_severe=.true.)
    endif
  case('euler_explicit')
    allocate(integrator_euler_explicit :: integrator)
  case('leapfrog')
    allocate(integrator_leapfrog :: integrator)
    select type(integrator)
    type is(integrator_leapfrog)
      call integrator%init(nu=nu, alpha=alpha)
    endselect
  case('lmm_ssp')
    allocate(integrator_lmm_ssp :: integrator)
    if (present(steps)) then
      select type(integrator)
      type is(integrator_lmm_ssp)
        call integrator%init(steps=steps)
      endselect
    else
      call integrator%trigger_error(error=ERROR_MISSING_STEPS_NUMBER,                                  &
                                    error_message='missing steps number for initializing integrator!', &
                                    is_severe=.true.)
    endif
  case('lmm_ssp_vss')
    allocate(integrator_lmm_ssp_vss :: integrator)
    if ((.not.present(steps)).or.(.not.present(order))) then
      call integrator%trigger_error(error=ERROR_MISSING_STEPS_NUMBER,                                  &
                                    error_message='missing steps number for initializing integrator!', &
                                    is_severe=.true.)
    else
      select type(integrator)
      type is(integrator_lmm_ssp_vss)
        call integrator%init(steps=steps, order=order)
      endselect
    endif
  case('runge_kutta_emd')
    allocate(integrator_runge_kutta_emd :: integrator)
    if (present(stages)) then
      select type(integrator)
      type is(integrator_runge_kutta_emd)
        call integrator%init(stages=stages, tolerance=tolerance)
      endselect
    else
      call integrator%trigger_error(error=ERROR_MISSING_STAGES_NUMBER,                                  &
                                    error_message='missing stages number for initializing integrator!', &
                                    is_severe=.true.)
    endif
  case('runge_kutta_ls')
    allocate(integrator_runge_kutta_ls :: integrator)
    if (present(stages)) then
      select type(integrator)
      type is(integrator_runge_kutta_ls)
        call integrator%init(stages=stages)
      endselect
    else
      call integrator%trigger_error(error=ERROR_MISSING_STAGES_NUMBER,                                  &
                                    error_message='missing stages number for initializing integrator!', &
                                    is_severe=.true.)
    endif
  case('runge_kutta_tvd')
    allocate(integrator_runge_kutta_tvd :: integrator)
    if (present(stages)) then
      select type(integrator)
      type is(integrator_runge_kutta_tvd)
        call integrator%init(stages=stages)
      endselect
    else
      call integrator%trigger_error(error=ERROR_MISSING_STAGES_NUMBER,                                  &
                                    error_message='missing stages number for initializing integrator!', &
                                    is_severe=.true.)
    endif
  endselect
  endfunction foodie_integrator
endmodule foodie
