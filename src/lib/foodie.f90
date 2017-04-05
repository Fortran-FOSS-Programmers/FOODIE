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

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foodie_adt_integrand, only : integrand
use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
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
use foodie_integrator_runge_kutta_ssp, only : integrator_runge_kutta_ssp

implicit none
private
public :: foodie_integrator
public :: foodie_integrator_class_names
public :: foodie_integrator_schemes
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
public :: integrator_runge_kutta_ssp
public :: is_available
public :: is_class_available
public :: is_scheme_available

contains
  function foodie_integrator(scheme, tolerance, nu, alpha) result(integrator)
  !< Return a concrete instance of [[integrator_object]] given a scheme selection.
  !<
  !< This is the FOODIE integrators factory.
  !<
  !< @note If an error occurs the error status of [[integrator_object]] is updated.
  character(*), intent(in)                 :: scheme                      !< Selected integrator given.
  real(R_P),    intent(in), optional       :: tolerance                   !< Tolerance on the local truncation error.
  real(R_P),    intent(in), optional       :: nu                          !< Williams-Robert-Asselin filter coefficient.
  real(R_P),    intent(in), optional       :: alpha                       !< Robert-Asselin filter coefficient.
  class(integrator_object), allocatable    :: integrator                  !< The FOODIE integrator.
  type(integrator_adams_bashforth)         :: int_adams_bashforth         !< Integrator Adams Bashforth.
  type(integrator_adams_bashforth_moulton) :: int_adams_bashforth_moulton !< Integrator Adams Bashforth Moulton.
  type(integrator_adams_moulton)           :: int_adams_moulton           !< Integrator Adams Moulton.
  type(integrator_back_df)                 :: int_back_df                 !< Integrator back differentiation formula.
  type(integrator_euler_explicit)          :: int_euler_explicit          !< Integrator euler explicit.
  type(integrator_leapfrog)                :: int_leapfrog                !< Integrator leapfrog.
  type(integrator_lmm_ssp)                 :: int_lmm_ssp                 !< Integrator lmm ssp.
  type(integrator_lmm_ssp_vss)             :: int_lmm_ssp_vss             !< Integrator lmm ssp_vss.
  type(integrator_runge_kutta_emd)         :: int_runge_kutta_emd         !< Integrator Runge Kutta_embdedded.
  type(integrator_runge_kutta_ls)          :: int_runge_kutta_ls          !< Integrator Runge Kutta low storage.
  type(integrator_runge_kutta_ssp)         :: int_runge_kutta_ssp         !< Integrator Runge Kutta ssp.

  if     (index(trim(adjustl(scheme)), trim(int_adams_bashforth_moulton%class_name())) > 0) then
    allocate(integrator_adams_bashforth_moulton :: integrator)
    select type(integrator)
    type is(integrator_adams_bashforth_moulton)
      call integrator%initialize(scheme=scheme)
    endselect
  elseif (index(trim(adjustl(scheme)), trim(int_adams_bashforth%class_name())) > 0) then
    allocate(integrator_adams_bashforth :: integrator)
    select type(integrator)
    type is(integrator_adams_bashforth)
      call integrator%initialize(scheme=scheme)
    endselect
  elseif (index(trim(adjustl(scheme)), trim(int_adams_moulton%class_name())) > 0) then
    allocate(integrator_adams_moulton :: integrator)
    select type(integrator)
    type is(integrator_adams_moulton)
      call integrator%initialize(scheme=scheme)
    endselect
  elseif (index(trim(adjustl(scheme)), trim(int_back_df%class_name())) > 0) then
    allocate(integrator_back_df :: integrator)
    select type(integrator)
    type is(integrator_back_df)
      call integrator%initialize(scheme=scheme)
    endselect
  elseif (index(trim(adjustl(scheme)), trim(int_euler_explicit%class_name())) > 0) then
    allocate(integrator_euler_explicit :: integrator)
  elseif (index(trim(adjustl(scheme)), trim(int_leapfrog%class_name())) > 0) then
    allocate(integrator_leapfrog :: integrator)
    select type(integrator)
    type is(integrator_leapfrog)
      call integrator%initialize(scheme=scheme, nu=nu, alpha=alpha)
    endselect
  elseif (index(trim(adjustl(scheme)), trim(int_lmm_ssp_vss%class_name())) > 0) then
    allocate(integrator_lmm_ssp_vss :: integrator)
    select type(integrator)
    type is(integrator_lmm_ssp_vss)
      call integrator%initialize(scheme=scheme)
    endselect
  elseif (index(trim(adjustl(scheme)), trim(int_lmm_ssp%class_name())) > 0) then
    allocate(integrator_lmm_ssp :: integrator)
    select type(integrator)
    type is(integrator_lmm_ssp)
      call integrator%initialize(scheme=scheme)
    endselect
  elseif (index(trim(adjustl(scheme)), trim(int_runge_kutta_emd%class_name())) > 0) then
    allocate(integrator_runge_kutta_emd :: integrator)
    select type(integrator)
    type is(integrator_runge_kutta_emd)
      call integrator%initialize(scheme=scheme, tolerance=tolerance)
    endselect
  elseif (index(trim(adjustl(scheme)), trim(int_runge_kutta_ls%class_name())) > 0) then
    allocate(integrator_runge_kutta_ls :: integrator)
    select type(integrator)
    type is(integrator_runge_kutta_ls)
      call integrator%initialize(scheme=scheme)
    endselect
  elseif (index(trim(adjustl(scheme)), trim(int_runge_kutta_ssp%class_name())) > 0) then
    allocate(integrator_runge_kutta_ssp :: integrator)
    select type(integrator)
    type is(integrator_runge_kutta_ssp)
      call integrator%initialize(scheme=scheme)
    endselect
  else
    write(stderr, '(A)')'error: "'//trim(adjustl(scheme))//'" scheme is unknown!'
    stop
  endif
  endfunction foodie_integrator

  pure function foodie_integrator_class_names() result(names)
  !< Return the list of available intergrator class of schemes names.
  character(len=99), allocatable           :: names(:)                    !< Available integrator class names.
  type(integrator_adams_bashforth)         :: int_adams_bashforth         !< Integrator Adams Bashforth.
  type(integrator_adams_bashforth_moulton) :: int_adams_bashforth_moulton !< Integrator Adams Bashforth Moulton.
  type(integrator_adams_moulton)           :: int_adams_moulton           !< Integrator Adams Moulton.
  type(integrator_back_df)                 :: int_back_df                 !< Integrator back differentiation formula.
  type(integrator_euler_explicit)          :: int_euler_explicit          !< Integrator euler explicit.
  type(integrator_leapfrog)                :: int_leapfrog                !< Integrator leapfrog.
  type(integrator_lmm_ssp)                 :: int_lmm_ssp                 !< Integrator lmm ssp.
  type(integrator_lmm_ssp_vss)             :: int_lmm_ssp_vss             !< Integrator lmm ssp_vss.
  type(integrator_runge_kutta_emd)         :: int_runge_kutta_emd         !< Integrator Runge Kutta_embdedded.
  type(integrator_runge_kutta_ls)          :: int_runge_kutta_ls          !< Integrator Runge Kutta low storage.
  type(integrator_runge_kutta_ssp)         :: int_runge_kutta_ssp         !< Integrator Runge Kutta ssp.

  names = [       int_adams_bashforth         % class_name()]
  names = [names, int_adams_bashforth_moulton % class_name()]
  names = [names, int_adams_moulton           % class_name()]
  names = [names, int_back_df                 % class_name()]
  names = [names, int_euler_explicit          % class_name()]
  names = [names, int_leapfrog                % class_name()]
  names = [names, int_lmm_ssp                 % class_name()]
  names = [names, int_lmm_ssp_vss             % class_name()]
  names = [names, int_runge_kutta_emd         % class_name()]
  names = [names, int_runge_kutta_ls          % class_name()]
  names = [names, int_runge_kutta_ssp         % class_name()]
  endfunction foodie_integrator_class_names

  pure function foodie_integrator_schemes(class_name) result(schemes)
  !< Return the list of all available intergrator schemes, or only the schemes belonging to the given class name.
  !<
  !< @note In the case the class name passed is not available the schemes array is returned not allocated and not warnings are
  !< given: it is up to the caller to check the status of the schemes array.
  character(*), intent(in), optional       :: class_name                  !< Return only the schemes belogn to the given class.
  character(len=99), allocatable           :: schemes(:)                  !< Available integrators.
  type(integrator_adams_bashforth)         :: int_adams_bashforth         !< Integrator Adams Bashforth.
  type(integrator_adams_bashforth_moulton) :: int_adams_bashforth_moulton !< Integrator Adams Bashforth Moulton.
  type(integrator_adams_moulton)           :: int_adams_moulton           !< Integrator Adams Moulton.
  type(integrator_back_df)                 :: int_back_df                 !< Integrator back differentiation formula.
  type(integrator_euler_explicit)          :: int_euler_explicit          !< Integrator euler explicit.
  type(integrator_leapfrog)                :: int_leapfrog                !< Integrator leapfrog.
  type(integrator_lmm_ssp)                 :: int_lmm_ssp                 !< Integrator lmm ssp.
  type(integrator_lmm_ssp_vss)             :: int_lmm_ssp_vss             !< Integrator lmm ssp_vss.
  type(integrator_runge_kutta_emd)         :: int_runge_kutta_emd         !< Integrator Runge Kutta_embdedded.
  type(integrator_runge_kutta_ls)          :: int_runge_kutta_ls          !< Integrator Runge Kutta low storage.
  type(integrator_runge_kutta_ssp)         :: int_runge_kutta_ssp         !< Integrator Runge Kutta ssp.

  if (present(class_name)) then
    if     (trim(int_adams_bashforth%class_name()) == trim(adjustl(class_name))) then
      schemes = int_adams_bashforth%supported_schemes()
    elseif (trim(int_adams_bashforth_moulton%class_name()) == trim(adjustl(class_name))) then
      schemes = int_adams_bashforth_moulton%supported_schemes()
    elseif (trim(int_adams_moulton%class_name()) == trim(adjustl(class_name))) then
      schemes = int_adams_moulton%supported_schemes()
    elseif (trim(int_back_df%class_name()) == trim(adjustl(class_name))) then
      schemes = int_back_df%supported_schemes()
    elseif (trim(int_euler_explicit%class_name()) == trim(adjustl(class_name))) then
      schemes = int_euler_explicit%supported_schemes()
    elseif (trim(int_leapfrog%class_name()) == trim(adjustl(class_name))) then
      schemes = int_leapfrog%supported_schemes()
    elseif (trim(int_lmm_ssp%class_name()) == trim(adjustl(class_name))) then
      schemes = int_lmm_ssp%supported_schemes()
    elseif (trim(int_lmm_ssp_vss%class_name()) == trim(adjustl(class_name))) then
      schemes = int_lmm_ssp_vss%supported_schemes()
    elseif (trim(int_runge_kutta_emd%class_name()) == trim(adjustl(class_name))) then
      schemes = int_runge_kutta_emd%supported_schemes()
    elseif (trim(int_runge_kutta_ls%class_name()) == trim(adjustl(class_name))) then
      schemes = int_runge_kutta_ls%supported_schemes()
    elseif (trim(int_runge_kutta_ssp%class_name()) == trim(adjustl(class_name))) then
      schemes = int_runge_kutta_ssp%supported_schemes()
    endif
  else
    schemes =           int_adams_bashforth         % supported_schemes()
    schemes = [schemes, int_adams_bashforth_moulton % supported_schemes()]
    schemes = [schemes, int_adams_moulton           % supported_schemes()]
    schemes = [schemes, int_back_df                 % supported_schemes()]
    schemes = [schemes, int_euler_explicit          % supported_schemes()]
    schemes = [schemes, int_leapfrog                % supported_schemes()]
    schemes = [schemes, int_lmm_ssp                 % supported_schemes()]
    schemes = [schemes, int_lmm_ssp_vss             % supported_schemes()]
    schemes = [schemes, int_runge_kutta_emd         % supported_schemes()]
    schemes = [schemes, int_runge_kutta_ls          % supported_schemes()]
    schemes = [schemes, int_runge_kutta_ssp         % supported_schemes()]
  endif
  endfunction foodie_integrator_schemes

  pure function is_available(scheme)
  !< Return .true. if the given scheme (or class of schemes name) is available in the FOODIE library.
  character(*), intent(in)   :: scheme       !< Selected integrator given.
  logical                    :: is_available !< Availability result.

  is_available = is_class_available(scheme=scheme)
  if (is_available) return
  is_available = is_scheme_available(scheme=scheme)
  if (is_available) return
  endfunction is_available

  pure function is_class_available(scheme)
  !< Return .true. if the given class of schemes name is available.
  character(*), intent(in)   :: scheme                    !< Selected integrator given.
  logical                    :: is_class_available        !< Availability result.
  character(99), allocatable :: integrator_class_names(:) !< Name of FOODIE integrator classes.
  integer(I_P)               :: i                         !< Counter.

  integrator_class_names = foodie_integrator_class_names()
  is_class_available = .false.
  do i=1, size(integrator_class_names, dim=1)
    is_class_available = (trim(adjustl(scheme))==trim(adjustl(integrator_class_names(i))))
    if (is_class_available) return
  enddo
  endfunction is_class_available

  pure function is_scheme_available(scheme)
  !< Return .true. if the given scheme (class name) is available in the FOODIE library.
  character(*), intent(in)   :: scheme                !< Selected integrator given.
  logical                    :: is_scheme_available   !< Availability result.
  character(99), allocatable :: integrator_schemes(:) !< Name of FOODIE integrator schemes.
  integer(I_P)               :: i                     !< Counter.

  integrator_schemes = foodie_integrator_schemes()
  is_scheme_available = .false.
  do i=1, size(integrator_schemes, dim=1)
    is_scheme_available = (trim(adjustl(scheme))==trim(adjustl(integrator_schemes(i))))
    if (is_scheme_available) return
  enddo
  endfunction is_scheme_available
endmodule foodie
