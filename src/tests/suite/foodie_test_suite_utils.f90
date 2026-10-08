!< Utilities of the FOODIE tests suite: test problems factory and expected orders of the schemes.

module foodie_test_suite_utils
!< Utilities of the FOODIE tests suite: test problems factory and expected orders of the schemes.
!<
!< The expected orders are the formal ones, with the following exceptions documented here once for all tests:
!<
!< + `runge_kutta_emd_stages_7_order_4` (Dormand-Prince) advances the solution with its 5th order formula;
!< + `leapfrog_raw` (Robert-Asselin-Williams filtered leapfrog, fixed filter coefficient) is 1st order accurate, the 3rd
!<   order accuracy of [Williams 2011] being referred to the amplitude of oscillations;
!< + `runge_kutta_lssp_*` (linear SSP) schemes are *s*, *s-1* order accurate only for linear autonomous problems, being
!<   (at most) 2nd order accurate for non autonomous ones.

use flap, only : command_line_interface
use foodie, only : foodie_integrator_factory, integrator_object
use foodie_test_integrand_lcce, only : integrand_lcce
use foodie_test_integrand_oscillation, only : integrand_oscillation
use foodie_test_integrand_polynomial, only : integrand_polynomial
use foodie_test_integrand_riccati, only : integrand_riccati
use foodie_test_integrand_tester_object, only : integrand_tester_object
use penf, only : I_P, R_P

implicit none
private
public :: convergence_order
public :: convergence_problem
public :: has_fast_mode
public :: is_linear_multistep
public :: make_problem
public :: polynomial_degree
public :: STAGES

integer(I_P), parameter :: STAGES = 3_I_P !< Stages number of the schemes with selectable stages (linear SSP).

contains
   function make_problem(problem, args) result(integrand)
   !< Return the initial conditions of a test problem, set by means of its own CLI (args are its options).
   character(*), intent(in)                    :: problem   !< Problem name: lcce, oscillation, polynomial, riccati.
   character(*), intent(in), optional          :: args      !< Problem options.
   class(integrand_tester_object), allocatable :: integrand !< Initial conditions.
   type(command_line_interface)                :: cli       !< Command line interface handler.
   character(:), allocatable                   :: args_     !< Problem options, local variable.
   integer(I_P)                                :: error     !< Error trapping flag.

   args_ = '' ; if (present(args)) args_ = ' '//args
   select case(problem)
   case('lcce')
      allocate(integrand_lcce :: integrand)
   case('oscillation')
      allocate(integrand_oscillation :: integrand)
   case('polynomial')
      allocate(integrand_polynomial :: integrand)
   case('riccati')
      allocate(integrand_riccati :: integrand)
   case default
      error stop 'make_problem: unknown problem "'//problem//'"'
   endselect
   call cli%init(progname='foodie_test_suite', description='FOODIE tests suite')
   call integrand%set_cli(cli=cli)
   call cli%parse(args=problem//args_, error=error)
   if (error /= 0) error stop 'make_problem: wrong options "'//problem//args_//'"'
   call integrand%parse_cli(cli=cli)
   endfunction make_problem

   function has_fast_mode(scheme)
   !< Return .true. if the scheme provides the fast mode integrate.
   character(*), intent(in)              :: scheme        !< Scheme name.
   logical                               :: has_fast_mode !< Inquire result.
   class(integrator_object), allocatable :: integrator    !< Integrator.

   call foodie_integrator_factory(scheme=scheme, integrator=integrator, stages=STAGES)
   has_fast_mode = integrator%has_fast_mode()
   endfunction has_fast_mode

   pure function is_linear_multistep(scheme)
   !< Return .true. if the scheme is a linear multistep one: its polynomial exactness degree is exactly its order.
   !<
   !< For the other schemes the quadrature defined by the weights/nodes can be more accurate than the scheme order.
   character(*), intent(in) :: scheme              !< Scheme name.
   logical                  :: is_linear_multistep !< Inquire result.

   is_linear_multistep = (index(scheme, 'adams_')  == 1 .or. &
                          index(scheme, 'back_df') == 1 .or. &
                          index(scheme, 'lmm_ssp') == 1 .or. &
                          trim(scheme) == 'leapfrog')
   endfunction is_linear_multistep

   pure function polynomial_degree(scheme) result(degree)
   !< Return the maximum degree of the polynomials (solution of a quadrature) integrated exactly by the scheme.
   character(*), intent(in) :: scheme !< Scheme name.
   integer(I_P)             :: degree !< Degree.

   select case(trim(scheme))
   case('leapfrog_raw')
      degree = 1
   case('runge_kutta_lssp_stages_s_order_s')
      degree = min(2_I_P, STAGES)
   case('runge_kutta_lssp_stages_s_order_s_1')
      degree = min(2_I_P, STAGES - 1)
   case default
      degree = convergence_order(scheme)
   endselect
   endfunction polynomial_degree

   pure function convergence_order(scheme) result(order)
   !< Return the expected order of accuracy of the scheme on its convergence test problem.
   character(*), intent(in) :: scheme !< Scheme name.
   integer(I_P)             :: order  !< Order of accuracy.

   select case(trim(scheme))
   case('euler_explicit')
      order = 1
   case('leapfrog')
      order = 2
   case('leapfrog_raw')
      order = 1
   case('runge_kutta_emd_stages_7_order_4')
      order = 5
   case('runge_kutta_lssp_stages_s_order_s')
      order = STAGES
   case('runge_kutta_lssp_stages_s_order_s_1')
      order = STAGES - 1
   case default
      if (index(scheme, '_order_') > 0) then
         order = trailing_integer(scheme)
      elseif (index(scheme, 'adams_moulton_') == 1) then
         order = trailing_integer(scheme) + 1
      else ! adams_bashforth_N, adams_bashforth_moulton_N, back_df_N
         order = trailing_integer(scheme)
      endif
   endselect
   contains
      pure function trailing_integer(string) result(n)
      !< Return the integer at the end of the string (after the last underscore).
      character(*), intent(in) :: string !< String.
      integer(I_P)             :: n      !< Integer.

      read(string(index(string, '_', back=.true.)+1:len_trim(string)), *) n
      endfunction trailing_integer
   endfunction convergence_order

   subroutine convergence_problem(scheme, problem, args, final_time, Dt0)
   !< Return the convergence test problem suitable for the scheme.
   !<
   !< The nonlinear non autonomous Riccati problem is the default; the leapfrog family is tested on the oscillation problem
   !< (leapfrog being unstable on dissipative problems) and the linear SSP family on the linear autonomous lcce problem.
   character(*),              intent(in)  :: scheme     !< Scheme name.
   character(:), allocatable, intent(out) :: problem    !< Problem name.
   character(:), allocatable, intent(out) :: args       !< Problem options.
   real(R_P),                 intent(out) :: final_time !< Final integration time.
   real(R_P),                 intent(out) :: Dt0        !< Coarsest time step.

   if (index(scheme, 'leapfrog') == 1) then
      problem = 'oscillation' ; args = '' ; final_time = 6e4_R_P ; Dt0 = 3e3_R_P
   elseif (index(scheme, 'runge_kutta_lssp') == 1) then
      problem = 'lcce' ; args = '' ; final_time = 2._R_P ; Dt0 = 0.5_R_P
   else
      problem = 'riccati' ; args = '' ; final_time = 2._R_P ; Dt0 = 0.5_R_P
   endif
   endsubroutine convergence_problem
endmodule foodie_test_suite_utils
