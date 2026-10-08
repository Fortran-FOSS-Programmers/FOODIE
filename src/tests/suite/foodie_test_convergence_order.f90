!< Test the observed order of accuracy of all FOODIE schemes.

program foodie_test_convergence_order
!< Test the observed order of accuracy of all FOODIE schemes.
!<
!< Each scheme integrates its convergence problem (the nonlinear non autonomous Riccati one by default) with time steps halved
!< from the coarsest one: the observed order is computed on the finest couple of solutions whose errors are both in the
!< asymptotic range \([E_{min}, E_{max}]\), above round-off. The error is the max over all time steps (L-inf norm in time):
!< the error at the final time alone is not reliable, the global error function can vanish there. The test fails if the observed order is lower than the expected
!< one minus a tolerance, or if no couple is in range. Orders higher than 6 are not verified (only reported): in double precision
!< the asymptotic range is not reached before round-off, their order conditions being verified by the polynomial exactness
!< test (multistep schemes) and by the audit of the tableaus (scripts/audit_rk_tableaus.py, Runge-Kutta schemes). Both standard and fast integrate modes are tested.

use foodie, only : foodie_integrator_schemes
use foodie_test_integrand_tester_object, only : integrand_tester_object
use foodie_test_integrate, only : integrate
use foodie_test_suite_utils, only : convergence_order, convergence_problem, has_fast_mode, make_problem, STAGES
use penf, only : I_P, R_P, str

implicit none
integer(I_P), parameter                     :: REFINEMENTS     = 10_I_P    !< Number of time step halvings.
integer(I_P), parameter                     :: ITERATIONS      = 20_I_P    !< Fixed point iterations of implicit schemes.
integer(I_P), parameter                     :: MAX_ORDER_CHECK = 6_I_P     !< Max order verifiable in double precision.
real(R_P),    parameter                     :: E_MIN           = 1e-12_R_P !< Min error of the asymptotic range.
real(R_P),    parameter                     :: E_MAX           = 1e-3_R_P  !< Max error of the asymptotic range.
real(R_P),    parameter                     :: ORDER_TOL       = 0.3_R_P   !< Tolerance on the observed order.
character(99), allocatable                  :: schemes(:)                  !< Schemes.
class(integrand_tester_object), allocatable :: integrand_0                 !< Initial conditions.
character(:), allocatable                   :: problem                     !< Convergence problem.
character(:), allocatable                   :: args                        !< Convergence problem options.
real(R_P),    allocatable                   :: error_(:)                   !< Max error over time steps of all variables.
real(R_P),    allocatable                   :: error_final(:)              !< Error at the final time of all variables.
real(R_P)                                   :: error(0:REFINEMENTS)        !< Error (max over variables).
real(R_P)                                   :: final_time                  !< Final time.
real(R_P)                                   :: Dt0                         !< Coarsest time step.
real(R_P)                                   :: order                       !< Observed order.
integer(I_P)                                :: p                           !< Expected order.
integer(I_P)                                :: k                           !< Refinements counter.
integer(I_P)                                :: k_obs                       !< Refinement of the observed order.
integer(I_P)                                :: s                           !< Schemes counter.
integer(I_P)                                :: m                           !< Modes counter.
integer(I_P)                                :: failures                    !< Failures counter.
character(5)                                :: outcome                     !< Test outcome.

failures = 0
schemes = foodie_integrator_schemes()
do s=1, size(schemes, dim=1)
   p = convergence_order(trim(schemes(s)))
   call convergence_problem(scheme=trim(schemes(s)), problem=problem, args=args, final_time=final_time, Dt0=Dt0)
   if (allocated(integrand_0)) deallocate(integrand_0)
   allocate(integrand_0, source=make_problem(problem=problem, args=args)) ! defined assignment needs an allocated lhs
   if (allocated(error_)) deallocate(error_)
   allocate(error_(1:size(integrand_0%error(t=0._R_P), dim=1)))
   error_final = error_
   do m=0, 1
      if (m == 1) then
         if (.not.has_fast_mode(trim(schemes(s)))) cycle
      endif
      do k=0, REFINEMENTS
         call integrate(scheme=trim(schemes(s)), integrand_0=integrand_0, Dt=Dt0 / 2._R_P ** k, final_time=final_time, &
                        iterations=ITERATIONS, stages=STAGES, is_fast=(m == 1), error=error_final, error_max=error_)
         error(k) = maxval(error_)
      enddo
      k_obs = -1
      do k=REFINEMENTS - 1, 0, -1
         if (in_range(error(k)) .and. in_range(error(k+1))) then
            k_obs = k
            exit
         endif
      enddo
      order = 0._R_P
      if (k_obs >= 0) order = log(error(k_obs) / error(k_obs+1)) / log(2._R_P)
      if (p > MAX_ORDER_CHECK) then
         outcome = 'SKIP ' ! observed order reported, verification left to the polynomial exactness test
      elseif (k_obs >= 0) then
         outcome = merge('PASS ', 'FAIL ', order >= p - ORDER_TOL)
      else
         outcome = 'FAIL '
      endif
      if (outcome == 'FAIL ') failures = failures + 1
      print '(A)', outcome//trim(schemes(s))//merge(' [fast]', '       ', m == 1)//' on '//problem//                  &
                   ': expected order '//trim(str(p, no_sign=.true.))//', observed '//trim(str(order, compact=.true.))// &
                   ' (errors: '//errors_string()//')'
   enddo
enddo
print '(A)', trim(str(failures, no_sign=.true.))//' failures'
if (failures > 0) stop 1
contains
   elemental function in_range(e)
   !< Return .true. if the error is in the asymptotic range.
   real(R_P), intent(in) :: e        !< Error.
   logical               :: in_range !< Inquire result.

   in_range = (e >= E_MIN .and. e <= E_MAX)
   endfunction in_range

   function errors_string() result(string)
   !< Return the errors as a compact string.
   character(:), allocatable :: string !< Errors string.
   integer(I_P)              :: k      !< Refinements counter.

   string = ''
   do k=0, REFINEMENTS
      string = string//trim(str(error(k), compact=.true.))
      if (k < REFINEMENTS) string = string//' '
   enddo
   endfunction errors_string
endprogram foodie_test_convergence_order
