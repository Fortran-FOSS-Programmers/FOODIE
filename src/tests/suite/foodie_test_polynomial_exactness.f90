!< Test the polynomial exactness of all FOODIE schemes.

program foodie_test_polynomial_exactness
!< Test the polynomial exactness of all FOODIE schemes.
!<
!< The quadrature \(U_t = q t^{q-1}\) has the polynomial solution \(t^q\): a scheme of order *p* must integrate it exactly for
!< \(q \le p\) (order conditions of multistep schemes, quadrature ones of multistage schemes, correct use of the stages
!< times); a linear multistep scheme must not for \(q = p+1\) (its error constant being not null), that verifies the claimed
!< order is not underestimated. Being the residual independent of *U* errors are not propagated (no stability constraint),
!< thus also the highest order schemes are verified. Both standard and fast integrate modes are tested.

use foodie, only : foodie_integrator_schemes
use foodie_test_integrand_tester_object, only : integrand_tester_object
use foodie_test_integrate, only : integrate
use foodie_test_suite_utils, only : has_fast_mode, is_linear_multistep, make_problem, polynomial_degree, STAGES
use penf, only : I_P, R_P, str

implicit none
real(R_P),    parameter                     :: DT          = 0.1_R_P   !< Time step.
real(R_P),    parameter                     :: FINAL_TIME  = 2._R_P    !< Final time (20 steps).
real(R_P),    parameter                     :: TOL_EXACT   = 1e-12_R_P !< Max relative error of exactly integrated polynomials.
real(R_P),    parameter                     :: TOL_INEXACT = 1e-10_R_P !< Min relative error of not exactly integrated ones.
integer(I_P), parameter                     :: ITERATIONS  = 20_I_P    !< Fixed point iterations of implicit schemes.
character(99), allocatable                  :: schemes(:)              !< Schemes.
class(integrand_tester_object), allocatable :: integrand_0             !< Initial conditions.
real(R_P)                                   :: error(1)                !< Error.
real(R_P)                                   :: rel_error               !< Relative error.
real(R_P)                                   :: worst_exact             !< Worst relative error of exact degrees.
integer(I_P)                                :: p                       !< Expected exactness degree.
integer(I_P)                                :: q                       !< Polynomial degree.
integer(I_P)                                :: s                       !< Schemes counter.
integer(I_P)                                :: m                       !< Modes counter.
integer(I_P)                                :: failures                !< Failures counter.
logical                                     :: is_passed               !< Test result.
character(:), allocatable                   :: reason                  !< Failure reason.

failures = 0
schemes = foodie_integrator_schemes()
do s=1, size(schemes, dim=1)
   p = polynomial_degree(trim(schemes(s)))
   do m=0, 1
      if (m == 1) then
         if (.not.has_fast_mode(trim(schemes(s)))) cycle
      endif
      is_passed = .true. ; reason = '' ; worst_exact = 0._R_P
      do q=0, p + merge(1, 0, is_linear_multistep(trim(schemes(s))))
         if (allocated(integrand_0)) deallocate(integrand_0)
         allocate(integrand_0, source=make_problem(problem='polynomial', &
                                                args='-q '//trim(str(q, no_sign=.true.))))
         call integrate(scheme=trim(schemes(s)), integrand_0=integrand_0, Dt=DT, final_time=FINAL_TIME, &
                        iterations=ITERATIONS, stages=STAGES, is_fast=(m == 1), error=error)
         rel_error = error(1) / max(1._R_P, FINAL_TIME ** q)
         if (q <= p) then
            worst_exact = max(worst_exact, rel_error)
            if (.not.(rel_error <= TOL_EXACT)) then
               is_passed = .false. ; reason = reason//' degree '//trim(str(q, no_sign=.true.))//' not exact ('// &
                                                       trim(str(rel_error, compact=.true.))//')'
            endif
         else
            if (.not.(rel_error >= TOL_INEXACT)) then
               is_passed = .false. ; reason = reason//' degree '//trim(str(q, no_sign=.true.))//' unexpectedly exact ('// &
                                                       trim(str(rel_error, compact=.true.))//'), order underestimated'
            endif
         endif
      enddo
      if (.not.is_passed) failures = failures + 1
      print '(A)', merge('PASS ', 'FAIL ', is_passed)//trim(schemes(s))//merge(' [fast]', '       ', m == 1)// &
                   ' exact up to degree '//trim(str(p, no_sign=.true.))//', max rel error '// &
                   trim(str(worst_exact, compact=.true.))//reason
   enddo
enddo
print '(A)', trim(str(failures, no_sign=.true.))//' failures'
if (failures > 0) stop 1
endprogram foodie_test_polynomial_exactness
