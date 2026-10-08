!< Integration driver shared by the CLI tester and the tests suite.

module foodie_test_integrate
!< Integration driver shared by the CLI tester and the tests suite.
!<
!< Multi-step schemes are started up with the exact solution, thus the measured error is due only to the scheme.

use foodie, only : foodie_integrator_factory,              &
                   integrator_multistage_object,           &
                   integrator_multistage_multistep_object, &
                   integrator_multistep_object,            &
                   integrator_object
use foodie_test_integrand_tester_object, only : integrand_tester_object
use penf, only : I_P, R_P

implicit none
private
public :: exporter_object
public :: integrate
public :: observed_order

type, abstract :: exporter_object
   !< Exporter of the integrand solution during the integration (hook for the CLI tester).
   contains
      procedure(export_interface), pass(self), deferred :: export !< Export the current integrand solution.
endtype exporter_object

abstract interface
   subroutine export_interface(self, integrand, integrator, step, time, is_last)
   !< Export the current integrand solution.
   import :: exporter_object, integrand_tester_object, integrator_object, I_P, R_P
   class(exporter_object),         intent(inout) :: self    !< Exporter.
   class(integrand_tester_object), intent(in) :: integrand  !< Integrand.
   class(integrator_object),       intent(in) :: integrator !< Integrator.
   integer(I_P),                   intent(in) :: step       !< Time step counter, 0 before the first step.
   real(R_P),                      intent(in) :: time       !< Time.
   logical,                        intent(in) :: is_last    !< Flag for the last export.
   endsubroutine export_interface
endinterface

contains
   subroutine integrate(scheme, integrand_0, Dt, final_time, iterations, stages, is_fast, error, exporter, error_max)
   !< Integrate integrand by means of the given scheme, returning the error at the final time.
   character(*),                   intent(in)           :: scheme      !< Selected scheme.
   class(integrand_tester_object), intent(in)           :: integrand_0 !< Initial conditions.
   real(R_P),                      intent(in)           :: Dt          !< Time step.
   real(R_P),                      intent(in)           :: final_time  !< Final integration time.
   integer(I_P),                   intent(in)           :: iterations  !< Number of fixed point iterations.
   integer(I_P),                   intent(in)           :: stages      !< Number of stages.
   logical,                        intent(in)           :: is_fast     !< Activate fast mode integration.
   real(R_P),                      intent(out)          :: error(:)    !< Error of integrand integration.
   class(exporter_object),         intent(inout), optional :: exporter !< Exporter.
   real(R_P),                      intent(out),   optional :: error_max(:) !< Max error over all time steps (L-inf in time).
   class(integrand_tester_object), allocatable          :: integrand   !< Integrand.
   class(integrator_object),       allocatable          :: integrator  !< The integrator.
   real(R_P)                                            :: time        !< Time.
   integer(I_P)                                         :: step        !< Time steps counter.

   allocate(integrand, mold=integrand_0) ; integrand = integrand_0

   call foodie_integrator_factory(scheme=scheme, integrator=integrator, stages=stages, &
                                  tolerance=1e2_R_P, iterations=iterations, autoupdate=.true., U=integrand_0)

   step = 0
   time = 0._R_P
   if (present(error_max)) error_max = 0._R_P
   if (present(exporter)) &
      call exporter%export(integrand=integrand, integrator=integrator, step=step, time=time, is_last=.false.)

   select type(integrator)
   class is(integrator_multistage_object)
      do
         step = step + 1
         if (is_fast) then
            call integrator%integrate_fast(U=integrand, Dt=Dt, t=time)
         else
            call integrator%integrate(U=integrand, Dt=Dt, t=time)
         endif
         time = time + Dt
         if (present(error_max)) error_max = max(error_max, integrand%error(t=time, U0=integrand_0))
         if ((time >= final_time)) exit
         if (present(exporter)) &
            call exporter%export(integrand=integrand, integrator=integrator, step=step, time=time, is_last=.false.)
      enddo

   class is(integrator_multistep_object)
      do
         step = step + 1
         if (integrator%steps_number() >= step) then
            time = time + Dt
            integrator%Dt(step) = Dt
            integrator%t(step) = time
            integrator%previous(step) = integrand%exact_solution(t=time, U0=integrand_0)
            integrand = integrator%previous(step)
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand, Dt=Dt, t=time)
            else
               call integrator%integrate(U=integrand, Dt=Dt, t=time)
            endif
            time = time + Dt
         endif
         if (present(error_max)) error_max = max(error_max, integrand%error(t=time, U0=integrand_0))
         if ((time >= final_time)) exit
         if (present(exporter)) &
            call exporter%export(integrand=integrand, integrator=integrator, step=step, time=time, is_last=.false.)
      enddo

   class is(integrator_multistage_multistep_object)
      do
         step = step + 1
         if (integrator%steps_number() >= step) then
            time = time + Dt
            integrator%Dt(step) = Dt
            integrator%t(step) = time
            integrator%previous(step) = integrand%exact_solution(t=time, U0=integrand_0)
            integrand = integrator%previous(step)
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand, Dt=Dt, t=time)
            else
               call integrator%integrate(U=integrand, Dt=Dt, t=time)
            endif
            time = time + Dt
         endif
         if (present(error_max)) error_max = max(error_max, integrand%error(t=time, U0=integrand_0))
         if ((time >= final_time)) exit
         if (present(exporter)) &
            call exporter%export(integrand=integrand, integrator=integrator, step=step, time=time, is_last=.false.)
      enddo
   endselect

   if (present(exporter)) &

      call exporter%export(integrand=integrand, integrator=integrator, step=step, time=time, is_last=.true.)
   error = integrand%error(t=time, U0=integrand_0)
   endsubroutine integrate

   pure function observed_order(error, Dt)
   !< Estimate the order of accuracy using 2 subsequent refined numerical solutions.
   real(R_P), intent(in) :: error(1:, 1:)                        !< Computed errors.
   real(R_P), intent(in) :: Dt(1:)                               !< Time steps used.
   real(R_P)             :: observed_order(1:size(error, dim=1)) !< Estimation of the order of accuracy.
   integer(I_P)          :: v                                    !< Variables counter.

   do v=1, size(error, dim=1)
     observed_order(v) = log(error(v, 1) / error(v, 2)) / log(Dt(1) / Dt(2))
   enddo
   endfunction observed_order
endmodule foodie_test_integrate
