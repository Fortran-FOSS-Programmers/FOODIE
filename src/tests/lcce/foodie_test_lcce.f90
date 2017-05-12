!< Test FOODIE with the integration of linear constant coefficients equation.

module foodie_test_lcce_test
!< Linear constant coefficients test handler definition.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use flap, only : command_line_interface
use foodie, only : foodie_integrator_class_names,      &
                   foodie_integrator_factory,          &
                   foodie_integrator_schemes,          &
                   integrator_adams_bashforth,         &
                   integrator_adams_bashforth_moulton, &
                   integrator_adams_moulton,           &
                   integrator_back_df,                 &
                   integrator_euler_explicit,          &
                   integrator_leapfrog,                &
                   integrator_lmm_ssp,                 &
                   integrator_lmm_ssp_vss,             &
                   integrator_ms_runge_kutta_ssp,      &
                   integrator_object,                  &
                   integrator_runge_kutta_emd,         &
                   integrator_runge_kutta_ls,          &
                   integrator_runge_kutta_lssp,        &
                   integrator_runge_kutta_ssp,         &
                   is_available, is_class_available
use foodie_test_integrand_lcce, only : integrand_lcce
use penf, only : I_P, R_P, FR_P, str, strz

implicit none
private
public :: lcce_test

type :: lcce_test
   !< Class to handle linear constant coefficients equation test(s).
   !<
   !< Test is driven by the Command Line Interface (CLI) options.
   !<
   !< Test has only 1 public method `execute`: it executes test(s) accordingly to cli options.
   private
   type(command_line_interface) :: cli                     !< Command line interface handler.
   integer(I_P)                 :: error=0                 !< Error handler.
   character(99)                :: scheme=''               !< Scheme used.
   logical                      :: is_fast=.false.         !< Flag for activating fast schemes.
   integer(I_P)                 :: implicit_iterations=0   !< Number of iterations (implicit solvers).
   integer(I_P)                 :: stages=0                !< Number of stages.
   real(R_P), allocatable       :: Dt(:)                   !< Time step(s) exercised.
   real(R_P)                    :: a=0._R_P                !< *a* coefficient.
   real(R_P)                    :: b=0._R_P                !< *b* coefficient.
   real(R_P)                    :: U0=0._R_P               !< Initial conditions.
   real(R_P)                    :: final_time=0.0_R_P      !< Final integration time.
   logical                      :: results=.false.         !< Flag for activating results saving.
   character(99)                :: output=''               !< Output files basename.
   logical                      :: exact_solution=.false.  !< Flag for activating exact solution saving.
   logical                      :: errors_analysis=.false. !< Flag for activating errors analysis.
   contains
      ! public methods
      procedure, pass(self) :: execute !< Execute selected test(s).
      ! private methods
      procedure, pass(self), private :: initialize !< Initialize test: set Command Line Interface, parse it and check its validity.
      procedure, pass(self), private :: test       !< Perform the test.
endtype lcce_test

contains
   ! public methods
   subroutine execute(self)
   !< Execute selected test(s).
   class(lcce_test), intent(inout) :: self                  !< Test.
   character(99), allocatable      :: integrator_schemes(:) !< Name of FOODIE integrator schemes.
   integer(I_P)                    :: s                     !< Counter.

   call self%initialize
   if (trim(adjustl(self%scheme))/='all') then
      if (is_class_available(scheme=self%scheme)) then
        integrator_schemes = foodie_integrator_schemes(class_name=self%scheme)
      elseif (is_available(scheme=self%scheme)) then
        integrator_schemes = [trim(adjustl(self%scheme))]
      endif
   else
     integrator_schemes = foodie_integrator_schemes()
   endif
   do s=1, size(integrator_schemes, dim=1)
     self%scheme = trim(integrator_schemes(s))
     call self%test(scheme=self%scheme)
   enddo
   endsubroutine execute

   ! private methods
   subroutine initialize(self)
   !< Initialize test: set Command Line Interface, parse it and check its validity.
   class(lcce_test), intent(inout) :: self !< Test.

   call set_cli
   call parse_cli
   contains
      subroutine set_cli()
      !< Set Command Line Interface.

      associate(cli => self%cli)
        call cli%init(progname    = 'foodie_test_lcce',                                                         &
                      authors     = 'Fortran-FOSS-Programmers',                                                 &
                      license     = 'GNU GPLv3',                                                                &
                      description = 'Test FOODIE library on linear constant coefficients equation integration', &
                      examples    = ["foodie_test_lcce --scheme euler_explicit --save_results  ",               &
                                     "foodie_test_lcce --scheme all -r                         "])
        call cli%add(switch='--scheme', switch_ab='-s', help='integrator scheme used', required=.false., def='all', act='store')
        call cli%add(switch='--fast', help='activate fast solvers', required=.false., act='store_true', def='.false.')
        call cli%add(switch='--iterations', help='iterations number for implicit schemes', required=.false., act='store', def='5')
        call cli%add(switch='--stages', help='stages number', required=.false., def='4', act='store')
        call cli%add(switch='--time_step', switch_ab='-Dt', nargs='+', help='time step', required=.false., def='0.01', act='store')
        call cli%add(switch='--a', switch_ab='-a', help='"a" coefficient', required=.false., def='-1.0', act='store')
        call cli%add(switch='--b', switch_ab='-b', help='"b" coefficient', required=.false., def='0.0',  act='store')
        call cli%add(switch='--U0', switch_ab='-U0', help='initial state', required=.false., def='1.0',  act='store')
        call cli%add(switch='--t_final', switch_ab='-tf', help='final integration time', required=.false., def='1.0', act='store')
        call cli%add(switch='--save_results', switch_ab='-r', help='save result', required=.false., act='store_true', def='.false.')
        call cli%add(switch='--output', help='output files basename', required=.false., act='store', def='foodie_test_lcce')
        call cli%add(switch='--exact_solution', help='save exact solutiion', required=.false., act='store_true', def='.false.')
        call cli%add(switch='--errors_analysis', help='peform errors analysis', required=.false., act='store_true', def='.false.')
      endassociate
      endsubroutine set_cli

      subroutine parse_cli()
      !< Parse Command Line Interface and check its validity.
      character(99), allocatable :: integrator_class_names(:) !< Name of FOODIE integrator classes.
      character(99), allocatable :: integrator_schemes(:)     !< Name of FOODIE integrator schemes.
      integer(I_P)               :: i                         !< Counter.

      call self%cli%parse(error=self%error)
      call self%cli%get(switch='-s', val=self%scheme, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--fast', val=self%is_fast, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--iterations', val=self%implicit_iterations, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--stages', val=self%stages, error=self%error) ; if (self%error/=0) stop
      call self%cli%get_varying(switch='-Dt', val=self%Dt, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='-a', val=self%a, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='-b', val=self%b, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='-U0', val=self%U0, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='-tf', val=self%final_time, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='-r', val=self%results, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--output', val=self%output, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--exact_solution', val=self%exact_solution, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--errors_analysis', val=self%errors_analysis, error=self%error) ; if (self%error/=0) stop

      if (trim(adjustl(self%scheme)) /= 'all') then
         if (.not.is_available(scheme=self%scheme)) then
            integrator_class_names = foodie_integrator_class_names()
            integrator_schemes     = foodie_integrator_schemes()
            print '(A)', 'Error: the scheme "'//trim(adjustl(self%scheme))//'" is unknown!'
            print '(A)', 'Supported classes of schemes are:'
            do i=1, size(integrator_class_names, dim=1)
               print '(A)', '  '//trim(integrator_class_names(i))
            enddo
            print '(A)', 'Supported schemes are:'
            do i=1, size(integrator_schemes, dim=1)
               print '(A)', '  '//trim(integrator_schemes(i))
            enddo
            stop
         endif
      endif

      if (.not.is_dt_valid()) then
         print "(A)", 'Error: the final integration time must be an exact multiple of the time step used!'
         print "(A)", 'Final integration time: '//str(self%final_time, .true.)
         print "(A)", 'Time step: '//str(self%Dt, .true.)
         stop
      endif

      if (self%errors_analysis) then
         if (.not.(size(self%Dt, dim=1) > 1)) self%errors_analysis = .false.
      endif
      endsubroutine parse_cli

      function is_dt_valid()
      !< Verify if the selected time step Dt is valid.
      logical      :: is_dt_valid !< Return true is the selected time step Dt is valid.
      integer(I_P) :: t           !< Counter.

      is_dt_valid = .true.
      do t=1, size(self%Dt)
         is_dt_valid = ((self%final_time - int(self%final_time/self%Dt(t), I_P)*self%Dt(t))==0)
         if (.not.is_dt_valid) exit
      enddo
      endfunction is_dt_valid
   endsubroutine initialize

   subroutine test(self, scheme)
   !< Perform the test.
   class(lcce_test), intent(in) :: self          !< Test.
   character(*),     intent(in) :: scheme        !< Selected scheme.
   real(R_P), allocatable       :: solution(:,:) !< Solution at each time step.
   real(R_P), allocatable       :: error(:)      !< Error (norm L2) with respect the exact solution.
   real(R_P), allocatable       :: order(:)      !< Observed order based on subsequent refined solutions.
   integer(I_P)                 :: last_step     !< Last time step computed.
   integer(I_P)                 :: t             !< Counter.

   print "(A)", trim(adjustl(scheme))

   allocate(error(1:size(self%Dt)))
   error = 0.0_R_P
   if (self%errors_analysis) then
      allocate(order(1:size(error, dim=1)-1))
      order = 0.0_R_P
   endif

   do t=1, size(self%Dt, dim=1)
      call integrate(scheme=scheme,                       &
                     a=self%a,                            &
                     b=self%b,                            &
                     U0=self%U0,                          &
                     final_time=self%final_time,          &
                     Dt=self%Dt(t),                       &
                     iterations=self%implicit_iterations, &
                     stages=self%stages,                  &
                     is_fast=self%is_fast,                &
                     solution=solution,                   &
                     error=error(t),                      &
                     last_step=last_step)
      if (allocated(solution)) then
         print "(A,G0,A,G0,A,G0)", "  steps: ", last_step, " Dt: ", self%Dt(t), ", error: ", error(t)
         if (self%errors_analysis.and.t>1) then
            order(t-1) = observed_order(error=error(t-1:t), Dt=self%Dt(t-1:t))
            print "(A,F7.3)", "    observed order: ", order(t-1)
         endif
         call save_results(results=self%results,                    &
                           output=self%output,                      &
                           scheme=trim(adjustl(scheme)),            &
                           a=self%a,                                &
                           b=self%b,                                &
                           U0=self%U0,                              &
                           save_exact_solution=self%exact_solution, &
                           solution=solution(0:1,0:last_step))
      endif
   enddo
   endsubroutine test

   ! non type bound procedures
   subroutine check_scheme_has_fast_mode(scheme, integrator)
   !< Check if a scheme support fast mode integrate.
   character(*),             intent(in) :: scheme     !< Scheme name.
   class(integrator_object), intent(in) :: integrator !< Integrator instance.

   if (.not.integrator%has_fast_mode()) then
      write(stderr, '(A)') 'error: '//trim(adjustl(scheme))//' has not fast integrate mode!'
      stop
   endif
   endsubroutine check_scheme_has_fast_mode

   subroutine integrate(scheme, a, b, U0, final_time, Dt, iterations, stages, is_fast, solution, error, last_step)
   !< Integrate domain by means of the given scheme.
   character(*),           intent(in)    :: scheme        !< Selected scheme.
   real(R_P),              intent(in)    :: a             !< *a* coefficient.
   real(R_P),              intent(in)    :: b             !< *b* coefficient.
   real(R_P),              intent(in)    :: U0            !< Initial state.
   real(R_P),              intent(in)    :: final_time    !< Final integration time.
   real(R_P),              intent(in)    :: Dt            !< Time step.
   integer(I_P),           intent(in)    :: iterations    !< Number of fixed point iterations.
   integer(I_P),           intent(in)    :: stages        !< Number of stages.
   logical,                intent(in)    :: is_fast       !< Activate fast mode integration.
   real(R_P), allocatable, intent(out)   :: solution(:,:) !< Solution at each time step, X-Y.
   real(R_P),              intent(out)   :: error         !< Error (norm L2) with respect the exact solution.
   integer(I_P),           intent(out)   :: last_step     !< Last time step computed.
   class(integrator_object), allocatable :: integrator    !< The integrator.
   type(integrand_lcce)                  :: domain        !< Linear constant coefficients equation field.
   type(integrand_lcce), allocatable     :: previous(:)   !< Previous time steps solutions.
   type(integrand_lcce), allocatable     :: stage(:)      !< Runge-Kutta stages.
   type(integrand_lcce)                  :: buffer        !< Buffer field.
   type(integrand_lcce)                  :: filter        !< Filter displacement.
   integer                               :: step          !< Time steps counter.
   integer                               :: step_offset   !< Time steps counter offset for slicing previous data array.
   real(R_P), allocatable                :: Dts(:)        !< Time steps for variable stepsize.
   real(R_P)                             :: Dt_a          !< Adaptive time step.

   call domain%initialize(a=a, b=b, U0=U0)

   if (allocated(solution)) deallocate(solution) ; allocate(solution(0:domain%integrand_dimension(), 0:int(final_time/Dt)))
   solution = 0.0_R_P
   solution(1:, 0) = domain%output()

   call foodie_integrator_factory(scheme=scheme, integrator=integrator, stages=stages, tolerance=1e2_R_P)
   if (is_fast) call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integrator=integrator)

   if (integrator%is_multistep()) allocate(previous(1:integrator%steps_number()))
   if (integrator%is_multistage()) allocate(stage(1:integrator%stages_number()))

   select type(integrator)
   type is(integrator_adams_bashforth)
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (integrator%steps_number() >= step) then
           domain = [domain%exact_solution(t=step * Dt)]
           previous(step) = domain
        else
           if (is_fast) then
              call integrator%integrate_fast(U=domain,          &
                                             previous=previous, &
                                             buffer=buffer,     &
                                             Dt=Dt,             &
                                             t=solution(0, step-integrator%steps_number():step-1))
           else
              call integrator%integrate(U=domain,          &
                                        previous=previous, &
                                        Dt=Dt,             &
                                        t=solution(0, step-integrator%steps_number():step-1))
           endif
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_adams_bashforth_moulton)
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (integrator%steps_number() >= step) then
           domain = [domain%exact_solution(t=step * Dt)]
           previous(step) = domain
        else
          if (is_fast) then
             call integrator%integrate_fast(U=domain,          &
                                            previous=previous, &
                                            buffer=buffer,     &
                                            Dt=Dt,             &
                                            t=solution(0, step-integrator%steps_number():step-1))
          else
             call integrator%integrate(U=domain,          &
                                       previous=previous, &
                                       Dt=Dt,             &
                                       t=solution(0, step-integrator%steps_number():step-1))
          endif
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_adams_moulton)
      if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps_number()+1))
      if (integrator%steps_number()==0) then
        step_offset = 1                         ! for 0 step-(a convention)-solver offset is 1
      else
        step_offset = integrator%steps_number() ! for >0 step-solver offset is steps
      endif
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (integrator%steps_number() >= step) then
           domain = [domain%exact_solution(t=step * Dt)]
           previous(step) = domain
        else
          if (is_fast) then
             if (iterations>1) then
               call integrator%integrate_fast(U=domain,                              &
                                              previous=previous,                     &
                                              buffer=buffer,                         &
                                              Dt=Dt,                                 &
                                              t=solution(0,step-step_offset:step-1), &
                                              iterations=iterations)
             else
               call integrator%integrate_fast(U=domain,          &
                                              previous=previous, &
                                              buffer=buffer,     &
                                              Dt=Dt,             &
                                              t=solution(0,step-step_offset:step-1))
             endif
          else
             if (iterations>1) then
               call integrator%integrate(U=domain,                              &
                                         previous=previous,                     &
                                         Dt=Dt,                                 &
                                         t=solution(0,step-step_offset:step-1), &
                                         iterations=iterations)
             else
               call integrator%integrate(U=domain,          &
                                         previous=previous, &
                                         Dt=Dt,             &
                                         t=solution(0,step-step_offset:step-1))
             endif
          endif
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_back_df)
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (integrator%steps_number() >= step) then
           domain = [domain%exact_solution(t=step * Dt)]
           previous(step) = domain
        else
          if (is_fast) then
             call integrator%integrate_fast(U=domain,          &
                                            previous=previous, &
                                            buffer=buffer,     &
                                            Dt=Dt,             &
                                            t=solution(0, step-integrator%steps_number():step-1))
          else
             call integrator%integrate(U=domain,          &
                                       previous=previous, &
                                       Dt=Dt,             &
                                       t=solution(0, step-integrator%steps_number():step-1))
          endif
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_euler_explicit)
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (is_fast) then
           call integrator%integrate_fast(U=domain, buffer=buffer, Dt=Dt, t=solution(0, step))
        else
           call integrator%integrate(U=domain, Dt=Dt, t=solution(0, step))
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step
   type is(integrator_leapfrog)
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (integrator%steps_number() >= step) then
           domain = [domain%exact_solution(t=step * Dt)]
           previous(step) = domain
        else
           if (index(scheme, 'raw') > 0 ) then
              if (is_fast) then
                call integrator%integrate_fast(U=domain, previous=previous, buffer=buffer, Dt=Dt, t=solution(0, step),filter=filter)
              else
                call integrator%integrate(U=domain, previous=previous, Dt=Dt, t=solution(0, step), filter=filter)
              endif
           else
              if (is_fast) then
                call integrator%integrate_fast(U=domain, previous=previous, buffer=buffer, Dt=Dt, t=solution(0, step))
              else
                call integrator%integrate(U=domain, previous=previous, Dt=Dt, t=solution(0, step))
              endif
           endif
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_lmm_ssp)
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (integrator%steps_number() >= step) then
           domain = [domain%exact_solution(t=step * Dt)]
           previous(step) = domain
        else
          if (is_fast) then
             call integrator%integrate_fast(U=domain,          &
                                            previous=previous, &
                                            buffer=buffer,     &
                                            Dt=Dt,             &
                                            t=solution(0, step-integrator%steps_number():step-1))
          else
             call integrator%integrate(U=domain,          &
                                       previous=previous, &
                                       Dt=Dt,             &
                                       t=solution(0, step-integrator%steps_number():step-1))
          endif
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_lmm_ssp_vss)
      if (allocated(Dts)) deallocate(Dts) ; allocate(Dts(1:integrator%steps_number())) ; Dts = Dt
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (integrator%steps_number() >= step) then
           domain = [domain%exact_solution(t=step * Dt)]
           previous(step) = domain
        else
          if (is_fast) then
             call integrator%integrate_fast(U=domain,          &
                                            previous=previous, &
                                            buffer=buffer,     &
                                            Dt=Dts,            &
                                            t=solution(0, step-integrator%steps_number():step-1))
          else
             call integrator%integrate(U=domain,          &
                                       previous=previous, &
                                       Dt=Dts,            &
                                       t=solution(0, step-integrator%steps_number():step-1))
          endif
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_ms_runge_kutta_ssp)
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (integrator%steps_number() >= step) then
           domain = [domain%exact_solution(t=step * Dt)]
           previous(step) = domain
        else
          if (is_fast) then
             call integrator%integrate_fast(U=domain,          &
                                            previous=previous, &
                                            stage=stage,       &
                                            buffer=buffer,     &
                                            Dt=Dt,             &
                                            t=solution(0, step-integrator%steps_number():step-1))
          else
             call integrator%integrate(U=domain,          &
                                       previous=previous, &
                                       stage=stage,       &
                                       Dt=Dt,             &
                                       t=solution(0, step-integrator%steps_number():step-1))
          endif
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_runge_kutta_emd)
       step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        Dt_a = Dt
        if (is_fast) then
           call integrator%integrate_fast(U=domain, stage=stage, buffer=buffer, Dt=Dt_a, t=solution(0, step))
        else
           call integrator%integrate(U=domain, stage=stage, Dt=Dt_a, t=solution(0, step))
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_runge_kutta_ls)
      if (allocated(stage)) deallocate(stage) ; allocate(stage(1:integrator%registers_number()))
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (is_fast) then
           call integrator%integrate_fast(U=domain, stage=stage, buffer=buffer, Dt=Dt, t=solution(0, step))
        else
           call integrator%integrate(U=domain, stage=stage, Dt=Dt, t=solution(0, step))
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_runge_kutta_lssp)
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (is_fast) then
           call integrator%integrate_fast(U=domain, stage=stage, buffer=buffer, Dt=Dt, t=solution(0, step))
        else
           call integrator%integrate(U=domain, stage=stage, Dt=Dt, t=solution(0, step))
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step

   type is(integrator_runge_kutta_ssp)
      step = 0
      do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
        step = step + 1
        if (is_fast) then
           call integrator%integrate_fast(U=domain, stage=stage, buffer=buffer, Dt=Dt, t=solution(0, step))
        else
           call integrator%integrate(U=domain, stage=stage, Dt=Dt, t=solution(0, step))
        endif
        solution(0, step) = step * Dt
        solution(1:, step) = domain%output()
      enddo
      last_step = step
   endselect

   error = 0._R_P
   do step=0, last_step
      domain = solution(1:, step)
      error = error + (domain%output() - domain%exact_solution(t=solution(0, step))) ** 2
   enddo
   error = sqrt(error)
   endsubroutine integrate

   subroutine save_results(results, output, scheme, a, b, U0, save_exact_solution, solution)
   !< Save results (and plots).
   logical,      intent(in)      :: results             !< Flag for activating results saving.
   character(*), intent(in)      :: output              !< Output files basename coming from CLI.
   character(*), intent(in)      :: scheme              !< Selected scheme: must be defined into *solvers*.
   real(R_P),    intent(in)      :: a                   !< *a* coefficient.
   real(R_P),    intent(in)      :: b                   !< *b* coefficient.
   real(R_P),    intent(in)      :: U0                  !< Initial state.
   logical,      intent(in)      :: save_exact_solution !< Flag for saving exact solution.
   real(R_P),    intent(in)      :: solution(0:,0:)     !< Solution at each time step.
   character(len=:), allocatable :: title               !< Output files title.
   character(len=:), allocatable :: basename            !< Output files basename.
   type(integrand_lcce)          :: lcce                !< Linear constant coefficients field.
   integer(I_P)                  :: rawfile             !< Raw file unit for saving results.
   integer(I_P)                  :: s                   !< Counter.

   basename = trim(adjustl(output))//'-'//trim(strz(ubound(solution, dim=2), 10))//'-time_steps-'//trim(adjustl(scheme))
   title = 'linear constant coefficients equation integration, solver='//trim(adjustl(scheme))
   if (results) then
      open(newunit=rawfile, file=basename//'.dat')
      write(rawfile, '(A)')'TITLE="'//title//'"'
      write(rawfile, '(A)')'VARIABLES="t" "u"'
      write(rawfile, '(A)')'ZONE T="'//trim(adjustl(scheme))//'"'
      do s=0, ubound(solution, dim=2)
         write(rawfile, '(2('//FR_P//',1X))')solution(:, s)
      enddo
      close(rawfile)
   endif
   if (save_exact_solution) then
      call lcce%initialize(a=a, b=b, U0=U0)
      basename = trim(adjustl(output))//'-'//trim(strz(ubound(solution, dim=2), 10))//'-time_steps-exact_solution'
      title = 'linear constant coefficients equation integration, solver=exact solution'
      open(newunit=rawfile, file=basename//'.dat')
      write(rawfile, '(A)')'TITLE="'//title//'"'
      write(rawfile, '(A)')'VARIABLES="t" "u"'
      write(rawfile, '(A)')'ZONE T="exact solution"'
      do s=0, ubound(solution, dim=2)
         write(rawfile, '(2('//FR_P//',1X))')solution(0, s), lcce%exact_solution(t=solution(0, s))
      enddo
      close(rawfile)
   endif
   endsubroutine save_results

   pure function observed_order(error, Dt)
   !< Estimate the order of accuracy using 2 subsequent refined numerical solutions.
   real(R_P), intent(IN) :: error(1:2)     !< Computed errors.
   real(R_P), intent(IN) :: Dt(1:2)        !< Time steps used.
   real(R_P)             :: observed_order !< Estimation of the order of accuracy.

   observed_order = log(error(1) / error(2)) / log(Dt(1) / Dt(2))
   endfunction observed_order
endmodule foodie_test_lcce_test

program foodie_test_lcce
!< Test FOODIE with the integration of linear constant coefficients equation.

use foodie_test_lcce_test, only : lcce_test

implicit none
type(lcce_test) :: test !< Linear constant coefficients equation test.

call test%execute
endprogram foodie_test_lcce
