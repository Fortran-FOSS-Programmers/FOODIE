!< Tester factory of FOODIE integrators.

module foodie_test_object
!< Definition of [[test_object]] for FOODIE tester factory.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use flap, only : command_line_interface
use foodie, only : foodie_integrator_class_names,          &
                   foodie_integrator_factory,              &
                   foodie_integrator_schemes,              &
                   integrand_object,                       &
                   integrator_adams_bashforth,             &
                   integrator_adams_bashforth_moulton,     &
                   integrator_adams_moulton,               &
                   integrator_back_df,                     &
                   integrator_euler_explicit,              &
                   integrator_leapfrog,                    &
                   integrator_lmm_ssp,                     &
                   integrator_lmm_ssp_vss,                 &
                   integrator_ms_runge_kutta_ssp,          &
                   integrator_multistage_object,           &
                   integrator_multistage_multistep_object, &
                   integrator_multistep_object,            &
                   integrator_object,                      &
                   integrator_runge_kutta_emd,             &
                   integrator_runge_kutta_ls,              &
                   integrator_runge_kutta_lssp,            &
                   integrator_runge_kutta_ssp,             &
                   is_available, is_class_available
use foodie_test_integrand_ladvection, only : integrand_ladvection
use foodie_test_integrand_lcce, only : integrand_lcce
use foodie_test_integrand_oscillation, only : integrand_oscillation
use foodie_test_integrand_polynomial, only : integrand_polynomial
use foodie_test_integrand_riccati, only : integrand_riccati
use foodie_test_integrate, only : exporter_object, integrate, observed_order
use foodie_test_integrand_tester_object, only : integrand_tester_object
use penf, only : I_P, R_P, FR_P, str, strz

implicit none
private
public :: test_object

type :: test_object
   !< Generic FOODIE test object.
   !<
   !< Test is driven by the Command Line Interface (CLI) options.
   !<
   !< Test has only 1 public method `execute`: it executes test(s) accordingly to cli options.
   private
   type(command_line_interface)                :: cli                 !< Command line interface handler.
   integer(I_P)                                :: error               !< Error handler.
   character(99)                               :: scheme              !< Scheme used.
   real(R_P), allocatable                      :: Dt(:)               !< Time step(s) exercised.
   logical                                     :: is_fast             !< Flag for activating fast schemes.
   integer(I_P)                                :: implicit_iterations !< Number of iterations (implicit solvers).
   integer(I_P)                                :: stages              !< Number of stages.
   real(R_P)                                   :: final_time          !< Final integration time.
   logical                                     :: save_results        !< Flag for activating results saving.
   character(99)                               :: output              !< Output files basename.
   integer(I_P)                                :: save_frequency      !< Save frequency.
   type(integrand_lcce)                        :: lcce_0              !< Initial conditions for linear constant coefficients eq.
   type(integrand_ladvection)                  :: ladvection_0        !< Initial conditions for linear advection test.
   type(integrand_oscillation)                 :: oscillation_0       !< Initial conditions for oscillation test.
   type(integrand_polynomial)                  :: polynomial_0        !< Initial conditions for polynomial test.
   type(integrand_riccati)                     :: riccati_0           !< Initial conditions for Riccati test.
   class(integrand_tester_object), allocatable :: integrand_0         !< Initial conditions.
   contains
      ! public methods
      procedure, pass(self) :: execute !< Execute selected test(s).
      ! private methods
      procedure, pass(self), private :: initialize !< Initialize test: set Command Line Interface, parse it and check its validity.
endtype test_object

type, extends(exporter_object) :: tecplot_exporter
   !< Exporter of results into tecplot files, checking also the fast mode availability before the first step.
   private
   character(99)                               :: output_base_name !< Output files basename.
   character(99)                               :: scheme           !< Scheme used.
   real(R_P)                                   :: Dt               !< Time step.
   real(R_P)                                   :: final_time       !< Final integration time.
   integer(I_P)                                :: save_frequency   !< Save frequency.
   logical                                     :: save_results     !< Flag for activating results saving.
   logical                                     :: is_fast          !< Flag for activating fast schemes.
   class(integrand_tester_object), allocatable :: integrand_0      !< Initial conditions.
   contains
      procedure, pass(self) :: export => exporter_export         !< Export the current integrand solution.
      procedure, pass(self) :: initialize => exporter_initialize !< Initialize the exporter for a (scheme, Dt) run.
endtype tecplot_exporter

contains
   ! public methods
   subroutine execute(self)
   !< Execute test(s).
   class(test_object), intent(inout) :: self                  !< Test.
   character(99), allocatable        :: integrator_schemes(:) !< Name of FOODIE integrator schemes.
   real(R_P), allocatable            :: error(:,:)            !< Error of integrand integration.
   real(R_P), allocatable            :: order(:,:)            !< Observed order of accuravy.
   integer(I_P)                      :: s                     !< Counter.
   integer(I_P)                      :: t                     !< Counter.
   type(tecplot_exporter)            :: exporter              !< Exporter of results into tecplot files.

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
   allocate(error(1:size(self%integrand_0%error(t=0._R_P), dim=1), 1:size(self%Dt, dim=1)))
   if (size(self%Dt, dim=1) > 1) allocate(order(1:size(error, dim=1), 1:size(error, dim=2)-1))
   do s=1, size(integrator_schemes, dim=1)
      print '(A)', trim(integrator_schemes(s))
      do t=1, size(self%Dt)
         call self%integrand_0%initialize(Dt=self%Dt(t))
         call exporter%initialize(test=self, scheme=trim(integrator_schemes(s)), Dt=self%Dt(t))
         print '(A)', self%integrand_0%description()
         call integrate(scheme=trim(integrator_schemes(s)),  &
                        integrand_0=self%integrand_0,        &
                        Dt=self%Dt(t),                       &
                        final_time=self%final_time,          &
                        iterations=self%implicit_iterations, &
                        stages=self%stages,                  &
                        is_fast=self%is_fast,                &
                        error=error(:,t),                    &
                        exporter=exporter)
         print*, 'Dt = ', self%Dt(t), ', error = ', error(:,t)
         if (t > 1) then
            order(:, t-1) = observed_order(error=error(:, t-1:t), Dt=self%Dt(t-1:t))
            print '(A,'//trim(str(size(order, dim=1), no_sign=.true.))//'F6.2)', ' Observed order =', order(:,t-1)
         endif
      enddo
   enddo
   endsubroutine execute

   ! private methods
   subroutine initialize(self)
   !< Initialize test: set Command Line Interface, parse it and check its validity.
   class(test_object), intent(inout) :: self !< Test.

   call set_cli
   call parse_cli
   contains
      subroutine set_cli()
      !< Set Command Line Interface.

      associate(cli => self%cli)
         call cli%init(progname    = 'foodie_tester',                                           &
                       authors     = 'Fortran-FOSS-Programmers',                                &
                       license     = 'GNU GPLv3',                                               &
                       description = 'Tester factory of FOODIE integrators',                    &
                       examples    = ["foodie_tester test --scheme euler_explicit --save_results  ", &
                                      "foodie_tester test --scheme all -r                         "])
         call cli%add_group(description='general test settings', group='test')
         call cli%add(group='test', switch='--scheme', switch_ab='-s', help='integrator scheme used', required=.false., def='all', &
                      act='store')
         call cli%add(group='test', switch='--time_step', switch_ab='-Dt', nargs='+', help='time step', required=.false., &
                      def='1e2', act='store')
         call cli%add(group='test', switch='--fast', help='activate fast solvers', required=.false., act='store_true', &
                      def='.false.')
         call cli%add(group='test', switch='--iterations', help='iterations number for implicit schemes', required=.false., &
                      act='store', def='5')
         call cli%add(group='test', switch='--stages', help='stages number', required=.false., def='2', act='store')
         call cli%add(group='test', switch='--final_time', switch_ab='-ft', help='integration time', required=.false., def='1', &
                      act='store')
         call cli%add(group='test', switch='--save_results', switch_ab='-r',help='save result', required=.false., &
                      act='store_true', def='.false.')
         call cli%add(group='test', switch='--output', help='output file basename', required=.false., act='store', &
                      def='foodie_test')
         call cli%add(group='test', switch='--save_frequency', help='save frequency', required=.false., act='store', &
                      def='1')
      endassociate
      call self%lcce_0%set_cli(cli=self%cli)
      call self%ladvection_0%set_cli(cli=self%cli)
      call self%oscillation_0%set_cli(cli=self%cli)
      call self%polynomial_0%set_cli(cli=self%cli)
      call self%riccati_0%set_cli(cli=self%cli)
      endsubroutine set_cli

      subroutine parse_cli()
      !< Parse Command Line Interface and check its validity.
      character(99), allocatable :: integrator_class_names(:) !< Name of FOODIE integrator classes.
      character(99), allocatable :: integrator_schemes(:)     !< Name of FOODIE integrator schemes.
      real(R_P), allocatable     :: initial_state(:)          !< Initial integrand state.
      integer(I_P)               :: i                         !< Counter.

      call self%cli%parse(error=self%error)
      call self%cli%get(group='test', switch='-s', val=self%scheme, error=self%error) ; if (self%error/=0) stop
      call self%cli%get_varying(group='test', switch='-Dt', val=self%Dt, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='test', switch='--fast', val=self%is_fast, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='test', switch='--iterations',val=self%implicit_iterations,error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='test', switch='--stages', val=self%stages, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='test', switch='-ft', val=self%final_time, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='test', switch='-r', val=self%save_results, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='test', switch='--output', val=self%output, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='test', switch='--save_frequency',val=self%save_frequency,error=self%error) ; if (self%error/=0) stop
      ! parse only the options of the invoked problem, the other problems being not used
      if (self%cli%run_command('lcce'))             call self%lcce_0%parse_cli(cli=self%cli)
      if (self%cli%run_command('linear_advection')) call self%ladvection_0%parse_cli(cli=self%cli)
      if (self%cli%run_command('oscillation'))      call self%oscillation_0%parse_cli(cli=self%cli)
      if (self%cli%run_command('polynomial'))       call self%polynomial_0%parse_cli(cli=self%cli)
      if (self%cli%run_command('riccati'))          call self%riccati_0%parse_cli(cli=self%cli)

      if     (self%cli%run_command('lcce')) then
         allocate(integrand_lcce :: self%integrand_0)
         select type(integrand_0=>self%integrand_0)
         type is(integrand_lcce)
            integrand_0 = self%lcce_0
         endselect
      elseif (self%cli%run_command('linear_advection')) then
         allocate(integrand_ladvection :: self%integrand_0)
         select type(integrand_0=>self%integrand_0)
         type is(integrand_ladvection)
            integrand_0 = self%ladvection_0
         endselect
      elseif (self%cli%run_command('oscillation')) then
         allocate(integrand_oscillation :: self%integrand_0)
         select type(integrand_0=>self%integrand_0)
         type is(integrand_oscillation)
            integrand_0 = self%oscillation_0
         endselect
      elseif (self%cli%run_command('polynomial')) then
         allocate(integrand_polynomial :: self%integrand_0)
         select type(integrand_0=>self%integrand_0)
         type is(integrand_polynomial)
            integrand_0 = self%polynomial_0
         endselect
      elseif (self%cli%run_command('riccati')) then
         allocate(integrand_riccati :: self%integrand_0)
         select type(integrand_0=>self%integrand_0)
         type is(integrand_riccati)
            integrand_0 = self%riccati_0
         endselect
      else
         write(stderr, '(A)') 'error: select the problem to integrate, one of "lcce", "linear_advection", "oscillation", '// &
                              '"polynomial", "riccati"'
         stop 1
      endif

      if (.not.is_dt_valid()) then
         write(stderr, '(A)') 'error: the final integration time must be an exact multiple of the time step used'
         write(stderr, '(A)') '       and time step must respect CFL condition, if any'
         write(stderr, '(A)') 'Final integration time: '//str(self%final_time, .true.)
         write(stderr, '(A)') 'Time step: '//str(self%Dt, .true.)
         stop
      endif

      if (trim(adjustl(self%scheme)) /= 'all') then
         if (.not.is_available(scheme=self%scheme)) then
            integrator_class_names = foodie_integrator_class_names()
            integrator_schemes     = foodie_integrator_schemes()
            write(stderr, '(A)') 'error: the scheme "'//trim(adjustl(self%scheme))//'" is unknown!'
            write(stderr, '(A)') 'Supported classes of schemes are:'
            do i=1, size(integrator_class_names, dim=1)
               write(stderr, '(A)') '  '//trim(integrator_class_names(i))
            enddo
            write(stderr, '(A)') 'Supported schemes are:'
            do i=1, size(integrator_schemes, dim=1)
               write(stderr, '(A)') '  '//trim(integrator_schemes(i))
            enddo
            stop
         endif
      endif
      endsubroutine parse_cli

      function is_dt_valid()
      !< Verify if the selected time step Dt is valid.
      logical      :: is_dt_valid !< Return true is the selected time step Dt is valid.
      integer(I_P) :: t           !< Counter.

      is_dt_valid = .true.
      do t=1, size(self%Dt)
         is_dt_valid = ((self%final_time - int(self%final_time/self%Dt(t), I_P)*self%Dt(t))==0)
         ! if (is_dt_valid) then
         !    select type(integrand=>self%integrand_0)
         !    type is(integrand_ladvection)
         !       is_dt_valid = is_dt_valid .and. self%Dt(t) <= integrand%Dt(final_time=self%final_time)
         !       if (.not.is_dt_valid) then
         !          write(stderr, '(A)') 'error: Dt violates CFL condition, Dt_max = '//str(integrand%Dt(final_time=self%final_time))
         !       endif
         !    endselect
         ! endif
         if (.not.is_dt_valid) exit
      enddo
      endfunction is_dt_valid
   endsubroutine initialize

   ! tecplot_exporter methods
   subroutine exporter_initialize(self, test, scheme, Dt)
   !< Initialize the exporter for a (scheme, Dt) run.
   class(tecplot_exporter), intent(inout) :: self   !< Exporter.
   type(test_object),       intent(in)    :: test   !< Test.
   character(*),            intent(in)    :: scheme !< Scheme used.
   real(R_P),               intent(in)    :: Dt     !< Time step.

   self%output_base_name = test%output
   self%scheme = scheme
   self%Dt = Dt
   self%final_time = test%final_time
   self%save_frequency = test%save_frequency
   self%save_results = test%save_results
   self%is_fast = test%is_fast
   if (allocated(self%integrand_0)) deallocate(self%integrand_0)
   allocate(self%integrand_0, source=test%integrand_0)
   endsubroutine exporter_initialize

   subroutine exporter_export(self, integrand, integrator, step, time, is_last)
   !< Export current integrand solution to tecplot file, check fast mode availability before the first step.
   class(tecplot_exporter),        intent(inout) :: self       !< Exporter.
   class(integrand_tester_object), intent(in)    :: integrand  !< Integrand.
   class(integrator_object),       intent(in)    :: integrator !< Integrator.
   integer(I_P),                   intent(in)    :: step       !< Time step counter, 0 before the first step.
   real(R_P),                      intent(in)    :: time       !< Time.
   logical,                        intent(in)    :: is_last    !< Flag for the last export.

   if (step == 0 .and. self%is_fast) call check_scheme_has_fast_mode(scheme=trim(self%scheme), integrator=integrator)
   if (.not.self%save_results) return
   if (step == 0) then
      call integrand%export_tecplot(file_name=trim(adjustl(self%output_base_name))//                           &
                                              integrand%description(prefix='-')//                              &
                                              integrator%description(prefix='-')//                             &
                                              '-steps_'//trim(strz(int(self%final_time/self%Dt), 10))//'.dat', &
                                    t=time,                                                                    &
                                    scheme=trim(self%scheme),                                                  &
                                    with_exact_solution=.true.,                                                &
                                    U0=self%integrand_0)
      return
   endif
   if (mod(step, self%save_frequency)==0) then
      select type(integrand)
      type is(integrand_ladvection)
         call integrand%export_tecplot(t=time, scheme=trim(self%scheme))
      type is(integrand_oscillation)
         if (is_last) then
            call integrand%export_tecplot(t=time)
         else
            call integrand%export_tecplot(t=time, with_exact_solution=.true.)
         endif
      class default
         if (.not.is_last) call integrand%export_tecplot(t=time, with_exact_solution=.true.)
      endselect
   endif
   if (is_last) call integrand%export_tecplot(close_file=.true.)
   endsubroutine exporter_export

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
endmodule foodie_test_object

program foodie_tester
!< Tester factory of FOODIE integrators.

use foodie_test_object, only : test_object

implicit none
type(test_object) :: test !< FOODIE test.

call test%execute
endprogram foodie_tester
