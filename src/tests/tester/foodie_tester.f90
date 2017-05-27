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
use foodie_test_integrand_oscillation, only : integrand_oscillation
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
   character(99)                               :: test                !< Test executed.
   character(99)                               :: scheme              !< Scheme used.
   real(R_P), allocatable                      :: Dt(:)               !< Time step(s) exercised.
   logical                                     :: is_fast             !< Flag for activating fast schemes.
   integer(I_P)                                :: implicit_iterations !< Number of iterations (implicit solvers).
   integer(I_P)                                :: stages              !< Number of stages.
   real(R_P)                                   :: final_time          !< Final integration time.
   logical                                     :: save_results        !< Flag for activating results saving.
   character(99)                               :: output              !< Output files basename.
   integer(I_P)                                :: save_frequency      !< Save frequency.
   type(integrand_ladvection)                  :: ladvection_0        !< Initial conditions for linear advection test.
   type(integrand_oscillation)                 :: oscillation_0       !< Initial conditions for oscillation test.
   class(integrand_tester_object), allocatable :: integrand_0         !< Initial conditions.
   contains
      ! public methods
      procedure, pass(self) :: execute !< Execute selected test(s).
      ! private methods
      procedure, pass(self), private :: initialize !< Initialize test: set Command Line Interface, parse it and check its validity.
endtype test_object

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
   allocate(error(1:self%integrand_0%integrand_dimension(), 1:size(self%Dt, dim=1)))
   if (size(self%Dt, dim=1) > 1) allocate(order(1:self%integrand_0%integrand_dimension(), 1:size(error, dim=2)-1))
   do s=1, size(integrator_schemes, dim=1)
      print '(A)', trim(integrator_schemes(s))
      do t=1, size(self%Dt)
         call integrate(scheme=trim(integrator_schemes(s)),          &
                        integrand_0=self%integrand_0,                &
                        Dt=self%Dt(t),                               &
                        final_time=self%final_time,                  &
                        iterations=self%implicit_iterations,         &
                        stages=self%stages,                          &
                        is_fast=self%is_fast,                        &
                        save_results=self%save_results,              &
                        output_base_name=trim(adjustl(self%output)), &
                        save_frequency=self%save_frequency,          &
                        error=error(:,t))
         print*, 'Dt = ', self%Dt(t), ', error (max) = ', maxval(error(:,t))
         if (t > 1) then
            order(:, t-1) = observed_order(error=error(:, t-1:t), Dt=self%Dt(t-1:t))
            print '(A,F6.2)', ' Observed order (min) =', minval(order(:,t-1))
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
                       examples    = ["foodie_tester --scheme euler_explicit --save_results  ", &
                                      "foodie_tester --scheme all -r                         "])
         call cli%add(switch='--test', switch_ab='-t', help='test executed', required=.false., def='oscillation', &
                      act='store', choices='linear_advection,oscillation')
         call cli%add(switch='--scheme', switch_ab='-s', help='integrator scheme used', required=.false., def='all', act='store')
         call cli%add(switch='--time_step', switch_ab='-Dt', nargs='+', help='time step', required=.false., def='1e2', act='store')
         call cli%add(switch='--fast', help='activate fast solvers', required=.false., act='store_true', def='.false.')
         call cli%add(switch='--iterations', help='iterations number for implicit schemes', required=.false., act='store', def='5')
         call cli%add(switch='--stages', help='stages number', required=.false., def='2', act='store')
         call cli%add(switch='--final_time', switch_ab='-ft', help='integration time', required=.false., def='1e6', act='store')
         call cli%add(switch='--save_results', switch_ab='-r',help='save result', required=.false., act='store_true', def='.false.')
         call cli%add(switch='--output', help='output file basename', required=.false., act='store', def='foodie_test')
         call cli%add(switch='--save_frequency', help='save frequency', required=.false., act='store', def='1')
      endassociate
      call self%ladvection_0%set_cli(cli=self%cli)
      call self%oscillation_0%set_cli(cli=self%cli)
      endsubroutine set_cli

      subroutine parse_cli()
      !< Parse Command Line Interface and check its validity.
      character(99), allocatable :: integrator_class_names(:) !< Name of FOODIE integrator classes.
      character(99), allocatable :: integrator_schemes(:)     !< Name of FOODIE integrator schemes.
      real(R_P), allocatable     :: initial_state(:)          !< Initial integrand state.
      integer(I_P)               :: i                         !< Counter.

      call self%cli%parse(error=self%error)
      call self%cli%get(switch='-t', val=self%test, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='-s', val=self%scheme, error=self%error) ; if (self%error/=0) stop
      call self%cli%get_varying(switch='-Dt', val=self%Dt, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--fast', val=self%is_fast, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--iterations', val=self%implicit_iterations, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--stages', val=self%stages, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='-ft', val=self%final_time, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='-r', val=self%save_results, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--output', val=self%output, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--save_frequency', val=self%save_frequency, error=self%error) ; if (self%error/=0) stop
      call self%ladvection_0%parse_cli(cli=self%cli)
      call self%oscillation_0%parse_cli(cli=self%cli)

      select case(trim(adjustl(self%test)))
      case('linear_advection')
         allocate(integrand_ladvection :: self%integrand_0)
         self%integrand_0 = self%ladvection_0
      case('oscillation')
         allocate(integrand_oscillation :: self%integrand_0)
         self%integrand_0 = self%oscillation_0
      endselect

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
         if (is_dt_valid) then
            select type(integrand=>self%integrand_0)
            type is(integrand_ladvection)
               is_dt_valid = is_dt_valid .and. self%Dt(t) <= integrand%Dt(final_time=self%final_time)
               if (.not.is_dt_valid) then
                  write(stderr, '(A)') 'error: Dt violates CFL condition, Dt_max = '//str(integrand%Dt(final_time=self%final_time))
               endif
            endselect
         endif
         if (.not.is_dt_valid) exit
      enddo
      endfunction is_dt_valid
   endsubroutine initialize

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

   subroutine integrate(scheme, integrand_0, Dt, final_time, iterations, stages, is_fast, save_results, output_base_name, &
                        save_frequency, error)
   !< Integrate integrand by means of the given scheme.
   character(*),                   intent(in)   :: scheme           !< Selected scheme.
   class(integrand_tester_object), intent(in)   :: integrand_0      !< Initial conditions.
   real(R_P),                      intent(in)   :: Dt               !< Time step.
   real(R_P),                      intent(in)   :: final_time       !< Final integration time.
   integer(I_P),                   intent(in)   :: iterations       !< Number of fixed point iterations.
   integer(I_P),                   intent(in)   :: stages           !< Number of stages.
   logical,                        intent(in)   :: is_fast          !< Activate fast mode integration.
   logical,                        intent(in)   :: save_results     !< Save results.
   character(*),                   intent(in)   :: output_base_name !< Base name of output results file.
   integer(I_P),                   intent(in)   :: save_frequency   !< Save frequency.
   real(R_P),                      intent(out)  :: error(:)         !< Error of integrand integration.
   class(integrand_tester_object) , allocatable :: integrand        !< Integrand.
   class(integrator_object), allocatable        :: integrator       !< The integrator.
   type(integrator_runge_kutta_ssp)             :: integrator_start !< The (auto) start integrator.
   real(R_P)                                    :: time             !< Time.
   integer(I_P)                                 :: step             !< Time steps counter.

   allocate(integrand, mold=integrand_0) ; integrand = integrand_0

   call foodie_integrator_factory(scheme=scheme, integrator=integrator, stages=stages, &
                                  tolerance=1e2_R_P, iterations=iterations, autoupdate=.true., U=integrand_0)
   if (is_fast) call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integrator=integrator)

   if (integrator%is_multistep()) call integrator_start%initialize(scheme='runge_kutta_ssp_stages_5_order_4', U=integrand_0)

   step = 0
   time = 0._R_P
   if (save_results) call integrand%export_tecplot(file_name=output_base_name//                                     &
                                                             integrand%description(prefix='-')//                    &
                                                             integrator%description(prefix='-')//                   &
                                                             '-steps_'//trim(strz(int(final_time/Dt), 10))//'.dat', &
                                                             t=time, scheme=scheme)

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
         if ((time >= final_time)) exit
         call integrand_export_tecplot
      enddo

   class is(integrator_multistep_object)
      do
         step = step + 1
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, Dt=Dt, t=time)
            integrator%previous(step) = integrand
            time = time + Dt
            integrator%Dt(step) = Dt
            integrator%t(step) = time
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand, Dt=Dt, t=time)
            else
               call integrator%integrate(U=integrand, Dt=Dt, t=time)
            endif
            time = time + Dt
         endif
         if ((time >= final_time)) exit
         call integrand_export_tecplot
      enddo

   class is(integrator_multistage_multistep_object)
      do
         step = step + 1
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, Dt=Dt, t=time)
            integrator%previous(step) = integrand
            time = time + Dt
            integrator%Dt(step) = Dt
            integrator%t(step) = time
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand, Dt=Dt,t=time)
            else
               call integrator%integrate(U=integrand, Dt=Dt, t=time)
            endif
            time = time + Dt
         endif
         if ((time >= final_time)) exit
         call integrand_export_tecplot
      enddo
   endselect

   call integrand_close_tecplot

   call integrand_compute_error
   contains
      subroutine integrand_close_tecplot
      !< Close current integrand tecplot file.

      select type(integrand)
      type is(integrand_ladvection)
         if (save_results .and. mod(step, save_frequency)==0) call integrand%export_tecplot(t=time, scheme=scheme)
      type is(integrand_oscillation)
         if (save_results .and. mod(step, save_frequency)==0) call integrand%export_tecplot(t=time)
      endselect
      if (save_results) call integrand%export_tecplot(close_file=.true.)
      endsubroutine integrand_close_tecplot

      subroutine integrand_compute_error
      !< Close current integrand tecplot file.

      select type(integrand)
      type is(integrand_ladvection)
      type is(integrand_oscillation)
         error = abs(integrand%output() - integrand%exact_solution(t=time))
      endselect
      endsubroutine integrand_compute_error

      subroutine integrand_export_tecplot
      !< Export current integrand solution to tecplot file.

      select type(integrand)
      type is(integrand_ladvection)
         if (save_results .and. mod(step, save_frequency)==0) call integrand%export_tecplot(t=time, scheme=scheme)
      type is(integrand_oscillation)
         if (save_results .and. mod(step, save_frequency)==0) call integrand%export_tecplot(t=time)
      endselect
      endsubroutine integrand_export_tecplot
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
endmodule foodie_test_object

program foodie_tester
!< Tester factory of FOODIE integrators.

use foodie_test_object, only : test_object

implicit none
type(test_object) :: test !< FOODIE test.

call test%execute
endprogram foodie_tester
