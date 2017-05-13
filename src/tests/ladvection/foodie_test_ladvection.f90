!< Test FOODIE with the integration of 1D linear advection PDE.

module foodie_test_ladvection_test
!< Oscillation test handler definition.

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
use foodie_test_integrand_ladvection, only : integrand_ladvection
use penf, only : I_P, R_P, FR_P, str, strz

implicit none
private
public :: ladvection_test

type :: ladvection_test
   !< Class to handle 1D linear advection test(s).
   !<
   !< Test is driven by the Command Line Interface (CLI) options.
   !<
   !< Test has only 1 public method `execute`: it executes test(s) accordingly to cli options.
   private
   type(command_line_interface) :: cli                    !< Command line interface handler.
   integer(I_P)                 :: error=0                !< Error handler.
   character(99)                :: scheme=''              !< Scheme used.
   real(R_P)                    :: CFL                    !< CFL value.
   logical                      :: is_fast=.false.        !< Flag for activating fast schemes.
   integer(I_P)                 :: implicit_iterations=0  !< Number of iterations (implicit solvers).
   integer(I_P)                 :: stages=0               !< Number of stages.
   integer(I_P)                 :: steps_max              !< Maximum number of time steps.
   real(R_P)                    :: final_time=0.0_R_P     !< Final integration time.
   character(99)                :: w_scheme=''            !< WENO Scheme used.
   integer(I_P)                 :: weno_order=0           !< WENO reconstruction order.
   real(R_P)                    :: weno_eps=0._R_P        !< WENO epsilon value.
   real(R_P)                    :: a=0._R_P               !< Advection coefficient.
   integer(I_P)                 :: Ni=0                   !< Number of grid cells.
   character(3)                 :: BC_L                   !< Left boundary condition type.
   character(3)                 :: BC_R                   !< Right boundary condition type.
   logical                      :: save_results=.false.   !< Flag for activating results saving.
   character(99)                :: output=''              !< Output files basename.
   logical                      :: verbose=.false.        !< Flag for activating verbose output.
   type(integrand_ladvection)   :: integrand_0            !< Initial conditions.
   contains
      ! public methods
      procedure, pass(self) :: execute !< Execute selected test(s).
      ! private methods
      procedure, pass(self), private :: initialize !< Initialize test: set Command Line Interface, parse it and check its validity.
endtype ladvection_test

contains
   ! public methods
   subroutine execute(self)
   !< Execute test(s).
   class(ladvection_test), intent(inout) :: self                  !< Test.
   character(99), allocatable            :: integrator_schemes(:) !< Name of FOODIE integrator schemes.
   integer(I_P)                          :: s                     !< Counter.

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
      call integrate(scheme=trim(integrator_schemes(s)),  &
                     CFL=self%CFL,                        &
                     integrand_0=integrand_0,             &
                     final_step=self%final_step,          &
                     final_time=self%final_time,          &
                     iterations=self%implicit_iterations, &
                     stages=self%stages,                  &
                     is_fast=self%is_fast,                &
                     save_results=self%save_results,      &
                     save_frequency=self%save_frequency,  &
                     output_base_name=self%output)
   enddo
   endsubroutine execute

   ! private methods
   subroutine initialize(self)
   !< Initialize test: set Command Line Interface, parse it and check its validity.
   class(ladvection_test), intent(inout) :: self !< Test.

   call set_cli
   call parse_cli
   contains
      subroutine set_cli()
      !< Set Command Line Interface.

      associate(cli => self%cli)
         call cli%init(progname    = 'foodie_test_ladvection',                                           &
                       authors     = 'Fortran-FOSS-Programmers',                                         &
                       license     = 'GNU GPLv3',                                                        &
                       description = 'Test FOODIE library on 1D linear advection PDE integration',       &
                       examples    = ["foodie_test_ladvection --scheme euler_explicit --save_results  ", &
                                      "foodie_test_ladvection --scheme all -r                         "])
         call cli%add(switch='--scheme', switch_ab='-s', help='integrator scheme used', required=.false., def='all', act='store')
         call cli%add(switch='--cfl', help='CFL value', required=.false., act='store', def='0.8')
         call cli%add(switch='--fast', help='activate fast solvers', required=.false., act='store_true', def='.false.')
         call cli%add(switch='--iterations', help='iterations number for implicit schemes', required=.false., act='store', def='5')
         call cli%add(switch='--stages', help='stages number', required=.false., def='2', act='store')
         call cli%add(switch='--steps', help='number time steps performed', required=.false., act='store', def='100')
         call cli%add(switch='--t_final', switch_ab='-tf', help='final integration time', required=.false., def='0', act='store')
         call cli%add(switch='--w-scheme', help='WENO scheme', required=.false., act='store', def='reconstructor-JS', &
           choices='reconstructor-JS,reconstructor-M-JS,reconstructor-M-Z,reconstructor-Z')
         call cli%add(switch='--weno-order', help='WENO order', required=.false., act='store', def='1')
         call cli%add(switch='--weno-eps', help='WENO epsilon parameter', required=.false., act='store', def='0.000001')
         call cli%add(switch='-a', help='advection coefficient', required=.false., act='store', def='1.0')
         call cli%add(switch='--Ni', help='number finite volumes used', required=.false., act='store', def='100')
         call cli%add(switch='--save_results', switch_ab='-r',help='save result', required=.false., act='store_true', def='.false.')
         call cli%add(switch='--output', help='output file basename', required=.false., act='store', def='foodie_test_ladvection')
         call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.')
      endassociate
      endsubroutine set_cli

      subroutine parse_cli()
      !< Parse Command Line Interface and check its validity.
      character(99), allocatable :: integrator_class_names(:) !< Name of FOODIE integrator classes.
      character(99), allocatable :: integrator_schemes(:)     !< Name of FOODIE integrator schemes.
      real(R_P), allocatable     :: initial_state(:)          !< Initial integrand state.
      integer(I_P)               :: i                         !< Counter.

      call self%cli%parse(error=self%error)
      call self%cli%get(switch='-s', val=self%scheme, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--cfl', val=self%CFL, error=error) ; if (error/=0) stop
      call self%cli%get(switch='--fast', val=self%is_fast, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--iterations', val=self%implicit_iterations, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--stages', val=self%stages, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--steps', val=self%steps_max, error=error) ; if (error/=0) stop
      call self%cli%get(switch='-tf', val=self%final_time, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--w-scheme', val=self%w_scheme, error=error) ; if (error/=0) stop
      call self%cli%get(switch='--weno-order', val=self%weno_order, error=error) ; if (error/=0) stop
      call self%cli%get(switch='--weno-eps', val=self%weno_eps, error=error) ; if (error/=0) stop
      call self%cli%get(switch='-a', val=self%a, error=error) ; if (error/=0) stop
      call self%cli%get(switch='--Ni', val=self%Ni, error=error) ; if (error/=0) stop
      call self%cli%get(switch='-r', val=self%results, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--output', val=self%output, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(switch='--verbose',    val=verbose,    error=error) ; if (error/=0) stop

      if (self%final_time > 0._R_P) self%steps_max = 0

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

      allocate(initial_state(1:self%Ni))
      call square_wave_initial_state(initial_state=initial_state, BC_L=self%BC_L, BC_R=self%BC_R, Dx=1._R_P / self%Ni)
      call self%advection%initialize(Ni=Ni, Dx=1._R_P / self%Ni,     &
                                     BC_L=self%BC_L, BC_R=self%BC_R, &
                                     initial_state=initial_state,    &
                                     w_scheme=self%w_scheme,         &
                                     a=self%a,                       &
                                     weno_order=self%weno_order)
      endsubroutine parse_cli
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

   subroutine integrate(scheme, integrand_0, CFL, final_step, final_time, iterations, stages, is_fast, &
                        save_results, save_frequency, output_base_name)
   !< Integrate integrand by means of the given scheme.
   character(*),               intent(in)  :: scheme           !< Selected scheme.
   type(integrand_ladvection), intent(in)  :: integrand_0      !< Initial conditions.
   real(R_P),                  intent(in)  :: CFL              !< CFL stability coefficient.
   integer(I_P),               intent(in)  :: final_step       !< Final integration step.
   real(R_P),                  intent(in)  :: final_time       !< Final integration time.
   integer(I_P),               intent(in)  :: iterations       !< Number of fixed point iterations.
   integer(I_P),               intent(in)  :: stages           !< Number of stages.
   logical,                    intent(in)  :: is_fast          !< Activate fast mode integration.
   logical,                    intent(in)  :: save_results     !< Save results.
   integer(I_P),               intent(in)  :: save_frequency   !< Save frequency as steps multiple.
   character(*),               intent(in)  :: output_base_name !< Base name of output results file.
   class(integrator_object), allocatable   :: integrator       !< The integrator.
   type(integrator_runge_kutta_ssp)        :: integrator_start !< The (auto) start integrator.
   type(integrand_ladvection)              :: integrand        !< Integrand.
   type(integrand_ladvection), allocatable :: previous(:)      !< Previous time steps solutions.
   type(integrand_ladvection), allocatable :: stage(:)         !< Runge-Kutta stages.
   type(integrand_ladvection), allocatable :: stage_start(:)   !< Runge-Kutta (autor) start stages.
   type(integrand_ladvection)              :: buffer           !< Buffer oscillation field.
   type(integrand_ladvection)              :: filter           !< Filter displacement.
   real(R_P), allocatable                  :: time(:)          !< Time.
   real(R_P), allocatable                  :: Dt(:)            !< Time steps.
   real(R_P)                               :: Dt_a             !< Autoadptive time step.
   integer                                 :: step             !< Time steps counter.
   integer                                 :: step_offset      !< Time steps counter offset for slicing previous data array.
   character(len=:), allocatable           :: output_file_name !< File name of output results file.

   integrand = integrand_0

   output_file_name = trim(adjustl(output_base_name))//'-'//trim(adjustl(scheme))//'-'//strz(n=integrand%Ni, zero_pad=10)//'.dat'

   call foodie_integrator_factory(scheme=scheme, integrator=integrator, stages=stages, tolerance=1e2_R_P)
   if (is_fast) call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integrator=integrator)

   if (integrator%is_multistep()) then
      allocate(previous(1:integrator%steps_number()))
      call integrator_start%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
      allocate(stage_start(1:integrator_start%stages_number()))
      allocate(time(0:1:integrator%steps_number()))
      allocate(Dt(1:1:integrator%steps_number()))
   endif
   if (integrator%is_multistage()) allocate(stage(1:integrator%stages_number()))

   step = 0
   time = 0._R_P
   Dt = 0._R_P
   if (save_results) call integrand%export_tecplot(file_name=output_file_name, t=time(0))
   select type(integrator)
   type is(integrator_adams_bashforth)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, stage=stage_start, Dt=Dt(step), t=time(step))
            previous(step) = integrand
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand,       &
                                              previous=previous, &
                                              buffer=buffer,     &
                                              Dt=Dt(step),       &
                                              t=time(step-integrator%steps_number():step-1))
            else
               call integrator%integrate(U=integrand,       &
                                         previous=previous, &
                                         Dt=Dt(step),       &
                                         t=time(step-integrator%steps_number():step-1))
            endif
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_adams_bashforth_moulton)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, stage=stage_start, Dt=Dt(step), t=time(step))
            previous(step) = integrand
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand,       &
                                              previous=previous, &
                                              buffer=buffer,     &
                                              Dt=Dt(step),       &
                                              t=time(step-integrator%steps_number():step-1))
            else
               call integrator%integrate(U=integrand,       &
                                         previous=previous, &
                                         Dt=Dt(step),       &
                                         t=time(step-integrator%steps_number():step-1))
            endif
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_adams_moulton)
      if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps_number()+1))
      if (integrator%steps_number()==0) then
         step_offset = 1                         ! for 0 step-(a convention)-solver offset is 1
      else
         step_offset = integrator%steps_number() ! for >0 step-solver offset is steps
      endif
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, stage=stage_start, Dt=Dt(step), t=time(step))
            previous(step) = integrand
         else
            if (is_fast) then
               if (iterations>1) then
                  call integrator%integrate_fast(U=integrand,                     &
                                                 previous=previous,               &
                                                 buffer=buffer,                   &
                                                 Dt=Dt(step),                     &
                                                 t=time(step-step_offset:step-1), &
                                                 iterations=iterations)
               else
                  call integrator%integrate_fast(U=integrand,       &
                                                 previous=previous, &
                                                 buffer=buffer,     &
                                                 Dt=Dt(step),       &
                                                 t=time(step-step_offset:step-1))
               endif
            else
               if (iterations>1) then
                  call integrator%integrate(U=integrand,                     &
                                            previous=previous,               &
                                            Dt=Dt(step),                     &
                                            t=time(step-step_offset:step-1), &
                                            iterations=iterations)
               else
                  call integrator%integrate(U=integrand,       &
                                            previous=previous, &
                                            Dt=Dt(step),       &
                                            t=time(step-step_offset:step-1))
               endif
            endif
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_back_df)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, stage=stage_start, Dt=Dt(step), t=time(step))
            previous(step) = integrand
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand,       &
                                              previous=previous, &
                                              buffer=buffer,     &
                                              Dt=Dt(step),       &
                                              t=time(step-integrator%steps_number():step-1))
            else
               call integrator%integrate(U=integrand,       &
                                         previous=previous, &
                                         Dt=Dt(step),       &
                                         t=time(step-integrator%steps_number():step-1))
            endif
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_euler_explicit)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (is_fast) then
            call integrator%integrate_fast(U=integrand, buffer=buffer, Dt=Dt(step), t=time(step))
         else
            call integrator%integrate(U=integrand, Dt=Dt(step), t=time(step))
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_leapfrog)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, stage=stage_start, Dt=Dt(step), t=time(step))
            previous(step) = integrand
         else
            if (index(scheme, 'raw') > 0 ) then
               if (is_fast) then
                 call integrator%integrate_fast(U=integrand, previous=previous, buffer=buffer, Dt=Dt(step), &
                                                t=time(step), filter=filter)
               else
                 call integrator%integrate(U=integrand, previous=previous, Dt=Dt(step), t=time(step), filter=filter)
               endif
            else
               if (is_fast) then
                 call integrator%integrate_fast(U=integrand, previous=previous, buffer=buffer, Dt=Dt(step), t=time(step))
               else
                 call integrator%integrate(U=integrand, previous=previous, Dt=Dt(step), t=time(step))
               endif
            endif
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_lmm_ssp)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, stage=stage_start, Dt=Dt(step), t=time(step))
            previous(step) = integrand
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand,       &
                                              previous=previous, &
                                              buffer=buffer,     &
                                              Dt=Dt(step),       &
                                              t=time(step-integrator%steps_number():step-1))
            else
               call integrator%integrate(U=integrand,       &
                                         previous=previous, &
                                         Dt=Dt(step),       &
                                         t=time(step-integrator%steps_number():step-1))
            endif
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_lmm_ssp_vss)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, stage=stage_start, Dt=Dt(step), t=time(step))
            previous(step) = integrand
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand,                                  &
                                              previous=previous,                            &
                                              buffer=buffer,                                &
                                              Dt=Dt(step-integrator%steps_number():step-1), &
                                              t=time(step-integrator%steps_number():step-1))
            else
               call integrator%integrate(U=integrand,                                  &
                                         previous=previous,                            &
                                         Dt=Dt(step-integrator%steps_number():step-1), &
                                         t=time(step-integrator%steps_number():step-1))
            endif
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_ms_runge_kutta_ssp)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (integrator%steps_number() >= step) then
            call integrator_start%integrate(U=integrand, stage=stage_start, Dt=Dt(step), t=time(step))
            previous(step) = integrand
         else
            if (is_fast) then
               call integrator%integrate_fast(U=integrand,       &
                                              previous=previous, &
                                              stage=stage,       &
                                              buffer=buffer,     &
                                              Dt=Dt(step),       &
                                              t=time(step-integrator%steps_number():step-1))
            else
               call integrator%integrate(U=integrand,       &
                                         previous=previous, &
                                         stage=stage,       &
                                         Dt=Dt(step),       &
                                         t=time(step-integrator%steps_number():step-1))
            endif
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_runge_kutta_emd)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         Dt_a = Dt(step)
         if (is_fast) then
            call integrator%integrate_fast(U=integrand, stage=stage, buffer=buffer, Dt=Dt_, t=time(step))
         else
            call integrator%integrate(U=integrand, stage=stage, Dt=Dt_a, t=time(step))
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_runge_kutta_ls)
      if (allocated(stage)) deallocate(stage) ; allocate(stage(1:integrator%registers_number()))
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (is_fast) then
            call integrator%integrate_fast(U=integrand, stage=stage, buffer=buffer, Dt=Dt(step), t=time(step))
         else
            call integrator%integrate(U=integrand, stage=stage, Dt=Dt(step), t=time(step))
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_runge_kutta_lssp)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (is_fast) then
            call integrator%integrate_fast(U=integrand, stage=stage, buffer=buffer, Dt=Dt(step), t=time(step))
         else
            call integrator%integrate(U=integrand, stage=stage, Dt=Dt(step), t=time(step))
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo

   type is(integrator_runge_kutta_ssp)
      do
         step = step + 1
         Dt(step) = integrand%dt(final_step=final_step, final_time=final_time, t=time(step), CFL=CFL)
         if (is_fast) then
            call integrator%integrate_fast(U=integrand, stage=stage, buffer=buffer, Dt=Dt(step), t=time(step))
         else
            call integrator%integrate(U=integrand, stage=stage, Dt=Dt(step), t=time(step))
         endif
         time(step) = time(step-1) + Dt(step)
         if (save_results.and.mod(step, save_frequency)==0) call integrand%export_tecplot(file_name=output_file_name, t=time(step))
         if ((time(step) == final_time).or.(step == final_step)) exit
      enddo
   endselect

   if (save_results) call integrand%export_tecplot(close_file=.true.)
   endsubroutine integrate
endmodule foodie_test_ladvection_test

program foodie_test_ladvection
!< Test FOODIE with the integration of 1D linear advection PDE.

use foodie_test_ladvection_test, only : ladvection_test

implicit none
type(ladvection_test) :: test !< Linear advection test.

call test%execute
endprogram foodie_test_ladvection
