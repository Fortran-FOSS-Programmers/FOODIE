!< Test FOODIE with the integration of Oscillation equations.

module oscillation_test_t
!< Oscillation test handler definition.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use flap, only : command_line_interface
use foodie, only : foodie_integrator,                  &
                   foodie_integrator_class_names,      &
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
use oscillation_oscillator, only : oscillator
use penf, only : I_P, R_P, FR_P, str, strz

implicit none
private
public :: oscillation_test

integer(I_P),  parameter :: space_dimension=2                                !< Space dimensions.
real(R_P),     parameter :: initial_state(1:space_dimension)=[0._R_P,1._R_P] !< Initial state.

type :: oscillation_test
  !< Class to handle oscillation test(s).
  !<
  !< Test is driven by the Command Line Interface (CLI) options.
  !<
  !< Test has only 1 public method `execute`: it executes test(s) accordingly to cli options.
  private
  type(command_line_interface) :: cli                      !< Command line interface handler.
  integer(I_P)                 :: error=0                  !< Error handler.
  real(R_P)                    :: frequency=0.0_R_P        !< Oscillation frequency.
  real(R_P)                    :: final_time=0.0_R_P       !< Final integration time.
  integer(I_P)                 :: implicit_iterations=0    !< Number of iterations (implicit solvers).
  integer(I_P)                 :: stages=0                 !< Number of stages.
  character(99)                :: output_cli='unset'       !< Output files basename.
  character(99)                :: scheme='adams_bashforth' !< Scheme used.
  real(R_P),    allocatable    :: Dt(:)                    !< Time step(s) exercised.
  real(R_P),    allocatable    :: tolerance(:)             !< Tolerance(s) exercised on local truncation error.
  logical                      :: errors_analysis=.false.  !< Flag for activating errors analysis.
  logical                      :: results=.false.          !< Flag for activating results saving.
  logical                      :: is_fast=.false.          !< Flag for activating fast schemes.
  contains
    ! public methods
    procedure, pass(self) :: execute !< Execute selected test(s).
    ! private methods
    procedure, pass(self), private :: initialize !< Initialize test: set Command Line Interface, parse it and check its validity.
    procedure, pass(self), private :: test       !< Perform the test.
endtype oscillation_test

abstract interface
  !< Interface for integrate procedures.
  subroutine integrate_interface(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  import :: I_P, R_P
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  endsubroutine integrate_interface
endinterface

contains
  ! public methods
  subroutine execute(self)
  !< Execute test(s).
  class(oscillation_test), intent(inout) :: self                  !< Test.
  character(99), allocatable             :: integrator_schemes(:) !< Name of FOODIE integrator schemes.
  integer(I_P)                           :: s                     !< Counter.

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
  class(oscillation_test), intent(inout) :: self !< Test.

  call set_cli
  call parse_cli
  contains
    subroutine set_cli()
    !< Set Command Line Interface.

    associate(cli => self%cli)
      call cli%init(progname    = 'oscillation',                                              &
                    authors     = 'Fortran-FOSS-Programmers',                                 &
                    license     = 'GNU GPLv3',                                                &
                    description = 'Test FOODIE library on Oscillation equations integration', &
                    examples    = ["oscillation --scheme euler_explicit --save_results  ",    &
                                   "oscillation --scheme all -r                         "])
      call cli%add(switch='--scheme', switch_ab='-s', help='integrator scheme used', required=.false., def='all', act='store')
      call cli%add(switch='--fast', help='activate fast solvers', required=.false., act='store_true', def='.false.')
      call cli%add(switch='--iterations', help='number of iterations for implicit schemes', required=.false., act='store', def='5')
      call cli%add(switch='--frequency', switch_ab='-f', help='oscillation frequency', required=.false., def='1e-4', act='store')
      call cli%add(switch='--time_step', switch_ab='-Dt', nargs='+', help='time step', required=.false., def='100.d0', act='store')
      call cli%add(switch='--tolerance', switch_ab='-tol', nargs='+', help='error tolerance', required=.false., def='0.001d0', &
                   act='store')
      call cli%add(switch='--stages', help='stages number', required=.false., def='2', act='store')
      call cli%add(switch='--t_final', switch_ab='-tf', help='final integration time', required=.false., def='1e6', act='store')
      call cli%add(switch='--save_results', switch_ab='-r', help='save results', required=.false., act='store_true', def='.false.')
      call cli%add(switch='--output', help='output files basename', required=.false., act='store', def='unset')
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
    call self%cli%get(switch='-f', val=self%frequency, error=self%error) ; if (self%error/=0) stop
    call self%cli%get_varying(switch='-Dt', val=self%Dt, error=self%error) ; if (self%error/=0) stop
    call self%cli%get_varying(switch='-tol', val=self%tolerance, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-tf', val=self%final_time, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-r', val=self%results, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='--output', val=self%output_cli, error=self%error) ; if (self%error/=0) stop
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
      if ((.not.(size(self%Dt, dim=1) > 1)).and.(.not.(size(self%tolerance, dim=1) > 1))) self%errors_analysis = .false.
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
  class(oscillation_test), intent(in)      :: self                        !< Test.
  character(*),            intent(in)      :: scheme                      !< Selected scheme.
  procedure(integrate_interface), pointer  :: integrate => null()         !< Actual integrate procedure.
  type(integrator_runge_kutta_emd)         :: int_runge_kutta_emd         !< Runge Kutta embedded integrator.
  real(R_P), allocatable                   :: solution(:,:)               !< Solution at each time step.
  real(R_P), allocatable                   :: error(:,:)                  !< Error (norm L2) with respect the exact solution.
  real(R_P), allocatable                   :: order(:,:)                  !< Observed order based on subsequent refined solutions.
  real(R_P), allocatable                   :: Dt_mean(:)                  !< Mean time steps used for adaptive solver solutions.
  integer(I_P)                             :: last_step                   !< Last time step computed.
  integer(I_P)                             :: t                           !< Counter.

  print "(A)", trim(adjustl(scheme))
  call initialize_test
  if (index(trim(adjustl(scheme)), trim(int_runge_kutta_emd%class_name())) > 0) then
    do t=1, size(self%tolerance, dim=1)
      call integrate(scheme=scheme,              &
                     frequency=self%frequency,   &
                     final_time=self%final_time, &
                     solution=solution,          &
                     error=error(:, t),          &
                     last_step=last_step,        &
                     tolerance=self%tolerance(t))
      Dt_mean(t) = self%final_time/real(last_step, kind=R_P)
      if (allocated(solution)) then
        print "(A,I10,A,F10.3,A,F10.3,A,E10.3,A,E10.3)", "    steps: ", last_step,              &
                                                         "    Dt: ", Dt_mean(t),                &
                                                         ", f*Dt: ", self%frequency*Dt_mean(t), &
                                                         ", E(x): ", error(1, t), ", E(y): ", error(2, t)
        if (self%errors_analysis.and.t>1) then
          order(:, t-1) = observed_order(error=error(:, t-1:t), Dt=Dt_mean(t-1:t))
          print "(A,F10.2,A,F10.2)", "      Observed order, O(x): ", order(1, t-1), ", O(y): " , order(2, t-1)
        endif
        call save_results(results=self%results,         &
                          output_cli=self%output_cli,   &
                          scheme=trim(adjustl(scheme)), &
                          frequency=self%frequency,     &
                          solution=solution(:, 0:last_step))
      endif
    enddo
  else
    do t=1, size(self%Dt, dim=1)
      call integrate(scheme=scheme,              &
                     frequency=self%frequency,   &
                     final_time=self%final_time, &
                     solution=solution,          &
                     error=error(:, t),          &
                     last_step=last_step,        &
                     Dt=self%Dt(t),              &
                     stages=self%stages,         &
                     iterations=self%implicit_iterations)
      if (allocated(solution)) then
        print "(A,I10,A,F10.3,A,F10.3,A,E10.3,A,E10.3)", "    steps: ", last_step,              &
                                                         "    Dt: ", self%Dt(t),                &
                                                         ", f*Dt: ", self%frequency*self%Dt(t), &
                                                         ", E(x): ", error(1, t), ", E(y): ", error(2, t)
        if (self%errors_analysis.and.t>1) then
          order(:, t-1) = observed_order(error=error(:, t-1:t), Dt=self%Dt(t-1:t))
          print "(A,F10.2,A,F10.2)", "      Observed order, O(x): ", order(1, t-1), ", O(y): " , order(2, t-1)
        endif
        call save_results(results=self%results,         &
                          output_cli=self%output_cli,   &
                          scheme=trim(adjustl(scheme)), &
                          frequency=self%frequency,     &
                          solution=solution(:, 0:last_step))
      endif
    enddo
  endif
  contains
    subroutine initialize_test
    !< Intialize test.
    type(integrator_adams_bashforth)         :: int_adams_bashforth         !< Adams Bashforth integrator.
    type(integrator_adams_bashforth_moulton) :: int_adams_bashforth_moulton !< Adams Bashforth Moulton  integrator.
    type(integrator_adams_moulton)           :: int_adams_moulton           !< Adams Moulton integrator.
    type(integrator_back_df)                 :: int_back_df                 !< Back differentiation formula integrator.
    type(integrator_euler_explicit)          :: int_euler_explicit          !< Euler explicit integrator.
    type(integrator_leapfrog)                :: int_leapfrog                !< Leapfrog integrator.
    type(integrator_lmm_ssp)                 :: int_lmm_ssp                 !< LMM SSP integrator.
    type(integrator_lmm_ssp_vss)             :: int_lmm_ssp_vss             !< LMM SSP VSS integrator.
    type(integrator_ms_runge_kutta_ssp)      :: int_ms_runge_kutta_ssp      !< Multistep Runge-Kutta SSP integrator.
    type(integrator_runge_kutta_ls)          :: int_runge_kutta_ls          !< Runge-Kutta low storage integrator.
    type(integrator_runge_kutta_lssp)        :: int_runge_kutta_lssp        !< Linear Runge-Kutta SSP integrator.
    type(integrator_runge_kutta_ssp)         :: int_runge_kutta_ssp         !< Runge Kutta-SSP integrator.

    if     (index(trim(adjustl(scheme)), trim(int_adams_bashforth_moulton%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_adams_bashforth_moulton)
          ! integrate => integrate_adams_bashforth_moulton_fast
       else
          integrate => integrate_adams_bashforth_moulton
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_adams_bashforth%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_adams_bashforth)
          integrate => integrate_adams_bashforth_fast
       else
          integrate => integrate_adams_bashforth
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_adams_moulton%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_adams_moulton)
          ! integrate => integrate_adams_moulton_fast
       else
          integrate => integrate_adams_moulton
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_back_df%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_back_df)
          integrate => integrate_back_df_fast
       else
          integrate => integrate_back_df
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_euler_explicit%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_euler_explicit)
          integrate => integrate_euler_explicit_fast
       else
          integrate => integrate_euler_explicit
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_leapfrog%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_leapfrog)
          integrate => integrate_leapfrog_fast
       else
          integrate => integrate_leapfrog
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_lmm_ssp_vss%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_lmm_ssp_vss)
          integrate => integrate_lmm_ssp_vss_fast
       else
          integrate => integrate_lmm_ssp_vss
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_lmm_ssp%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_lmm_ssp)
          integrate => integrate_lmm_ssp_fast
       else
          integrate => integrate_lmm_ssp
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_ms_runge_kutta_ssp%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_ms_runge_kutta_ssp)
          integrate => integrate_ms_runge_kutta_ssp_fast
       else
          integrate => integrate_ms_runge_kutta_ssp
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_runge_kutta_emd%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_runge_kutta_emd)
          integrate => integrate_runge_kutta_emd_fast
       else
          integrate => integrate_runge_kutta_emd
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_runge_kutta_lssp%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_runge_kutta_lssp)
          integrate => integrate_runge_kutta_lssp_fast
       else
          integrate => integrate_runge_kutta_lssp
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_runge_kutta_ls%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_runge_kutta_ls)
          integrate => integrate_runge_kutta_ls_fast
       else
          integrate => integrate_runge_kutta_ls
       endif
    elseif (index(trim(adjustl(scheme)), trim(int_runge_kutta_ssp%class_name())) > 0) then
       if (self%is_fast) then
          call check_scheme_has_fast_mode(scheme=trim(adjustl(scheme)), integr=int_runge_kutta_ssp)
          integrate => integrate_runge_kutta_ssp_fast
       else
          integrate => integrate_runge_kutta_ssp
       endif
    endif

    if (index(trim(adjustl(scheme)), trim(int_runge_kutta_emd%class_name())) > 0) then
      if (allocated(error)) deallocate(error) ; allocate(error(1:space_dimension, 1:size(self%tolerance)))
      if (allocated(Dt_mean)) deallocate(Dt_mean) ; allocate(Dt_mean(1:size(error, dim=2)))
    else
      if (allocated(error)) deallocate(error) ; allocate(error(1:space_dimension, 1:size(self%Dt)))
    endif

    error = 0.0_R_P
    if (self%errors_analysis) then
      if (allocated(order)) deallocate(order) ; allocate(order(1:space_dimension, 1:size(error, dim=2)-1))
    endif
    endsubroutine initialize_test
  endsubroutine test

  ! non type bound procedures
  subroutine check_scheme_has_fast_mode(scheme, integr)
  !< Check if a scheme support fast mode integrate.
  character(*),             intent(in) :: scheme !< Scheme name.
  class(integrator_object), intent(in) :: integr !< Integrator instance.

  if (.not.integr%has_fast_mode()) then
     write(stderr, '(A)') 'error: '//trim(adjustl(scheme))//' has not fast integrate mode!'
     stop
  endif
  endsubroutine check_scheme_has_fast_mode

  subroutine integrate_adams_bashforth(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, &
                                       stages)
  !< Integrate domain by means of the Adams-Bashforth scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_adams_bashforth)             :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      call integrator%integrate(U=domain,          &
                                previous=previous, &
                                Dt=Dt,             &
                                t=solution(0, step-integrator%steps:step-1))
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_adams_bashforth

  subroutine integrate_adams_bashforth_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, &
                                       stages)
  !< Integrate domain by means of the Adams-Bashforth scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_adams_bashforth)             :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      call integrator%integrate_fast(U=domain,          &
                                     previous=previous, &
                                     buffer=buffer,     &
                                     Dt=Dt,             &
                                     t=solution(0, step-integrator%steps:step-1))
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_adams_bashforth_fast

  subroutine integrate_adams_bashforth_moulton(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, &
                                               tolerance, stages)
  !< Integrate domain by means of the Adams-Bashforth-Moulton scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_adams_bashforth_moulton)     :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      call integrator%integrate(U=domain,          &
                                previous=previous, &
                                Dt=Dt,             &
                                t=solution(0, step-integrator%steps:step-1))
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_adams_bashforth_moulton

  subroutine integrate_adams_moulton(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the Adams-Moulton scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_adams_moulton)               :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  integer                                      :: step          !< Time steps counter.
  integer                                      :: step_offset   !< Time steps counter offset for slicing previous data array.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps+1))
  if (integrator%steps==0) then
    step_offset = 1                ! for 0 step-(a convention)-solver offset is 1
  else
    step_offset = integrator%steps ! for >0 step-solver offset is steps
  endif

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = domain
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

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_adams_moulton

  subroutine integrate_back_df(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the back differentiation formula scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_back_df)                     :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      if (iterations>1) then
        call integrator%integrate(U=domain,                                   &
                                  previous=previous,                          &
                                  Dt=Dt,                                      &
                                  t=solution(0,step-integrator%steps:step-1), &
                                  iterations=iterations)
      else
        call integrator%integrate(U=domain,                              &
                                  previous=previous,                     &
                                  Dt=Dt,                                 &
                                  t=solution(0,step-integrator%steps:step-1))
      endif
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_back_df

  subroutine integrate_back_df_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the back differentiation formula scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_back_df)                     :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate_fast(U=domain, stage=rk_stage, buffer=buffer, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      if (iterations>1) then
        call integrator%integrate_fast(U=domain,                                   &
                                       previous=previous,                          &
                                       buffer=buffer,                              &
                                       Dt=Dt,                                      &
                                       t=solution(0,step-integrator%steps:step-1), &
                                       iterations=iterations)
      else
        call integrator%integrate_fast(U=domain,                              &
                                       previous=previous,                     &
                                       buffer=buffer,                         &
                                       Dt=Dt,                                 &
                                       t=solution(0,step-integrator%steps:step-1))
      endif
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_back_df_fast

  subroutine integrate_euler_explicit(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the Euler explicit scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_euler_explicit)              :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate(U=domain, Dt=Dt, t=solution(0, step))

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_euler_explicit

  subroutine integrate_euler_explicit_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, &
                                           tolerance, stages)
  !< Integrate domain by means of the Euler explicit scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_euler_explicit)              :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate_fast(U=domain, buffer=buffer, Dt=Dt, t=solution(0, step))

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_euler_explicit_fast

  subroutine integrate_leapfrog(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the leapfrog scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_leapfrog)                    :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  integer                                      :: step          !< Time steps counter.
  type(oscillator)                             :: filter        !< Filter displacement.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  if (index(scheme, 'raw') > 0 ) then
    do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
      step = step + 1

      if (integrator%steps >= step) then
        call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
        previous(step) = domain
      else
        call integrator%integrate(U=domain, previous=previous, Dt=Dt, t=solution(0, step), filter=filter)
      endif

      solution(0, step) = step * Dt

      solution(1:, step) = domain%output()
    enddo
  else
    do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
      step = step + 1

      if (integrator%steps >= step) then
        call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
        previous(step) = domain
      else
        call integrator%integrate(U=domain, previous=previous, Dt=Dt, t=solution(0, step))
      endif

      solution(0, step) = step * Dt

      solution(1:, step) = domain%output()
    enddo
  endif
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_leapfrog

  subroutine integrate_leapfrog_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the leapfrog scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_leapfrog)                    :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  type(oscillator)                             :: filter        !< Filter displacement.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  if (index(scheme, 'raw') > 0 ) then
    do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
      step = step + 1

      if (integrator%steps >= step) then
        call integrator_rk%integrate_fast(U=domain, stage=rk_stage, buffer=buffer, Dt=Dt, t=solution(0, step))
        previous(step) = domain
      else
        call integrator%integrate_fast(U=domain, previous=previous, buffer=buffer, Dt=Dt, t=solution(0, step), filter=filter)
      endif

      solution(0, step) = step * Dt

      solution(1:, step) = domain%output()
    enddo
  else
    do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
      step = step + 1

      if (integrator%steps >= step) then
        call integrator_rk%integrate_fast(U=domain, stage=rk_stage, buffer=buffer, Dt=Dt, t=solution(0, step))
        previous(step) = domain
      else
        call integrator%integrate_fast(U=domain, previous=previous, buffer=buffer, Dt=Dt, t=solution(0, step))
      endif

      solution(0, step) = step * Dt

      solution(1:, step) = domain%output()
    enddo
  endif
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_leapfrog_fast

  subroutine integrate_lmm_ssp(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the LLM SSP scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_lmm_ssp)                     :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      call integrator%integrate(U=domain,          &
                                previous=previous, &
                                Dt=Dt,             &
                                t=solution(0, step-integrator%steps:step-1))
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_lmm_ssp

  subroutine integrate_lmm_ssp_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the LLM SSP scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_lmm_ssp)                     :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate_fast(U=domain, stage=rk_stage, buffer=buffer, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      call integrator%integrate_fast(U=domain,          &
                                     previous=previous, &
                                     buffer=buffer,     &
                                     Dt=Dt,             &
                                     t=solution(0, step-integrator%steps:step-1))
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_lmm_ssp_fast

  subroutine integrate_lmm_ssp_vss(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the LLM SSP variable step size scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_lmm_ssp_vss)                 :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  integer                                      :: step          !< Time steps counter.
  real(R_P), allocatable                       :: Dts(:)        !< Time steps for variable stepsize.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))
  if (allocated(Dts)) deallocate(Dts) ; allocate(Dts(1:integrator%steps)) ; Dts = Dt

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      call integrator%integrate(U=domain,          &
                                previous=previous, &
                                Dt=Dts,            &
                                t=solution(0, step-integrator%steps:step-1))
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_lmm_ssp_vss

  subroutine integrate_lmm_ssp_vss_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, &
                                        tolerance, stages)
  !< Integrate domain by means of the LLM SSP variable step size scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_lmm_ssp_vss)                 :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.
  real(R_P), allocatable                       :: Dts(:)        !< Time steps for variable stepsize.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))
  if (allocated(Dts)) deallocate(Dts) ; allocate(Dts(1:integrator%steps)) ; Dts = Dt

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate_fast(U=domain, stage=rk_stage, buffer=buffer, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      call integrator%integrate_fast(U=domain,          &
                                     previous=previous, &
                                     buffer=buffer,     &
                                     Dt=Dts,            &
                                     t=solution(0, step-integrator%steps:step-1))
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_lmm_ssp_vss_fast

  subroutine integrate_ms_runge_kutta_ssp(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, &
                                          stages)
  !< Integrate domain by means of the Multistep Runge-Kutta SSP scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_ms_runge_kutta_ssp)          :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages for the RK initial integrator.
  type(oscillator), allocatable                :: stage(:)      !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))
  if (allocated(stage)) deallocate(stage) ; allocate(stage(1:integrator%stages))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      call integrator%integrate(U=domain,          &
                                previous=previous, &
                                stage=stage,       &
                                Dt=Dt,             &
                                t=solution(0, step-integrator%steps:step-1))
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_ms_runge_kutta_ssp

  subroutine integrate_ms_runge_kutta_ssp_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, &
                                               tolerance, stages)
  !< Integrate domain by means of the Multistep Runge-Kutta SSP scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_ms_runge_kutta_ssp)          :: integrator    !< The integrator.
  type(integrator_runge_kutta_ssp)             :: integrator_rk !< RK integrator for starting non self-starting integrators.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages for the RK initial integrator.
  type(oscillator), allocatable                :: stage(:)      !< Runge-Kutta stages.
  type(oscillator), allocatable                :: previous(:)   !< Previous time steps solutions.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(previous)) deallocate(previous) ; allocate(previous(1:integrator%steps))
  if (allocated(stage)) deallocate(stage) ; allocate(stage(1:integrator%stages))

  call integrator_rk%initialize(scheme='runge_kutta_ssp_stages_5_order_4')
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator_rk%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    if (integrator%steps >= step) then
      call integrator_rk%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = domain
    else
      call integrator%integrate_fast(U=domain,          &
                                     previous=previous, &
                                     stage=stage,       &
                                     buffer=buffer,     &
                                     Dt=Dt,             &
                                     t=solution(0, step-integrator%steps:step-1))
    endif

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_ms_runge_kutta_ssp_fast

  subroutine integrate_runge_kutta_emd(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the Runge Kutta embedded scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_runge_kutta_emd)             :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  integer                                      :: step          !< Time steps counter.
  real(R_P)                                    :: Dt_a          !< Adaptive time step.

  call domain%init(initial_state=initial_state, frequency=frequency)

  Dt_a = 10000._R_P ! initial step very large to trigger adaptation
  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/10))) ! hope is enough

  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme, tolerance=tolerance)
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate(U=domain, stage=rk_stage, Dt=Dt_a, t=solution(0, step))

    solution(0, step) = solution(0, step - 1) + Dt_a

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_runge_kutta_emd

  subroutine integrate_runge_kutta_emd_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, &
                                            tolerance, stages)
  !< Integrate domain by means of the Runge Kutta embedded scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_runge_kutta_emd)             :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.
  real(R_P)                                    :: Dt_a          !< Adaptive time step.

  call domain%init(initial_state=initial_state, frequency=frequency)

  Dt_a = 10000._R_P ! initial step very large to trigger adaptation
  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/10))) ! hope is enough

  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme, tolerance=tolerance)
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate_fast(U=domain, stage=rk_stage, buffer=buffer, Dt=Dt_a, t=solution(0, step))

    solution(0, step) = solution(0, step - 1) + Dt_a

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_runge_kutta_emd_fast

  subroutine integrate_runge_kutta_ls(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the Runge Kutta low storage scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_runge_kutta_ls)              :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator%registers_number()))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_runge_kutta_ls

  subroutine integrate_runge_kutta_ls_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, &
                                           tolerance, stages)
  !< Integrate domain by means of the Runge Kutta low storage scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_runge_kutta_ls)              :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator%registers_number()))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate_fast(U=domain, stage=rk_stage, buffer=buffer, Dt=Dt, t=solution(0, step))

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_runge_kutta_ls_fast

  subroutine integrate_runge_kutta_lssp(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, &
                                        stages)
  !< Integrate domain by means of the Linear Runge Kutta SSP scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_runge_kutta_lssp)            :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme, stages=stages)
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_runge_kutta_lssp

  subroutine integrate_runge_kutta_lssp_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, &
                                             tolerance, stages)
  !< Integrate domain by means of the Linear Runge Kutta SSP scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_runge_kutta_lssp)            :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme, stages=stages)
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate_fast(U=domain, stage=rk_stage, buffer=buffer, Dt=Dt, t=solution(0, step))

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_runge_kutta_lssp_fast

  subroutine integrate_runge_kutta_ssp(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, tolerance, stages)
  !< Integrate domain by means of the Runge Kutta SSP scheme.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_runge_kutta_ssp)             :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate(U=domain, stage=rk_stage, Dt=Dt, t=solution(0, step))

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_runge_kutta_ssp

  subroutine integrate_runge_kutta_ssp_fast(scheme, frequency, final_time, solution, error, last_step, iterations, Dt, &
                                            tolerance, stages)
  !< Integrate domain by means of the Runge Kutta SSP scheme, fast mode.
  character(*),           intent(in)           :: scheme        !< Selected scheme.
  real(R_P),              intent(in)           :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)           :: final_time    !< Final integration time.
  real(R_P), allocatable, intent(out)          :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out)          :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out)          :: last_step     !< Last time step computed.
  integer(I_P),           intent(in), optional :: iterations    !< Number of fixed point iterations.
  real(R_P),              intent(in), optional :: Dt            !< Time step.
  real(R_P),              intent(in), optional :: tolerance     !< Local error tolerance.
  integer(I_P),           intent(in), optional :: stages        !< Number of stages.
  type(integrator_runge_kutta_ssp)             :: integrator    !< The integrator.
  type(oscillator)                             :: domain        !< Oscillation field.
  type(oscillator), allocatable                :: rk_stage(:)   !< Runge-Kutta stages.
  type(oscillator)                             :: buffer        !< Buffer oscillation field.
  integer                                      :: step          !< Time steps counter.

  call domain%init(initial_state=initial_state, frequency=frequency)

  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  solution = 0.0_R_P
  solution(1:, 0) = domain%output()

  call integrator%initialize(scheme=scheme)
  if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:integrator%stages))

  step = 0
  do while(solution(0, step) < final_time .and. step < ubound(solution, dim=2))
    step = step + 1

    call integrator%integrate_fast(U=domain, stage=rk_stage, buffer=buffer, Dt=Dt, t=solution(0, step))

    solution(0, step) = step * Dt

    solution(1:, step) = domain%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  endsubroutine integrate_runge_kutta_ssp_fast

  subroutine save_results(results, output_cli, scheme, frequency, solution)
  !< Save results (and plots).
  logical,      intent(in)      :: results                   !< Flag for activating results saving.
  character(*), intent(in)      :: output_cli                !< Output files basename coming from CLI.
  character(*), intent(in)      :: scheme                    !< Selected scheme: must be defined into *solvers*.
  real(R_P),    intent(in)      :: frequency                 !< Oscillation frequency.
  real(R_P),    intent(IN)      :: solution(0:,0:)           !< Solution at each time step.
  character(len=:), allocatable :: title                     !< Output files title.
  character(len=:), allocatable :: basename                  !< Output files basename.
  integer(I_P)                  :: rawfile                   !< Raw file unit for saving results.
  real(R_P)                     :: ex_sol(1:space_dimension) !< Exact solution.
  integer(I_P)                  :: s                         !< Counter.

  title = 'Oscillation equations integration, solver='//trim(adjustl(scheme))
  if (trim(adjustl(output_cli))/='unset') then
    basename = trim(adjustl(output_cli))//'-'//trim(strz(ubound(solution, dim=2), 10))//'-time_steps-'//trim(adjustl(scheme))
  else
    basename = 'oscillation_test-'//trim(strz(ubound(solution, dim=2), 10))//'-time_steps-'//trim(adjustl(scheme))
  endif
  if (results) then
    open(newunit=rawfile, file=basename//'.dat')
    write(rawfile, '(A)')'TITLE="'//title//'"'
    write(rawfile, '(A)')'VARIABLES="t" "x" "y" "amplitude" "phase"'
    write(rawfile, '(A)')'ZONE T="'//trim(adjustl(scheme))//'"'
    do s=0, ubound(solution, dim=2)
      write(rawfile, '(5('//FR_P//',1X))')solution(:, s), amplitude_phase(solution(1:2, s))
    enddo
    write(rawfile, '(A)')'ZONE T="Exact solution"'
    do s=0, ubound(solution, dim=2)
      ex_sol = exact_solution_xy(frequency=frequency, t=solution(0, s))
      write(rawfile, '(5('//FR_P//',1X))')solution(0, s), ex_sol, amplitude_phase(ex_sol)
    enddo
    close(rawfile)
  endif
  contains
    function amplitude_phase(sol) result(ap)
    !< Compute amplitude and phase of the solution provided in X-Y domain.
    real(R_P), intent(IN) :: sol(1:) !< Solution in X-Y domain.
    real(R_P)             :: ap(1:2) !< Amplitude and phase solution.

    ap(1) = sqrt(sol(1)**2 + sol(2)**2)
    ap(2) = atan(-sol(1) / sol(2))
    endfunction amplitude_phase
  endsubroutine save_results

  pure function exact_solution_xy(frequency, t) result(ex_sol)
  !< Compute the exact solution on X-Y domain.
  real(R_P), intent(in) :: frequency   !< Oscillation frequency.
  real(R_P), intent(IN) :: t           !< Time.
  real(R_P)             :: ex_sol(1:2) !< Exact solution.

  ex_sol(1) = initial_state(1) * cos(frequency * t) - initial_state(2) * sin(frequency * t)
  ex_sol(2) = initial_state(1) * sin(frequency * t) + initial_state(2) * cos(frequency * t)
  endfunction exact_solution_xy

  pure function error_L2(frequency, solution)
  !< Compute the L2 norm of numerical error with respect the exact solution.
  real(R_P), intent(in) :: frequency                   !< Oscillation frequency.
  real(R_P), intent(IN) :: solution(0:,0:)             !< Solution at each time step.
  real(R_P)             :: error_L2(1:space_dimension) !< L2 norm of the numerical error.
  real(R_P)             :: ex_sol(1:space_dimension)   !< Exact solution.
  integer(I_P)          :: s                           !< Steps/stages counter.
  integer(I_P)          :: v                           !< Variables counter.

  error_L2 = 0._R_P
  do s=0, ubound(solution, dim=2)
    ex_sol = exact_solution_xy(frequency=frequency, t=solution(0, s))
    do v=1, space_dimension
      error_L2(v) = error_L2(v) + (solution(v, s) - ex_sol(v))*(solution(v, s) - ex_sol(v))
    enddo
  enddo
  error_L2 = sqrt(error_L2)
  endfunction error_L2

  pure function observed_order(error, Dt)
  !< Estimate the order of accuracy using 2 subsequent refined numerical solutions.
  real(R_P), intent(IN) :: error(1:space_dimension, 1:2)     !< Computed errors.
  real(R_P), intent(IN) :: Dt(1:2)                           !< Time steps used.
  real(R_P)             :: observed_order(1:space_dimension) !< Estimation of the order of accuracy.
  integer(I_P)          :: v                                 !< Variables counter.

  do v=1, space_dimension
    observed_order(v) = log(error(v, 1) / error(v, 2)) / log(Dt(1) / Dt(2))
  enddo
  endfunction observed_order
endmodule oscillation_test_t

program integrate_oscillation
!< Test FOODIE with the integration of Oscillation equations.

use oscillation_test_t, only : oscillation_test

implicit none
type(oscillation_test) :: test !< Oscillation test.

call test%execute
endprogram integrate_oscillation
