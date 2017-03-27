!< Test FOODIE with the integration of Oscillation equations.
module oscillation_test_t
!-----------------------------------------------------------------------------------------------------------------------------------
!< Oscillation test handler definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use flap, only : command_line_interface
use foodie, only : integrator_adams_bashforth,         &
                   integrator_adams_bashforth_moulton, &
                   integrator_adams_moulton,           &
                   integrator_back_df,                 &
                   integrator_euler_explicit,          &
                   integrator_leapfrog,                &
                   integrator_runge_kutta_emd,         &
                   integrator_runge_kutta_ls,          &
                   integrator_runge_kutta_tvd
use oscillation_t, only : oscillation
use penf, only : I_P, R_P, FR_P, str, strz
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: oscillation_test
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
integer,       parameter :: space_dimension=2                                !< Space dimensions.
real(R_P),     parameter :: initial_state(1:space_dimension)=[0._R_P,1._R_P] !< Initial state.
character(99), parameter :: solvers(1:11) = ["all                    ", &
                                             "adams-bashforth        ", &
                                             "adams-bashforth-moulton", &
                                             "adams-moulton          ", &
                                             "backward-diff-formula  ", &
                                             "emd-runge-kutta        ", &
                                             "euler                  ", &
                                             "leapfrog               ", &
                                             "leapfrog-raw           ", &
                                             "ls-runge-kutta         ", &
                                             "tvd-runge-kutta        "]      !< List of available solvers.

type :: oscillation_test
  !< Class to handle oscillation test(s).
  !<
  !< Test is driven by the Command Line Interface (CLI) options.
  !<
  !< Test has only 1 public method `execute`: it executes test(s) accordingly to cli options.
  private
  type(command_line_interface) :: cli                      !< Command line interface handler.
  integer(I_P)                 :: error=0                  !< Error handler.
  logical                      :: errors_analysis=.false.  !< Flag for activating errors analysis.
  real(R_P)                    :: frequency=0.0_R_P        !< Oscillation frequency.
  real(R_P)                    :: final_time=0.0_R_P       !< Final integration time.
  integer(I_P)                 :: implicit_iterations=0    !< Number of iterations (implicit solvers).
  character(99)                :: output_cli='unset'       !< Output files basename.
  logical                      :: plots=.false.            !< Flag for activating plots saving.
  logical                      :: results=.false.          !< Flag for activating results saving.
  character(99)                :: solver='adams-bashforth' !< Solver used.
  integer(I_P), allocatable    :: stages_steps(:)          !< Number of stages/steps used.
  real(R_P),    allocatable    :: Dt(:)                    !< Time step(s) exercised.
  real(R_P),    allocatable    :: tolerance(:)             !< Tolerance(s) exercised on local truncation error.
  contains
    private
    ! Public methods
    procedure, pass(self), public  :: execute !< Execute selected test(s).
    ! Private methods
    procedure, pass(self), private :: init    !< Initialize test: set Command Line Interface, parse it and check its validity.
    procedure, pass(self), private :: test    !< Perform the test.
endtype oscillation_test
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! Public methods
  subroutine execute(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Execute test(s).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation_test), intent(inout) :: self !< Test.
  integer(I_P)                           :: s    !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%init
  if (trim(adjustl(self%solver))/='all') then
    call self%test(solver=self%solver)
  else
    do s=2, ubound(solvers, dim=1)
      self%solver = solvers(s)
      call self%test(solver=self%solver)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine execute

  ! Private methods
  subroutine init(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Initialize test: set Command Line Interface, parse it and check its validity.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation_test), intent(inout) :: self !< Test.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call set_cli
  call parse_cli
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    subroutine set_cli()
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Set Command Line Interface.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    associate(cli => self%cli)
      call cli%init(progname    = 'oscillation',                                              &
                    authors     = 'Fortran-FOSS-Programmers',                                 &
                    license     = 'GNU GPLv3',                                                &
                    description = 'Test FOODIE library on Oscillation equations integration', &
                    examples    = ["oscillation --solver euler --results  ",                  &
                                   "oscillation --solver ls-runge-kutta -r",                  &
                                   "oscillation --solver adams-bashforth  ",                  &
                                   "oscillation --solver all --plots -r   "])
      call cli%add(switch='--solver', switch_ab='-s', help='ODE solver used', required=.false., def='all', act='store')
      call cli%add(switch='--iterations', help='Number of iterations for implicit solvers', required=.false., act='store', def='5')
      call cli%add(switch='--frequency', switch_ab='-f', help='Oscillation frequency', required=.false., def='1e-4', act='store')
      call cli%add(switch='--ss', nargs='+', help='Stages/steps used', required=.false., def='-1', act='store')
      call cli%add(switch='--time_step', switch_ab='-Dt', nargs='+', help='Time step', required=.false., def='100.d0', act='store')
      call cli%add(switch='--tolerance', switch_ab='-tol', nargs='+', help='Error Tolerance', required=.false., def='0.001d0', &
                   act='store')
      call cli%add(switch='--t_final', switch_ab='-tf', help='Final integration time', required=.false., def='1e6', act='store')
      call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.')
      call cli%add(switch='--plots', switch_ab='-p', help='Save plots', required=.false., act='store_true', def='.false.')
      call cli%add(switch='--output', help='Output files basename', required=.false., act='store', def='unset')
      call cli%add(switch='--errors_analysis', help='Peform errors analysis', required=.false., act='store_true', def='.false.')
    endassociate
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_cli

    subroutine parse_cli()
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Parse Command Line Interface and check its validity.
    !-------------------------------------------------------------------------------------------------------------------------------
    character(len=:), allocatable :: valid_solvers_list !< Pretty printed list of available solvers.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    call self%cli%parse(error=self%error)
    call self%cli%get(switch='-s', val=self%solver, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='--iterations', val=self%implicit_iterations, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-f', val=self%frequency, error=self%error) ; if (self%error/=0) stop
    call self%cli%get_varying(switch='--ss', val=self%stages_steps, error=self%error) ; if (self%error/=0) stop
    call self%cli%get_varying(switch='-Dt', val=self%Dt, error=self%error) ; if (self%error/=0) stop
    call self%cli%get_varying(switch='-tol', val=self%tolerance, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-tf', val=self%final_time, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-r', val=self%results, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-p', val=self%plots, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='--output', val=self%output_cli, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='--errors_analysis', val=self%errors_analysis, error=self%error) ; if (self%error/=0) stop

    if (.not.is_solver_valid()) then
      print "(A)", 'Error: the solver "'//trim(adjustl(self%solver))//'" is unknown!'
      valid_solvers_list = list_solvers()
      print "(A)", valid_solvers_list
      stop
    endif
    if (.not.is_dt_valid()) then
      print "(A)", 'Error: the final integration time must be an exact multiple of the time step used!'
      print "(A)", 'Final integration time: '//str(self%final_time, .true.)
      print "(A)", 'Time step: '//str(self%Dt, .true.)
      stop
    endif
    if (size(self%stages_steps)==2) then
      if (.not.(self%stages_steps(2)>self%stages_steps(1).and.self%stages_steps(1)>=0)) then
        print "(A)", 'Error: when passing a range of stages/steps the valid format must be lower-upper (both positive)!'
        print "(A)", 'Range passed: '//trim(str(self%stages_steps(1), .true.))//'-'//trim(str(self%stages_steps(2), .true.))
        stop
      endif
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine parse_cli

    function is_solver_valid()
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Verify if the selected solver is valid.
    !-------------------------------------------------------------------------------------------------------------------------------
    logical      :: is_solver_valid !< Return true is the selected solver is available.
    integer(I_P) :: s               !< Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    is_solver_valid = .false.
    do s=1, ubound(solvers, dim=1)
      is_solver_valid = (trim(adjustl(self%solver))==trim(adjustl(solvers(s))))
      if (is_solver_valid) exit
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction is_solver_valid

    function is_dt_valid()
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Verify if the selected time step Dt is valid.
    !-------------------------------------------------------------------------------------------------------------------------------
    logical      :: is_dt_valid !< Return true is the selected time step Dt is valid.
    integer(I_P) :: t           !< Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    is_dt_valid = .true.
    do t=1, size(self%Dt)
      is_dt_valid = ((self%final_time - int(self%final_time/self%Dt(t), I_P)*self%Dt(t))==0)
      if (.not.is_dt_valid) exit
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction is_dt_valid

    function list_solvers() result(list)
    !-------------------------------------------------------------------------------------------------------------------------------
    !< List available solvers.
    !-------------------------------------------------------------------------------------------------------------------------------
    character(len=:), allocatable :: list !< Pretty printed list of available solvers.
    integer(I_P)                  :: s    !< Solvers counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    list = 'Valid solver names are:' // new_line('a')
    do s=1, ubound(solvers, dim=1)
      list = list // '  + ' // trim(adjustl(solvers(s))) // new_line('a')
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction list_solvers
  endsubroutine init

  subroutine test(self, solver)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Perform the test.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation_test), intent(in) :: self   !< Test.
  character(*),            intent(in) :: solver !< Selected solver.
  ! FOODIE integrators
  type(integrator_adams_bashforth)         :: ab_integrator     !< Adams-Bashforth integrator.
  type(integrator_adams_bashforth_moulton) :: abm_integrator    !< Adams-Bashforth-Moulton integrator.
  type(integrator_adams_moulton)           :: am_integrator     !< Adams-Moulton integrator.
  type(integrator_back_df)                 :: bdf_integrator    !< BDF integrator.
  type(integrator_runge_kutta_emd)         :: emd_rk_integrator !< Runge-Kutta integrator.
  type(integrator_euler_explicit)          :: euler_integrator  !< Euler integrator.
  type(integrator_leapfrog)                :: lf_integrator     !< Leapfrog integrator.
  type(integrator_runge_kutta_ls)          :: ls_rk_integrator  !< Low Storage Runge-Kutta integrator.
  type(integrator_runge_kutta_tvd)         :: tvd_rk_integrator !< TVD Runge-Kutta integrator.
  ! Auxiliary variables
  real(R_P), allocatable :: solution(:,:)           !< Solution at each time step.
  real(R_P), allocatable :: error(:,:)              !< Error (norm L2) with respect the exact solution.
  real(R_P), allocatable :: order(:,:)              !< Observed order based on subsequent refined solutions.
  real(R_P), allocatable :: Dt_mean(:)              !< Mean time steps used for adaptive solver solutions.
  logical                :: analyze_errors          !< Flag for activating errors analysis for the current tests set.
  integer(I_P)           :: stages_steps_range(1:2) !< Stages/Steps used.
  integer(I_P)           :: last_step               !< Last time step computed.
  integer(I_P)           :: s                       !< Counter.
  integer(I_P)           :: t                       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initialize stages/steps range
  if (size(self%stages_steps)==2) then
    if (self%stages_steps(2)>self%stages_steps(1).and.self%stages_steps(1)>=0) then
      stages_steps_range = [self%stages_steps(1), self%stages_steps(2)]
    endif
  elseif (self%stages_steps(1)>=0) then
    stages_steps_range = self%stages_steps(1)
  else
    select case(trim(adjustl(solver)))
    case('adams-bashforth')
      stages_steps_range = [ab_integrator%min_steps(), ab_integrator%max_steps()]
    case('adams-bashforth-moulton')
      stages_steps_range = [abm_integrator%min_steps(), abm_integrator%max_steps()]
    case('adams-moulton')
      stages_steps_range = [am_integrator%min_steps(), am_integrator%max_steps()]
    case('backward-diff-formula')
      stages_steps_range = [bdf_integrator%min_steps(), bdf_integrator%max_steps()]
    case('emd-runge-kutta')
      stages_steps_range = [emd_rk_integrator%min_stages(), emd_rk_integrator%max_stages()]
    case('euler')
      stages_steps_range = [euler_integrator%min_stages_steps(), euler_integrator%max_stages_steps()]
    case('leapfrog')
      stages_steps_range = [lf_integrator%min_steps(), lf_integrator%max_steps()]
    case('leapfrog-raw')
      stages_steps_range = [lf_integrator%min_steps(), lf_integrator%max_steps()]
    case('ls-runge-kutta')
      stages_steps_range = [ls_rk_integrator%min_stages(), ls_rk_integrator%max_stages()]
    case('tvd-runge-kutta')
      stages_steps_range = [tvd_rk_integrator%min_stages(), tvd_rk_integrator%max_stages()]
    endselect
  endif

  if (trim(adjustl(solver))=='emd-runge-kutta') then
    if (allocated(error)) deallocate(error) ; allocate(error(1:space_dimension, 1:size(self%tolerance)))
    if (allocated(Dt_mean)) deallocate(Dt_mean) ; allocate(Dt_mean(1:size(error, dim=2)))
  else
    if (allocated(error)) deallocate(error) ; allocate(error(1:space_dimension, 1:size(self%Dt)))
  endif
  error = 0.0_R_P

  analyze_errors = .false.
  if (trim(adjustl(solver))=='emd-runge-kutta') then
    if (size(self%tolerance)>1.and.self%errors_analysis) then
      if (allocated(order)) deallocate(order) ; allocate(order(1:space_dimension, 1:size(self%tolerance)-1))
      analyze_errors = .true.
    endif
  else
    if (size(self%Dt)>1.and.self%errors_analysis) then
      if (allocated(order)) deallocate(order) ; allocate(order(1:space_dimension, 1:size(self%Dt)-1))
      analyze_errors = .true.
    endif
  endif

  ! test(s)
  print "(A)", trim(adjustl(solver))
  do s=stages_steps_range(1), stages_steps_range(2)
    print "(A)", '  stages/steps '//trim(str(s,.true.))
    if (trim(adjustl(solver))=='emd-runge-kutta') then
      do t=1, size(self%tolerance)
        call solve(solver=solver,                       &
                   frequency=self%frequency,            &
                   final_time=self%final_time,          &
                   stages_steps=s,                      &
                   iterations=self%implicit_iterations, &
                   solution=solution,                   &
                   error=error(:, t),                   &
                   last_step=last_step,                 &
                   tolerance=self%tolerance(t))
        Dt_mean(t) = self%final_time/real(last_step, kind=R_P)

        if (allocated(solution)) then
          print "(A,I10,A,F10.3,A,F10.3,A,E10.3,A,E10.3)", "    steps: ", last_step,              &
                                                           "    Dt: ", Dt_mean(t),                &
                                                           ", f*Dt: ", self%frequency*Dt_mean(t), &
                                                           ", E(x): ", error(1, t), ", E(y): ", error(2, t)
          if (analyze_errors.and.t>1) then
            order(:, t-1) = observed_order(error=error(:, t-1:t), Dt=Dt_mean(t-1:t))
            print "(A,F10.2,A,F10.2)", "      Observed order, O(x): ", order(1, t-1), ", O(y): " , order(2, t-1)
          endif
          call save_results(results=self%results,                                    &
                            plots=self%plots,                                        &
                            output_cli=self%output_cli,                              &
                            solver=trim(adjustl(solver))//'-'//trim(str(s, .true.)), &
                            frequency=self%frequency,                                &
                            solution=solution(:, 0:last_step))
        endif
      enddo
    else
      do t=1, size(self%Dt)
        call solve(solver=solver,                       &
                   frequency=self%frequency,            &
                   final_time=self%final_time,          &
                   stages_steps=s,                      &
                   iterations=self%implicit_iterations, &
                   solution=solution,                   &
                   error=error(:, t),                   &
                   last_step=last_step,                 &
                   Dt=self%Dt(t))

        if (allocated(solution)) then
          print "(A,I10,A,F10.3,A,F10.3,A,E10.3,A,E10.3)", "    steps: ", last_step,              &
                                                           "    Dt: ", self%Dt(t), &
                                                           ", f*Dt: ", self%frequency*self%Dt(t), &
                                                           ", E(x): ", error(1, t), ", E(y): ", error(2, t)
          if (analyze_errors.and.t>1) then
            order(:, t-1) = observed_order(error=error(:, t-1:t), Dt=self%Dt(t-1:t))
            print "(A,F10.2,A,F10.2)", "      Observed order, O(x): ", order(1, t-1), ", O(y): " , order(2, t-1)
          endif
          call save_results(results=self%results,                                    &
                            plots=self%plots,                                        &
                            output_cli=self%output_cli,                              &
                            solver=trim(adjustl(solver))//'-'//trim(str(s, .true.)), &
                            frequency=self%frequency,                                &
                            solution=solution(:, 0:last_step))
        endif
      enddo
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test

  ! non type bound procedures
  subroutine solve(solver, frequency, final_time, stages_steps, iterations, solution, error, last_step, Dt, tolerance)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Rune the solver selected.
  !<
  !< The actual solver is selected by means of the *solver* input string that must be a valid string as defined into *solvers*
  !< parameter list.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*),           intent(in)  :: solver        !< Selected solver: must be defined into *solvers*.
  real(R_P),              intent(in)  :: frequency     !< Oscillation frequency.
  real(R_P),              intent(in)  :: final_time    !< Final integration time.
  integer(I_P),           intent(in)  :: stages_steps  !< Number of stages/steps used.
  integer(I_P),           intent(in)  :: iterations    !< Number of fixed point iterations.
  real(R_P), allocatable, intent(out) :: solution(:,:) !< Solution at each time step, X-Y.
  real(R_P),              intent(out) :: error(1:)     !< Error (norm L2) with respect the exact solution.
  integer(I_P),           intent(out) :: last_step     !< Last time step computed.
  real(R_P), optional,    intent(in)  :: Dt            !< Time step.
  real(R_P), optional,    intent(in)  :: tolerance     !< Local error tolerance.
  ! FOODIE integrators
  type(integrator_adams_bashforth)         :: ab_integrator     !< Adams-Bashforth integrator.
  type(integrator_adams_bashforth_moulton) :: abm_integrator    !< Adams-Bashforth-Moulton integrator.
  type(integrator_adams_moulton)           :: am_integrator     !< Adams-Moulton integrator.
  type(integrator_back_df)                 :: bdf_integrator    !< BDF integrator.
  type(integrator_runge_kutta_emd)         :: emd_rk_integrator !< Runge-Kutta integrator.
  type(integrator_euler_explicit)          :: euler_integrator  !< Euler integrator.
  type(integrator_leapfrog)                :: lf_integrator     !< Leapfrog integrator.
  type(integrator_runge_kutta_ls)          :: ls_rk_integrator  !< Low Storage Runge-Kutta integrator.
  type(integrator_runge_kutta_tvd)         :: tvd_rk_integrator !< TVD Runge-Kutta integrator.
  ! Auxiliary variables
  integer(I_P), parameter        :: max_rk_stages=5 !< Max RK stages used to init high order multi-step solver.
  type(oscillation)              :: oscillator      !< Oscillation field.
  logical                        :: supported       !< Flag for checking if the selected solver is supported.
  logical                        :: multistep       !< Flag for tagging multi-step class of solvers.
  integer                        :: step_offset     !< Time steps counter offset for slicing previous data array.
  logical                        :: adaptive        !< Flag for tagging time step adaptive class of solvers.
  real(R_P)                      :: Dt_a            !< Adaptive time step.
  type(oscillation), allocatable :: rk_stage(:)     !< Runge-Kutta stages.
  type(oscillation), allocatable :: previous(:)     !< Previous time steps solutions.
  type(oscillation)              :: filter          !< Filter displacement.
  integer                        :: step            !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initizialize
  call oscillator%init(initial_state=initial_state, frequency=frequency)

  supported = .true.           ! assume that selected solver is supported
  multistep = .false.          ! assume multi-stage solver by default
  if (stages_steps==0) then
    step_offset = 1            ! for 0 step-(a convention)-solver offset is 1
  else
    step_offset = stages_steps ! for >0 step-solver offset is steps
  endif
  adaptive = .false.           ! assume fixed time step by default
  Dt_a = 10000._R_P            ! initial step very large to trigger adaptation

  if (present(Dt)) then
    if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/Dt)))
  else
    if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(final_time/10))) ! hope is enough
  endif
  solution = 0.0_R_P
  solution(1:, 0) = oscillator%output()

  select case(trim(adjustl(solver)))
  case("adams-bashforth")
    supported = ab_integrator%is_supported(stages_steps)
    multistep = .true.
    call ab_integrator%init(steps=stages_steps)
    if (allocated(previous)) deallocate(previous) ; allocate(previous(1:stages_steps))

  case("adams-bashforth-moulton")
    supported = abm_integrator%is_supported(stages_steps)
    multistep = .true.
    call abm_integrator%init(steps=stages_steps)
    if (allocated(previous)) deallocate(previous) ; allocate(previous(1:stages_steps))

  case("adams-moulton")
    supported = am_integrator%is_supported(stages_steps)
    multistep = .true.
    call am_integrator%init(steps=stages_steps)
    if (allocated(previous)) deallocate(previous) ; allocate(previous(1:stages_steps+1))

  case("backward-diff-formula")
    supported = bdf_integrator%is_supported(stages_steps)
    multistep = .true.
    call bdf_integrator%init(steps=stages_steps)
    if (allocated(previous)) deallocate(previous) ; allocate(previous(1:stages_steps+1))

  case("emd-runge-kutta")
    supported = emd_rk_integrator%is_supported(stages_steps)
    adaptive = .true.
    call emd_rk_integrator%init(stages=stages_steps, tolerance=tolerance)
    if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:stages_steps))

  case("euler")
    supported = euler_integrator%is_supported(stages_steps)

  case("leapfrog")
    supported = lf_integrator%is_supported(stages_steps)
    multistep = .true.
    call lf_integrator%init()
    if (allocated(previous)) deallocate(previous) ; allocate(previous(1:stages_steps))

  case("leapfrog-raw")
    supported = lf_integrator%is_supported(stages_steps)
    multistep = .true.
    call lf_integrator%init()
    if (allocated(previous)) deallocate(previous) ; allocate(previous(1:stages_steps))

  case("ls-runge-kutta")
    supported = ls_rk_integrator%is_supported(stages_steps)
    call ls_rk_integrator%init(stages=stages_steps)
    if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:ls_rk_integrator%used_registers()))

  case("tvd-runge-kutta")
    supported = tvd_rk_integrator%is_supported(stages_steps)
    call tvd_rk_integrator%init(stages=stages_steps)
    if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:stages_steps))
  endselect

  if (multistep) then
    call tvd_rk_integrator%init(stages=max_rk_stages)
    if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:max_rk_stages))
  endif

  if (.not.supported) then
    print "(A)", 'The solver '//trim(adjustl(solver))//' does not support '//trim(str(stages_steps, .true.))//' stages/steps'
    if (allocated(solution)) deallocate(solution)
    return
  endif

  ! integrate
  step = 0
  do while(solution(0, step)<final_time.and.step<ubound(solution, dim=2))
    step = step + 1

    ! advance oscillator
    if (multistep) then
      if (stages_steps>=step) then
        call tvd_rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
        previous(step) = oscillator
      else
        select case(trim(adjustl(solver)))
        case("adams-bashforth")
          call ab_integrator%integrate(U=oscillator,      &
                                       previous=previous, &
                                       Dt=Dt,             &
                                       t=solution(0, step-step_offset:step-1))

        case("adams-bashforth-moulton")
          call abm_integrator%integrate(U=oscillator,      &
                                        previous=previous, &
                                        Dt=Dt,             &
                                        t=solution(0, step-step_offset:step-1))

        case("adams-moulton")
          if (iterations>1) then
            call am_integrator%integrate(U=oscillator,                          &
                                         previous=previous,                     &
                                         Dt=Dt,                                 &
                                         t=solution(0,step-step_offset:step-1), &
                                         iterations=iterations)
          else
            call am_integrator%integrate(U=oscillator,                          &
                                         previous=previous,                     &
                                         Dt=Dt,                                 &
                                         t=solution(0,step-step_offset:step-1))
          endif

        case("backward-diff-formula")
          if (iterations>1) then
            call bdf_integrator%integrate(U=oscillator,                          &
                                          previous=previous,                     &
                                          Dt=Dt,                                 &
                                          t=solution(0,step-step_offset:step-1), &
                                          iterations=iterations)
          else
            call bdf_integrator%integrate(U=oscillator,                          &
                                          previous=previous,                     &
                                          Dt=Dt,                                 &
                                          t=solution(0,step-step_offset:step-1))
          endif

        case("leapfrog")
          call lf_integrator%integrate(U=oscillator,      &
                                       previous=previous, &
                                       Dt=Dt,             &
                                       t=solution(0, step))

        case("leapfrog-raw")
          call lf_integrator%integrate(U=oscillator,        &
                                       previous=previous,   &
                                       Dt=Dt,               &
                                       t=solution(0, step), &
                                       filter=filter)
        endselect
      endif
    else
      select case(trim(adjustl(solver)))
      case("emd-runge-kutta")
        call emd_rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt_a, t=solution(0, step))

      case("euler")
        call euler_integrator%integrate(U=oscillator, Dt=Dt, t=solution(0, step))

      case("ls-runge-kutta")
        call ls_rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))

      case("tvd-runge-kutta")
        call tvd_rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))

      endselect
    endif

    ! update time
    if (adaptive) then
      solution(0, step) = solution(0, step - 1) + Dt_a
    else
      solution(0, step) = step * Dt
    endif

    ! store step solution
    solution(1:, step) = oscillator%output()
  enddo
  last_step = step

  error = error_L2(frequency=frequency, solution=solution(:, 0:last_step))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine solve

  subroutine save_results(results, plots, output_cli, solver, frequency, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save results (and plots).
  !---------------------------------------------------------------------------------------------------------------------------------
  logical,      intent(in)      :: results                   !< Flag for activating results saving.
  logical,      intent(in)      :: plots                     !< Flag for activating plots saving.
  character(*), intent(in)      :: output_cli                !< Output files basename coming from CLI.
  character(*), intent(in)      :: solver                    !< Selected solver: must be defined into *solvers*.
  real(R_P),    intent(in)      :: frequency                 !< Oscillation frequency.
  real(R_P),    intent(IN)      :: solution(0:,0:)           !< Solution at each time step.
  character(len=:), allocatable :: title                     !< Output files title.
  character(len=:), allocatable :: basename                  !< Output files basename.
  integer(I_P)                  :: rawfile                   !< Raw file unit for saving results.
  real(R_P)                     :: ex_sol(1:space_dimension) !< Exact solution.
  type(pyplot)                  :: plt                       !< Plot file handler.
  integer(I_P)                  :: s                         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  title = 'Oscillation equations integration, solver='//trim(adjustl(solver))
  if (trim(adjustl(output_cli))/='unset') then
    basename = trim(adjustl(output_cli))//'-'//trim(strz(ubound(solution, dim=2), 10))//'-time_steps-'//trim(adjustl(solver))
  else
    basename = 'oscillation_test-'//trim(strz(ubound(solution, dim=2), 10))//'-time_steps-'//trim(adjustl(solver))
  endif
  if (results) then
    open(newunit=rawfile, file=basename//'.dat')
    write(rawfile, '(A)')'TITLE="'//title//'"'
    write(rawfile, '(A)')'VARIABLES="t" "x" "y" "amplitude" "phase"'
    write(rawfile, '(A)')'ZONE T="FOODIE time serie"'
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
  if (plots) then
    call plt%initialize(grid=.true., xlabel='time', title=title, legend=.true.)
    call plt%add_plot(x=solution(0, :), y=solution(1, :), label='u', linestyle='r-', linewidth=1)
    call plt%add_plot(x=solution(0, :), y=solution(2, :), label='v', linestyle='b-', linewidth=1)
    call plt%savefig(basename//'.png')

    call plt%initialize(grid=.true., xlabel='v1', ylabel='v2', title=title//'-path', legend=.true.)
    call plt%add_plot(x=solution(1, :), y=solution(2, :), label='path', linestyle='r-', linewidth=1)
    call plt%savefig('path-'//basename//'.png')
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    function amplitude_phase(sol) result(ap)
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Compute amplitude and phase of the solution provided in X-Y domain.
    !-------------------------------------------------------------------------------------------------------------------------------
    real(R_P), intent(IN) :: sol(1:) !< Solution in X-Y domain.
    real(R_P)             :: ap(1:2) !< Amplitude and phase solution.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    ap(1) = sqrt(sol(1)**2 + sol(2)**2)
    ap(2) = atan(-sol(1) / sol(2))
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction amplitude_phase
  endsubroutine save_results

  pure function exact_solution_xy(frequency, t) result(ex_sol)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the exact solution on X-Y domain.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(in) :: frequency   !< Oscillation frequency.
  real(R_P), intent(IN) :: t           !< Time.
  real(R_P)             :: ex_sol(1:2) !< Exact solution.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ex_sol(1) = initial_state(1) * cos(frequency * t) - initial_state(2) * sin(frequency * t)
  ex_sol(2) = initial_state(1) * sin(frequency * t) + initial_state(2) * cos(frequency * t)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction exact_solution_xy

  pure function error_L2(frequency, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the L2 norm of numerical error with respect the exact solution.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(in) :: frequency                   !< Oscillation frequency.
  real(R_P), intent(IN) :: solution(0:,0:)             !< Solution at each time step.
  real(R_P)             :: error_L2(1:space_dimension) !< L2 norm of the numerical error.
  real(R_P)             :: ex_sol(1:space_dimension)   !< Exact solution.
  integer(I_P)          :: s                           !< Steps/stages counter.
  integer(I_P)          :: v                           !< Variables counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  error_L2 = 0._R_P
  do s=0, ubound(solution, dim=2)
    ex_sol = exact_solution_xy(frequency=frequency, t=solution(0, s))
    do v=1, space_dimension
      error_L2(v) = error_L2(v) + (solution(v, s) - ex_sol(v))*(solution(v, s) - ex_sol(v))
    enddo
  enddo
  error_L2 = sqrt(error_L2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction error_L2

  pure function observed_order(error, Dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Estimate the order of accuracy using 2 subsequent refined numerical solutions.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN) :: error(1:space_dimension, 1:2)     !< Computed errors.
  real(R_P), intent(IN) :: Dt(1:2)                           !< Time steps used.
  real(R_P)             :: observed_order(1:space_dimension) !< Estimation of the order of accuracy.
  integer(I_P)          :: v                                 !< Variables counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do v=1, space_dimension
    observed_order(v) = log(error(v, 1) / error(v, 2)) / log(Dt(1) / Dt(2))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction observed_order
endmodule oscillation_test_t

program integrate_oscillation
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODIE with the integration of Oscillation equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use oscillation_test_t, only : oscillation_test
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(oscillation_test) :: test !< Oscillation test.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
call test%execute
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram integrate_oscillation
