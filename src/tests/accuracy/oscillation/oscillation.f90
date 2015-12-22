!< Test FOODIE with the integration of Oscillation equations.
program integrate_oscillation
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODIE with the integration of Oscillation equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, FR_P, str, strz
use type_oscillation, only : oscillation
use Data_Type_Command_Line_Interface, only : Type_Command_Line_Interface
use foodie, only : adams_bashforth_integrator, &
                   adams_bashforth_moulton_integrator, &
                   adams_moulton_integrator, &
                   bdf_integrator, &
                   emd_runge_kutta_integrator, &
                   euler_explicit_integrator, &
                   leapfrog_integrator, &
                   ls_runge_kutta_integrator, &
                   tvd_runge_kutta_integrator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Command_Line_Interface) :: cli                                               !< Command line interface handler.
integer,   parameter              :: space_dimension=2                                 !< Space dimensions.
real(R_P), parameter              :: initial_state(1:space_dimension)=[0._R_P, 1._R_P] !< Initial state.
real(R_P)                         :: f                                                 !< Oscillation frequency.
real(R_P)                         :: t_final                                           !< Final time.
real(R_P)                         :: Dt                                                !< Time step.
real(R_P)                         :: tolerance                                         !< Tolerance on local truncation error.
integer(I_P)                      :: error                                             !< Error handler.
integer(I_P)                      :: stages_steps                                      !< Number of stages/steps used.
integer(I_P)                      :: implicit_iterations                               !< Number of iterations for implicit solvers.
character(99)                     :: solver                                            !< Solver used.
logical                           :: plots                                             !< Flag for activating plots saving.
logical                           :: results                                           !< Flag for activating results saving.
character(99)                     :: output_cli                                        !< Output files basename.
logical                           :: errors_analysis                                   !< Flag for activating errors analysis.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'oscillation',                                              &
              authors     = 'Fortran-FOSS-Programmers',                                 &
              license     = 'GNU GPLv3',                                                &
              description = 'Test FOODIE library on Oscillation equations integration', &
              examples    = ["oscillation --solver euler --results  ",                  &
                             "oscillation --solver ls-runge-kutta -r",                  &
                             "oscillation --solver adams-bashforth  ",                  &
                             "oscillation --solver all --plots -r   "])
call cli%add(switch='--solver', switch_ab='-s', help='ODE solver used', required=.false., def='all', act='store')
call cli%add(switch='--ss', help='Stages/steps used', required=.false., act='store', def='-1')
call cli%add(switch='--iterations', help='Number of iterations for implicit solvers', required=.false., act='store', def='5')
call cli%add(switch='--frequency', switch_ab='-f', help='Oscillation frequency', required=.false., def='1e-4', act='store')
call cli%add(switch='--time_step', switch_ab='-Dt', help='Integration time step', required=.false., def='100.d0', act='store')
call cli%add(switch='--tolerance', switch_ab='-tol', help='Tolerance on local error', required=.false., def='0.001d0', act='store')
call cli%add(switch='--t_final', switch_ab='-tf', help='Final integration time', required=.false., def='1e6', act='store')
call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.')
call cli%add(switch='--plots', switch_ab='-p', help='Save plots of results', required=.false., act='store_true', def='.false.')
call cli%add(switch='--output', help='Output files basename', required=.false., act='store', def='unset')
call cli%add(switch='--errors_analysis', help='Peform errors analysis', required=.false., act='store_true', def='.false.')
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='-s', val=solver, error=error) ; if (error/=0) stop
call cli%get(switch='--ss', val=stages_steps, error=error) ; if (error/=0) stop
call cli%get(switch='--iterations', val=implicit_iterations, error=error) ; if (error/=0) stop
call cli%get(switch='-f', val=f, error=error) ; if (error/=0) stop
call cli%get(switch='-Dt', val=Dt, error=error) ; if (error/=0) stop
call cli%get(switch='-tol', val=tolerance, error=error) ; if (error/=0) stop
call cli%get(switch='-tf', val=t_final, error=error) ; if (error/=0) stop
call cli%get(switch='-r', val=results, error=error) ; if (error/=0) stop
call cli%get(switch='-p', val=plots, error=error) ; if (error/=0) stop
call cli%get(switch='--output', val=output_cli, error=error) ; if (error/=0) stop
call cli%get(switch='--errors_analysis', val=errors_analysis, error=error) ; if (error/=0) stop
if (errors_analysis) then
  call perform_errors_analysis
else
  select case(trim(adjustl(solver)))
  case('adams-bashforth')
    call test_ab(steps=stages_steps)
  case('adams-bashforth-moulton')
    call test_abm(steps=stages_steps)
  case('adams-moulton')
    call test_am(steps=stages_steps)
  case('backward-diff-formula')
    call test_bdf(steps=stages_steps)
  case('emd-runge-kutta')
    call test_emd_rk(stages=stages_steps)
  case('euler')
    call test_euler
  case('leapfrog')
    call test_leapfrog
  case('leapfrog-raw')
    call test_leapfrog(raw=.true.)
  case('ls-runge-kutta')
    call test_ls_rk(stages=stages_steps)
  case('tvd-runge-kutta')
    call test_tvd_rk(stages=stages_steps)
  case('all')
    call test_ab(steps=stages_steps)
    call test_abm(steps=stages_steps)
    call test_am(steps=stages_steps)
    call test_emd_rk(stages=stages_steps)
    call test_euler
    call test_leapfrog
    call test_ls_rk(stages=stages_steps)
    call test_tvd_rk(stages=stages_steps)
  case default
    print "(A)", 'Error: unknown solver "'//trim(adjustl(solver))//'"'
    print "(A)", 'Valid solver names are:'
    print "(A)", '  + adams-bashforth'
    print "(A)", '  + adams-bashforth-moulton'
    print "(A)", '  + adams-moulton'
    print "(A)", '  + backward-diff-formula'
    print "(A)", '  + emd-runge-kutta'
    print "(A)", '  + euler'
    print "(A)", '  + leapfrog'
    print "(A)", '  + leapfrog-raw'
    print "(A)", '  + ls-runge-kutta'
    print "(A)", '  + tvd-runge-kutta'
    print "(A)", '  + all'
  endselect
endif
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine init(output, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Initialize solver local variable from global variables.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(len=:), allocatable, intent(OUT) :: output        !< Output files basename.
  real(R_P),        allocatable, intent(OUT) :: solution(:,:) !< Solution at each time step.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (mod(t_final, Dt)/=0) then
    print "(A)", 'Error: the final integration time must be an exact multiple of the time step used!'
    print "(A)", 'Final integration time: '//str(.true.,t_final)
    print "(A)", 'Time step: '//str(.true.,Dt)
    stop
  endif
  if (trim(adjustl(output_cli))/='unset') then
    if (trim(adjustl(solver))/='emd-runge-kutta') then
      output = trim(adjustl(output_cli))//'-'//trim(strz(10,int(t_final/Dt)))//'-time_steps-'//trim(adjustl(solver))
    else
      output = trim(adjustl(output_cli))//'-'//trim(str(.true.,n=tolerance))//'-tolerance-'//trim(adjustl(solver))
    endif
  else
    if (trim(adjustl(solver))/='emd-runge-kutta') then
      output = 'oscillation_integration-'//trim(strz(10,int(t_final/Dt)))//'-time_steps-'//trim(adjustl(solver))
    else
      output = 'oscillation_integration-'//trim(str(.true.,n=tolerance))//'-tolerance-'//trim(adjustl(solver))
    endif
  endif
  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(t_final/Dt)))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  function exact_solution(t) result(ex_sol)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the exact solution.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN) :: t                         !< Time.
  real(R_P)             :: ex_sol(1:space_dimension) !< Exact solution.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ex_sol(1) = initial_state(1) * cos(f * t) - initial_state(2) * sin(f * t)
  ex_sol(2) = initial_state(1) * sin(f * t) + initial_state(2) * cos(f * t)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction exact_solution

  subroutine save_results(title, basename, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save plots of results.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN) :: title                     !< Plot title.
  character(*), intent(IN) :: basename                  !< Output files basename.
  real(R_P),    intent(IN) :: solution(0:,0:)           !< Solution at each time step.
  integer(I_P)             :: rawfile                   !< Raw file unit for saving results.
  real(R_P)                :: ex_sol(1:space_dimension) !< Exact solution.
  type(pyplot)             :: plt                       !< Plot file handler.
  integer(I_P)             :: i                         !< Counter.
  integer(I_P)             :: s                         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (results) then
    open(newunit=rawfile, file=basename//'.dat')
    write(rawfile, '(A)')'TITLE="'//title//'"'
    write(rawfile, '(A)')'VARIABLES="t" "x" "y"'
    write(rawfile, '(A)')'ZONE T="FOODIE time serie"'
    do s=0, ubound(solution, dim=2)
      write(rawfile, '(3('//FR_P//',1X))')(solution(i, s), i=0, space_dimension)
    enddo
    write(rawfile, '(A)')'ZONE T="Exact solution"'
    do s=0, ubound(solution, dim=2)
      ex_sol = exact_solution(t=solution(0, s))
      write(rawfile, '(3('//FR_P//',1X))')solution(0, s), (ex_sol(i), i=1, space_dimension)
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
  endsubroutine save_results

  subroutine test_ab(steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit Adams-Bashforth class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)      :: steps            !< Number of steps used: if negative all AB solvers are used.
  integer, parameter            :: ab_steps=4       !< Adams-Bashforth steps number.
  real(R_P), allocatable        :: solution(:,:)    !< Solution at each time step.
  integer(I_P)                  :: s                !< AB steps counter.
  integer(I_P)                  :: steps_range(1:2) !< Steps used.
  character(len=:), allocatable :: output           !< Output files basename.
  character(len=:), allocatable :: title            !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  solver = 'adams-bashforth'
  call init(output=output, solution=solution)
  print "(A)", 'Integrating Oscillation equations by means of Adams-Bashforth class of solvers'
  steps_range = [1, ab_steps] ; if (steps>0) steps_range = [steps, steps]
  do s=steps_range(1), steps_range(2)
    print "(A)", ' AB-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, explicit Adams-Bashforth, t='//str(n=t_final)//' steps='//trim(str(.true., s))
    call ab_solver(steps=s, solution=solution)
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_ab

  subroutine test_abm(steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test predictor-corrector Adams-Bashforth-Moulton class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)      :: steps            !< Number of steps used: if negative all AB solvers are used.
  integer, parameter            :: abm_steps=4      !< Adams-Bashforth-Moulton steps number.
  real(R_P), allocatable        :: solution(:,:)    !< Solution at each time step.
  integer(I_P)                  :: s                !< ABM steps counter.
  integer(I_P)                  :: steps_range(1:2) !< Steps used.
  character(len=:), allocatable :: output           !< Output files basename.
  character(len=:), allocatable :: title            !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  solver = 'adams-bashforth-moulton'
  call init(output=output, solution=solution)
  print "(A)", 'Integrating Oscillation equations by means of Adams-Bashforth-Moulton class of solvers'
  steps_range = [1, abm_steps] ; if (steps>0) steps_range = [steps, steps]
  do s=steps_range(1), steps_range(2)
    print "(A)", ' ABM-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, predictor/corrector Adams-Bashforth-Moulton, t='//str(n=t_final)//&
            ' steps='//trim(str(.true., s))
    call abm_solver(steps=s, solution=solution)
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_abm

  subroutine test_am(steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test implicit Adams-Moulton class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)      :: steps            !< Number of steps used: if negative all AB solvers are used.
  integer, parameter            :: am_steps=3       !< Adams-Moulton steps number.
  real(R_P), allocatable        :: solution(:,:)    !< Solution at each time step.
  integer(I_P)                  :: s                !< AM steps counter.
  integer(I_P)                  :: steps_range(1:2) !< Steps used.
  character(len=:), allocatable :: output           !< Output files basename.
  character(len=:), allocatable :: title            !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  solver = 'adams-moulton'
  call init(output=output, solution=solution)
  print "(A)", 'Integrating Oscillation equations by means of Adams-Moulton class of solvers'
  steps_range = [0, am_steps] ; if (steps>0) steps_range = [steps, steps]
  do s=steps_range(1), steps_range(2)
    print "(A)", ' AM-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, implicit Adams-Moulton, t='//str(n=t_final)//&
            ' steps='//trim(str(.true., s))
    if (implicit_iterations>1) then
      call am_solver(steps=s, solution=solution, iterations=implicit_iterations)
    else
      call am_solver(steps=s, solution=solution)
    endif
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_am

  subroutine test_bdf(steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test implicit Backward-Differentiation-Formula class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)      :: steps            !< Number of steps used: if negative all BDF solvers are used.
  integer, parameter            :: bdf_steps=3      !< BDF steps number.
  real(R_P), allocatable        :: solution(:,:)    !< Solution at each time step.
  integer(I_P)                  :: s                !< Steps counter.
  integer(I_P)                  :: steps_range(1:2) !< Steps used.
  character(len=:), allocatable :: output           !< Output files basename.
  character(len=:), allocatable :: title            !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  solver = 'bdf'
  call init(output=output, solution=solution)
  print "(A)", 'Integrating Oscillation equations by means of Backward-Differentiation-Formula class of solvers'
  steps_range = [0, bdf_steps] ; if (steps>0) steps_range = [steps, steps]
  do s=steps_range(1), steps_range(2)
    print "(A)", ' BDF-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, implicit Backward-Differentation-Formula, t='//str(n=t_final)//&
            ' steps='//trim(str(.true., s))
    if (implicit_iterations>1) then
      call bdf_solver(steps=s, solution=solution, iterations=implicit_iterations)
    else
      call bdf_solver(steps=s, solution=solution)
    endif
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_bdf

  subroutine test_emd_rk(stages)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit embedded Runge-Kutta class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)      :: stages            !< Number of stages used: if negative all RK solvers are used.
  integer, parameter            :: rk_stages=17      !< Runge-Kutta stages number.
  real(R_P), allocatable        :: solution(:,:)     !< Solution at each time step.
  integer(I_P)                  :: s                 !< RK stages counter.
  integer(I_P)                  :: stages_range(1:2) !< Stages used.
  character(len=:), allocatable :: output            !< Output files basename.
  character(len=:), allocatable :: title             !< Output files title.
  integer(I_P)                  :: last_step         !< Last time step computed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  solver = 'emd-runge-kutta'
  call init(output=output, solution=solution)
  print "(A)", 'Integrating Oscillation equations by means of embedded Runge-Kutta class of solvers'
  stages_range = [1, rk_stages] ; if (stages>0) stages_range = [stages, stages]
  do s=stages_range(1), stages_range(2)
    if (s==1) cycle ! 1 stages not allowed
    if (s==2) cycle ! 2 stages not suited for test
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    if (s==5) cycle ! 5 stages not yet implemented
    if (s==8) cycle ! 8 stages not yet implemented
    if (s==10) cycle ! 10 stages not yet implemented
    if (s==11) cycle ! 11 stages not yet implemented
    if (s==12) cycle ! 12 stages not yet implemented
    if (s==13) cycle ! 13 stages not yet implemented
    if (s==14) cycle ! 14 stages not yet implemented
    if (s==15) cycle ! 15 stages not yet implemented
    if (s==16) cycle ! 16 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, explicit embedded Runge-Kutta, t='//str(n=t_final)//' steps='//trim(str(.true., s))
    call emd_rk_solver(stages=s, tol=tolerance, solution=solution, last_step=last_step)
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution(:, 0:last_step))
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_emd_rk

  subroutine test_euler()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit forward Euler ODE solver.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), allocatable        :: solution(:,:) !< Solution at each time step.
  character(len=:), allocatable :: output        !< Output files basename.
  character(len=:), allocatable :: title         !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  solver = 'euler'
  call init(output=output, solution=solution)
  print "(A)", 'Integrating Oscillation equations by means of explicit Euler solver'
  title = 'Oscillation equation integration, explicit Euler, t='//str(n=t_final)
  call euler_solver(solution=solution)
  call save_results(title=title, basename=output, solution=solution)
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_euler

  subroutine test_leapfrog(raw)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit leapfrog class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  logical, optional, intent(IN) :: raw           !< Activate RAW filter.
  real(R_P), allocatable        :: solution(:,:) !< Solution at each time step.
  character(len=:), allocatable :: output        !< Output files basename.
  character(len=:), allocatable :: title         !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  solver = 'leapfrog'
  call init(output=output, solution=solution)
  print "(A)", 'Integrating Oscillation equations by means of leapfrog (RAW filtered) class of solvers'
  title = 'Oscillation equation integration, leapfrog (RAW filtered), t='//str(n=t_final)
  call leapfrog_solver(solution=solution, raw=raw)
  call save_results(title=title, basename=output, solution=solution)
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_leapfrog

  subroutine test_ls_rk(stages)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit low storage Runge-Kutta class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)      :: stages            !< Number of stages used: if negative all RK solvers are used.
  integer, parameter            :: rk_stages=14      !< Runge-Kutta stages number.
  real(R_P), allocatable        :: solution(:,:)     !< Solution at each time step.
  integer(I_P)                  :: s                 !< RK stages counter.
  integer(I_P)                  :: stages_range(1:2) !< Stages used.
  character(len=:), allocatable :: output            !< Output files basename.
  character(len=:), allocatable :: title             !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  solver = 'ls-runge-kutta'
  call init(output=output, solution=solution)
  print "(A)", 'Integrating Oscillation equations by means of low storage (2N) Runge-Kutta class of solvers'
  stages_range = [1, rk_stages] ; if (stages>0) stages_range = [stages, stages]
  do s=stages_range(1), stages_range(2)
    if (s==2) cycle ! 2 stages not yet implemented
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    if (s==8) cycle ! 8 stages not yet implemented
    if (s==9) cycle ! 9 stages not yet implemented
    if (s==10) cycle ! 10 stages not yet implemented
    if (s==11) cycle ! 11 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, explicit low storage Runge-Kutta, t='//str(n=t_final)//&
      ' steps='//trim(str(.true., s))
    call ls_rk_solver(stages=s, solution=solution)
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_ls_rk

  subroutine test_tvd_rk(stages)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit TVD/SSP Runge-Kutta class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)      :: stages            !< Number of stages used: if negative all RK solvers are used.
  integer, parameter            :: rk_stages=5       !< Runge-Kutta stages number.
  real(R_P), allocatable        :: solution(:,:)     !< Solution at each time step.
  integer(I_P)                  :: s                 !< RK stages counter.
  integer(I_P)                  :: stages_range(1:2) !< Stages used.
  character(len=:), allocatable :: output            !< Output files basename.
  character(len=:), allocatable :: title             !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  solver = 'tvd-runge-kutta'
  call init(output=output, solution=solution)
  print "(A)", 'Integrating Oscillation equations by means of TVD/SSP Runge-Kutta class of solvers'
  stages_range = [1, rk_stages] ; if (stages>0) stages_range = [stages, stages]
  do s=stages_range(1), stages_range(2)
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, explicit TVD/SSP Runge-Kutta, t='//str(n=t_final)//' steps='//trim(str(.true., s))
    call tvd_rk_solver(stages=s, solution=solution)
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_tvd_rk

  subroutine ab_solver(steps, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve problem with explicit Adams-Bashforth class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)         :: steps                 !< Number of steps used: if negative all AB solvers are used.
  real(R_P),    intent(OUT)        :: solution(0:,0:)       !< Solution at each time step.
  type(tvd_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter               :: rk_stages=5           !< Runge-Kutta stages number.
  type(oscillation)                :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(adams_bashforth_integrator) :: ab_integrator         !< Adams-Bashforth integrator.
  type(oscillation)                :: previous(1:steps)     !< Previous time steps solutions.
  type(oscillation)                :: oscillator            !< Oscillation field.
  integer                          :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call ab_integrator%init(steps=steps)
  select case(steps)
  case(1, 2, 3)
    call rk_integrator%init(stages=steps)
  case(4)
    call rk_integrator%init(stages=5)
  endselect
  call oscillator%init(initial_state=initial_state, f=f)
  solution = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    if (steps>=step) then
      ! time steps from 1 to s - 1 must be computed with other scheme
      call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = oscillator
    else
      call ab_integrator%integrate(U=oscillator, previous=previous, Dt=Dt, t=solution(0, step-steps:step-1))
    endif
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine ab_solver

  subroutine abm_solver(steps, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve problem with predictor-corrector Adams-Bashforth-Moulton class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)                 :: steps                 !< Number of steps used: if negative all ABM solvers are used.
  real(R_P),    intent(OUT)                :: solution(0:,0:)       !< Solution at each time step.
  type(tvd_runge_kutta_integrator)         :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter                       :: rk_stages=5           !< Runge-Kutta stages number.
  type(oscillation)                        :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(adams_bashforth_moulton_integrator) :: abm_integrator        !< Adams-Bashforth-Moulton integrator.
  type(oscillation)                        :: previous(1:steps)     !< Previous time steps solutions.
  type(oscillation)                        :: oscillator            !< Oscillation field.
  integer                                  :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call abm_integrator%init(steps=steps)
  select case(steps)
  case(1, 2, 3)
    call rk_integrator%init(stages=steps)
  case(4)
    call rk_integrator%init(stages=5)
  endselect
  call oscillator%init(initial_state=initial_state, f=f)
  solution = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    if (steps>=step) then
      ! time steps from 1 to s - 1 must be computed with other scheme
      call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = oscillator
    else
      call abm_integrator%integrate(U=oscillator, previous=previous, Dt=Dt, t=solution(0, step-steps:step-1))
    endif
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abm_solver

  subroutine am_solver(steps, iterations, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve problem with implicit Adams-Moulton class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P),           intent(IN)  :: steps                 !< Number of steps used: if negative all AM solvers are used.
  integer(I_P), optional, intent(IN)  :: iterations            !< Number of fixed point iterations.
  real(R_P),              intent(OUT) :: solution(0:,0:)       !< Solution at each time step.
  type(tvd_runge_kutta_integrator)    :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter                  :: rk_stages=5           !< Runge-Kutta stages number.
  type(oscillation)                   :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(adams_moulton_integrator)      :: am_integrator         !< Adams-Moulton integrator.
  type(oscillation)                   :: previous(1:steps)     !< Previous time steps solutions.
  type(oscillation)                   :: oscillator            !< Oscillation field.
  integer                             :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call am_integrator%init(steps=steps)
  select case(steps)
  case(0, 1, 2)
    call rk_integrator%init(stages=steps+1)
  case(3)
    call rk_integrator%init(stages=5)
  endselect
  call oscillator%init(initial_state=initial_state, f=f)
  solution = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    if (steps>=step) then
      ! time steps from 1 to s - 1 must be computed with other scheme
      call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = oscillator
    else
      if (steps==0) then
        call am_integrator%integrate(U=oscillator, previous=previous, Dt=Dt, t=solution(0,step-1:step-1), iterations=iterations)
      else
        call am_integrator%integrate(U=oscillator, previous=previous, Dt=Dt, t=solution(0,step-steps:step-1), iterations=iterations)
      endif
    endif
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine am_solver

  subroutine bdf_solver(steps, iterations, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve problem with implicit Bacward-Differentiation-Formula class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P),           intent(IN)  :: steps                 !< Number of steps used: if negative all BDF solvers are used.
  integer(I_P), optional, intent(IN)  :: iterations            !< Number of fixed point iterations.
  real(R_P),              intent(OUT) :: solution(0:,0:)       !< Solution at each time step.
  type(tvd_runge_kutta_integrator)    :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter                  :: rk_stages=5           !< Runge-Kutta stages number.
  type(oscillation)                   :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(bdf_integrator)                :: integrator            !< BDF integrator.
  type(oscillation)                   :: previous(1:steps)     !< Previous time steps solutions.
  type(oscillation)                   :: oscillator            !< Oscillation field.
  integer                             :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call integrator%init(steps=steps)
  select case(steps)
  case(1, 2)
    call rk_integrator%init(stages=steps+1)
  case(3, 4, 5, 6)
    call rk_integrator%init(stages=5)
  endselect
  call oscillator%init(initial_state=initial_state, f=f)
  solution = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    if (steps>=step) then
      ! time steps from 1 to s - 1 must be computed with other scheme
      call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = oscillator
    else
      if (steps==0) then
        call integrator%integrate(U=oscillator, previous=previous, Dt=Dt, t=solution(0,step-1:step-1), iterations=iterations)
      else
        call integrator%integrate(U=oscillator, previous=previous, Dt=Dt, t=solution(0,step-steps:step-1), iterations=iterations)
      endif
    endif
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine bdf_solver

  subroutine emd_rk_solver(stages, tol, solution, last_step)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve problem with explicit embedded Runge-Kutta class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)         :: stages             !< Number of stages used: if negative all RK solvers are used.
  real(R_P),    intent(IN)         :: tol                !< Tolerance on local truncation error.
  real(R_P),    intent(OUT)        :: solution(0:,0:)    !< Solution at each time step.
  integer(I_P), intent(OUT)        :: last_step          !< Last time step number.
  type(emd_runge_kutta_integrator) :: rk_integrator      !< Runge-Kutta integrator.
  type(oscillation)                :: rk_stage(1:stages) !< Runge-Kutta stages.
  type(oscillation)                :: oscillator         !< Oscillation field.
  integer(I_P)                     :: step               !< Time steps counter.
  real(R_P)                        :: Dt_a               !< Adaptive time step.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Dt_a = 10000._R_P ! initial step very large to trigger adaptation
  call rk_integrator%init(stages=stages, tolerance=tol)
  call oscillator%init(initial_state=initial_state, f=f)
  solution = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt_a, t=solution(0, step))
    solution(0, step) = solution(0, step - 1) + Dt_a
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  last_step = step
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine emd_rk_solver

  subroutine euler_solver(solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve problem with explicit forward Euler ODE solver.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(OUT)          :: solution(0:,0:)  !< Solution at each time step.
  type(euler_explicit_integrator) :: euler_integrator !< Euler integrator.
  type(oscillation)               :: oscillator       !< Oscillation field.
  integer(I_P)                    :: step             !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call oscillator%init(initial_state=initial_state, f=f)
  solution = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    call euler_integrator%integrate(U=oscillator, Dt=Dt, t=solution(0, step))
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine euler_solver

  subroutine leapfrog_solver(solution, raw)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve problem with explicit leapfrog class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P),        intent(OUT)    :: solution(0:,0:)       !< Solution at each time step.
  logical, optional, intent(IN)    :: raw                   !< Activate RAW filter.
  type(tvd_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter               :: rk_stages=5           !< Runge-Kutta stages number.
  type(oscillation)                :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(leapfrog_integrator)        :: lf_integrator         !< Leapfrog integrator.
  type(oscillation)                :: previous(1:2)         !< Previous time steps solutions.
  type(oscillation)                :: filter                !< Filter displacement.
  type(oscillation)                :: oscillator            !< Oscillation field.
  integer                          :: step                  !< Time steps counter.
  logical                          :: filtered              !< Activate RAW filter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  filtered = .false. ; if (present(raw)) filtered = raw
  call lf_integrator%init()
  call rk_integrator%init(stages=rk_stages)
  call oscillator%init(initial_state=initial_state, f=f)
  solution = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    if (2>=step) then
      ! time steps from 1 to 2 must be computed with other scheme
      call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
      previous(step) = oscillator
    else
      if (filtered) then
        call lf_integrator%integrate(U=oscillator, previous=previous, Dt=Dt, t=solution(0, step), filter=filter)
      else
        call lf_integrator%integrate(U=oscillator, previous=previous, Dt=Dt, t=solution(0, step))
      endif
    endif
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine leapfrog_solver

  subroutine ls_rk_solver(stages, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve problem with explicit low storage Runge-Kutta class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)        :: stages                !< Number of stages used: if negative all RK solvers are used.
  real(R_P),    intent(OUT)       :: solution(0:,0:)       !< Solution at each time step.
  type(ls_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter              :: registers=2           !< Runge-Kutta stages registers number.
  type(oscillation)               :: rk_stage(1:registers) !< Runge-Kutta registers.
  type(oscillation)               :: oscillator            !< Oscillation field.
  integer(I_P)                    :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call rk_integrator%init(stages=stages)
  call oscillator%init(initial_state=initial_state, f=f)
  solution = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine ls_rk_solver

  subroutine tvd_rk_solver(stages, solution)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve problem with explicit TVD/SSP Runge-Kutta class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)         :: stages             !< Number of stages used: if negative all RK solvers are used.
  real(R_P),    intent(OUT)        :: solution(0:,0:)    !< Solution at each time step.
  type(tvd_runge_kutta_integrator) :: rk_integrator      !< Runge-Kutta integrator.
  type(oscillation)                :: rk_stage(1:stages) !< Runge-Kutta stages.
  type(oscillation)                :: oscillator         !< Oscillation field.
  integer(I_P)                     :: step               !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call rk_integrator%init(stages=stages)
  call oscillator%init(initial_state=initial_state, f=f)
  solution = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine tvd_rk_solver

  subroutine perform_errors_analysis()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Peform errors analysis.
  !<
  !< The analysis is based on the integration of the Oscillation equations with different, decreasing time steps in the range
  !< *[5000, 2500, 1250, 625, 320, 100]*. The error is estimated by the L2 norm of the difference between the exact (\(U_e\)) and
  !< the numerical (\(U_{\Delta t}\)) solutions for each time step:
  !<
  !<$$
  !<\varepsilon (\Delta t) = || U_e - U_{\Delta t} ||_2 = \sqrt{ \sum_{s=1}^{N_s} { \left(U_e(t_0 + s * \Delta t) -
  !<U_{\Delta t}(t_0 + s * \Delta t) \right)^2 }}
  !<$$
  !<
  !< Using two pairs of subsequent-decreasing time steps solution is possible to estimate the order of accuracy of the solver
  !< employed computing the \emph{observed order} of accuracy:
  !<
  !<$$
  !<p = \frac{log10 \left( \frac{\varepsilon (\Delta t_1)}{\varepsilon (\Delta t_2)} \right)}
  !<{log10 \left( \frac{\Delta t_1}{\Delta t_2} \right)}
  !<$$
  !<
  !<where \(\frac{\Delta t_1}{\Delta t_2}>1\).
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), parameter       :: NDt=6                               !< Number of time steps exercised.
  real(R_P),    parameter       :: time_steps(1:NDt)=[5000._R_P, &
                                                      2500._R_P, &
                                                      1250._R_P, &
                                                      625._R_P,  &
                                                      320._R_P,  &
                                                      100._R_P]      !< Time steps exercised.
  real(R_P),    parameter       :: tolerances(1:NDt)=[0.001_R_P,     &
                                                      0.0001_R_P,    &
                                                      0.00001_R_P,   &
                                                      0.000001_R_P,  &
                                                      0.0000001_R_P, &
                                                      0.00000001_R_P]  !< Tolerances for embedded RK solvers.
  integer, parameter            :: ab_steps=4                          !< Adams-Bashforth steps number.
  integer, parameter            :: bdf_steps=6                         !< BDF steps number.
  integer, parameter            :: emd_rk_stages=17                    !< Embedded Runge-Kutta stages number.
  integer, parameter            :: tvd_rk_stages=5                     !< TVD/SSP Runge-Kutta stages number.
  integer, parameter            :: ls_rk_stages=14                     !< Low storage Runge-Kutta stages number.
  integer(I_P)                  :: s                                   !< Steps/stages counter.
  integer(I_P)                  :: d                                   !< Time steps-exercised counter.
  character(len=:), allocatable :: output                              !< Output files basename.
  character(len=:), allocatable :: title                               !< Output files title.
  real(R_P), allocatable        :: solution(:,:)                       !< Solution at each time step.
  integer(I_P)                  :: last_step                           !< Last time step computed.
  real(R_P)                     :: emd_Dt_mean(1:emd_rk_stages, 1:NDt) !< Mean time resolution for adaptive step methods.
  ! errors and orders
  real(R_P) :: ab_errors(1:space_dimension, 1:ab_steps, 1:NDt)            !< Adams-Bashforth errors.
  real(R_P) :: abm_errors(1:space_dimension, 1:ab_steps, 1:NDt)           !< Adams-Bashforth-Moulton errors.
  real(R_P) :: am_errors(1:space_dimension, 1:ab_steps, 1:NDt)            !< Adams-Moulton errors.
  real(R_P) :: bdf_errors(1:space_dimension, 1:bdf_steps, 1:NDt)          !< BDF errors.
  real(R_P) :: emd_rk_errors(1:space_dimension, 1:emd_rk_stages, 1:NDt)   !< TVD/SSP Runge-Kutta errors.
  real(R_P) :: lf_errors(1:space_dimension, 1:2, 1:NDt)                   !< Leapfrog errors.
  real(R_P) :: ls_rk_errors(1:space_dimension, 1:ls_rk_stages, 1:NDt)     !< Low storage Runge-Kutta errors.
  real(R_P) :: tvd_rk_errors(1:space_dimension, 1:tvd_rk_stages, 1:NDt)   !< TVD/SSP Runge-Kutta errors.
  real(R_P) :: ab_orders(1:space_dimension, 1:ab_steps, 1:NDt-1)          !< Adams-Bashforth orders.
  real(R_P) :: abm_orders(1:space_dimension, 1:ab_steps, 1:NDt-1)         !< Adams-Bashforth-Moulton orders.
  real(R_P) :: am_orders(1:space_dimension, 1:ab_steps, 1:NDt-1)          !< Adams-Moulton orders.
  real(R_P) :: bdf_orders(1:space_dimension, 1:bdf_steps, 1:NDt-1)        !< BDF orders.
  real(R_P) :: emd_rk_orders(1:space_dimension, 1:emd_rk_stages, 1:NDt-1) !< Embedded Runge-Kutta orders.
  real(R_P) :: lf_orders(1:space_dimension, 1:2, 1:NDt-1)                 !< Leapfrog orders.
  real(R_P) :: ls_rk_orders(1:space_dimension, 1:ls_rk_stages, 1:NDt-1)   !< Low storage Runge-Kutta orders.
  real(R_P) :: tvd_rk_orders(1:space_dimension, 1:tvd_rk_stages, 1:NDt-1) !< TVD/SSP Runge-Kutta orders.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ab_errors = -1._R_P     ; ab_orders = 0._R_P
  abm_errors = -1._R_P    ; abm_orders = 0._R_P
  am_errors = -1._R_P     ; am_orders = 0._R_P
  bdf_errors = -1._R_P    ; bdf_orders = 0._R_P
  emd_rk_errors = -1._R_P ; emd_rk_orders = 0._R_P
  lf_errors = -1._R_P     ; lf_orders = 0._R_P
  ls_rk_errors = -1._R_P  ; ls_rk_orders = 0._R_P
  tvd_rk_errors = -1._R_P ; tvd_rk_orders = 0._R_P
  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(t_final/Dt)))
  results = .true. ! activate results saving
  ! Analyze errors of Adams-Bashforth solvers
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'adams-bashforth'
    call init(output=output, solution=solution)
    do s=1, ab_steps
      title = 'Oscillation equation integration, explicit Adams-Bashforth, t='//str(n=t_final)//' steps='//trim(str(.true., s))
      call ab_solver(steps=s, solution=solution)
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
      ab_errors(:, s, d) = compute_errors(solution=solution)
    enddo
  enddo
  do s=1, ab_steps
    do d=1, NDt-1
      ab_orders(:, s, d) = estimate_orders(solver_error=ab_errors(:, s, d:d+1), Dt_used=time_steps(d:d+1))
    enddo
  enddo
  ! Analyze errors of Adams-Bashforth-Moulton solvers
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'adams-bashforth-moulton'
    call init(output=output, solution=solution)
    do s=1, ab_steps
      title = 'Oscillation equation integration, predictor-corrector Adams-Bashforth-Moulton, t='//str(n=t_final)//&
              ' steps='//trim(str(.true., s))
      call abm_solver(steps=s, solution=solution)
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
      abm_errors(:, s, d) = compute_errors(solution=solution)
    enddo
  enddo
  do s=1, ab_steps
    do d=1, NDt-1
      abm_orders(:, s, d) = estimate_orders(solver_error=abm_errors(:, s, d:d+1), Dt_used=time_steps(d:d+1))
    enddo
  enddo
  ! Analyze errors of Adams-Moulton solvers
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'adams-moulton'
    call init(output=output, solution=solution)
    do s=0, ab_steps-1
      title = 'Oscillation equation integration, implicit Adams-Moulton, t='//str(n=t_final)//&
              ' steps='//trim(str(.true., s))
      if (implicit_iterations>1) then
        call am_solver(steps=s, solution=solution, iterations=implicit_iterations)
      else
        call am_solver(steps=s, solution=solution)
      endif
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
      am_errors(:, s+1, d) = compute_errors(solution=solution)
    enddo
  enddo
  do s=0, ab_steps-1
    do d=1, NDt-1
      am_orders(:, s+1, d) = estimate_orders(solver_error=am_errors(:, s+1, d:d+1), Dt_used=time_steps(d:d+1))
    enddo
  enddo
  ! Analyze errors of Backward-Differentiation-Formul solvers
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'bdf'
    call init(output=output, solution=solution)
    do s=1, bdf_steps
      title = 'Oscillation equation integration, implicit Backward-Differentiation-Formula, t='//str(n=t_final)//&
              ' steps='//trim(str(.true., s))
      if (implicit_iterations>1) then
        call bdf_solver(steps=s, solution=solution, iterations=implicit_iterations)
      else
        call bdf_solver(steps=s, solution=solution)
      endif
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
      bdf_errors(:, s, d) = compute_errors(solution=solution)
    enddo
  enddo
  do s=1, bdf_steps
    do d=1, NDt-1
      bdf_orders(:, s, d) = estimate_orders(solver_error=bdf_errors(:, s, d:d+1), Dt_used=time_steps(d:d+1))
    enddo
  enddo
  ! Analyze errors of embedded Runge-Kutta solvers
  emd_Dt_mean = 0._R_P
  do d=1, NDt ! loop over exercised tolerances
    Dt = time_steps(NDt)
    tolerance = tolerances(d)
    solver = 'emd-runge-kutta'
    call init(output=output, solution=solution)
    do s=1, emd_rk_stages
      if (s==1) cycle ! 1 stages not allowed
      if (s==2) cycle ! 2 stages not suited for test
      if (s==3) cycle ! 3 stages not yet implemented
      if (s==4) cycle ! 4 stages not yet implemented
      if (s==5) cycle ! 5 stages not yet implemented
      if (s==8) cycle ! 8 stages not yet implemented
      if (s==10) cycle ! 10 stages not yet implemented
      if (s==11) cycle ! 11 stages not yet implemented
      if (s==12) cycle ! 12 stages not yet implemented
      if (s==13) cycle ! 13 stages not yet implemented
      if (s==14) cycle ! 14 stages not yet implemented
      if (s==15) cycle ! 15 stages not yet implemented
      if (s==16) cycle ! 16 stages not yet implemented
      title = 'Oscillation equation integration, explicit embedded Runge-Kutta t='//str(n=t_final)//' stages='//trim(str(.true., s))
      if (s==17) then
        call emd_rk_solver(stages=s, tol=tolerances(d)/1000._R_P, solution=solution, last_step=last_step)
      else
        call emd_rk_solver(stages=s, tol=tolerances(d), solution=solution, last_step=last_step)
      endif
      emd_Dt_mean(s, d) = t_final/real(last_step, kind=R_P)
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution(:, 0:last_step))
      emd_rk_errors(:, s, d) = compute_errors(solution=solution(:, 0:last_step))
    enddo
  enddo
  do s=1, emd_rk_stages
    if (s==1) cycle ! 1 stages not allowed
    if (s==2) cycle ! 2 stages not suited for test
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    if (s==5) cycle ! 5 stages not yet implemented
    if (s==8) cycle ! 8 stages not yet implemented
    if (s==10) cycle ! 10 stages not yet implemented
    if (s==11) cycle ! 11 stages not yet implemented
    if (s==12) cycle ! 12 stages not yet implemented
    if (s==13) cycle ! 13 stages not yet implemented
    if (s==14) cycle ! 14 stages not yet implemented
    if (s==15) cycle ! 15 stages not yet implemented
    if (s==16) cycle ! 16 stages not yet implemented
    do d=1, NDt-1
      emd_rk_orders(:, s, d) = estimate_orders(solver_error=emd_rk_errors(:, s, d:d+1), Dt_used=emd_Dt_mean(s, d:d+1))
    enddo
  enddo
  ! Analyze errors of leapfrog solver
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'leapfrog'
    call init(output=output, solution=solution)
    title = 'Oscillation equation integration, leapfrog (unfiltered), t='//str(n=t_final)
    call leapfrog_solver(solution=solution)
    call save_results(title=title, basename=output, solution=solution)
    lf_errors(:, 1, d) = compute_errors(solution=solution)
    Dt = time_steps(d)
    solver = 'leapfrog'
    call init(output=output, solution=solution)
    title = 'Oscillation equation integration, leapfrog (RAW filtered), t='//str(n=t_final)
    call leapfrog_solver(solution=solution, raw=.true.)
    call save_results(title=title, basename=output//'-raw', solution=solution)
    lf_errors(:, 2, d) = compute_errors(solution=solution)
  enddo
  do d=1, NDt-1
    lf_orders(:, 1, d) = estimate_orders(solver_error=lf_errors(:, 1, d:d+1), Dt_used=time_steps(d:d+1))
    lf_orders(:, 2, d) = estimate_orders(solver_error=lf_errors(:, 2, d:d+1), Dt_used=time_steps(d:d+1))
  enddo
  ! Analyze errors of low storage Runge-Kutta solvers
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'ls-runge-kutta'
    call init(output=output, solution=solution)
    do s=1, ls_rk_stages
      if (s==2) cycle ! 2 stages not yet implemented
      if (s==3) cycle ! 3 stages not yet implemented
      if (s==4) cycle ! 4 stages not yet implemented
      if (s==8) cycle ! 8 stages not yet implemented
      if (s==9) cycle ! 9 stages not yet implemented
      if (s==10) cycle ! 10 stages not yet implemented
      if (s==11) cycle ! 11 stages not yet implemented
      title = 'Oscillation equation integration, explicit low storage Runge-Kutta, t='//str(n=t_final)//&
        ' steps='//trim(str(.true., s))
      call ls_rk_solver(stages=s, solution=solution)
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
      ls_rk_errors(:, s, d) = compute_errors(solution=solution)
    enddo
  enddo
  do s=1, ls_rk_stages
    if (s==2) cycle ! 2 stages not yet implemented
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    if (s==8) cycle ! 8 stages not yet implemented
    if (s==9) cycle ! 9 stages not yet implemented
    if (s==10) cycle ! 10 stages not yet implemented
    if (s==11) cycle ! 11 stages not yet implemented
    do d=1, NDt-1
      ls_rk_orders(:, s, d) = estimate_orders(solver_error=ls_rk_errors(:, s, d:d+1), Dt_used=time_steps(d:d+1))
    enddo
  enddo
  ! Analyze errors of TVD/SSP Runge-Kutta solvers
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'tvd-runge-kutta'
    call init(output=output, solution=solution)
    do s=1, tvd_rk_stages
      if (s==4) cycle ! 4 stages not yet implemented
      title = 'Oscillation equation integration, explicit TVD/SSP Runge-Kutta, t='//str(n=t_final)//' steps='//trim(str(.true., s))
      call tvd_rk_solver(stages=s, solution=solution)
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)), solution=solution)
      tvd_rk_errors(:, s, d) = compute_errors(solution=solution)
    enddo
  enddo
  do s=1, tvd_rk_stages
    if (s==4) cycle ! 4 stages not yet implemented
    do d=1, NDt-1
      tvd_rk_orders(:, s, d) = estimate_orders(solver_error=tvd_rk_errors(:, s, d:d+1), Dt_used=time_steps(d:d+1))
    enddo
  enddo
  print "(A)", "Solver & Time Step & f*Dt & Error X & Error Y & Observed Order of Accuracy X & Observed Order of Accuracy Y"
  call print_analysis(solver='adams-bashforth', time_steps=time_steps, errors=ab_errors, orders=ab_orders)
  call print_analysis(solver='adams-bashforth-moulton', time_steps=time_steps, errors=abm_errors, orders=abm_orders)
  call print_analysis(solver='adams-moulton', time_steps=time_steps, errors=am_errors, orders=am_orders)
  call print_analysis(solver='bdf', time_steps=time_steps, errors=bdf_errors, orders=bdf_orders)
  call print_analysis(solver='emd-runge-kutta', mean_time_steps=emd_Dt_mean, errors=emd_rk_errors, orders=emd_rk_orders)
  call print_analysis(solver='leapfrog', time_steps=time_steps, errors=lf_errors, orders=lf_orders)
  call print_analysis(solver='ls-runge-kutta', time_steps=time_steps, errors=ls_rk_errors, orders=ls_rk_orders)
  call print_analysis(solver='tvd-runge-kutta', time_steps=time_steps, errors=tvd_rk_errors, orders=tvd_rk_orders)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine perform_errors_analysis

  function compute_errors(solution) result(error_L2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the L2 norm of numerical error with respect the exact solution.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN) :: solution(0:,0:)             !< Solution at each time step.
  real(R_P)             :: error_L2(1:space_dimension) !< L2 norm of the numerical error.
  real(R_P)             :: ex_sol(1:space_dimension)   !< Exact solution.
  integer(I_P)          :: s                           !< Steps/stages counter.
  integer(I_P)          :: v                           !< Variables counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  error_L2 = 0._R_P
  do s=0, ubound(solution, dim=2)
    ex_sol = exact_solution(t=solution(0, s))
    do v=1, space_dimension
      error_L2(v) = error_L2(v) + (solution(v, s) - ex_sol(v))**2
    enddo
  enddo
  error_L2 = sqrt(error_L2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_errors

  function estimate_orders(solver_error, Dt_used) result(order)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Estimate the order of accuracy using 2 subsequent refined numerical solutions.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN) :: solver_error(1:space_dimension, 1:2) !< Computed errors.
  real(R_P), intent(IN) :: Dt_used(1:2)                         !< Time steps used.
  real(R_P)             :: order(1:space_dimension)             !< Estimation of the order of accuracy.
  integer(I_P)          :: v                                    !< Variables counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do v=1, space_dimension
    order(v) = log(solver_error(v, 1) / solver_error(v, 2)) / log(Dt_used(1) / Dt_used(2))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction estimate_orders

  subroutine print_analysis(solver, time_steps, mean_time_steps, errors, orders)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Print summary of the error analysis.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*),        intent(IN) :: solver               !< Solver name.
  real(R_P), optional, intent(IN) :: time_steps(:)        !< Time steps exercised.
  real(R_P), optional, intent(IN) :: mean_time_steps(:,:) !< Mean time steps exercised.
  real(R_P),           intent(IN) :: errors(:,:,:)        !< Errors of the solver.
  real(R_P),           intent(IN) :: orders(:,:,:)        !< Observed orders of the solver.
  integer(I_P)                    :: s                    !< Solver index.
  integer(I_P)                    :: d                    !< Time steps-exercised counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(time_steps)) then
    do s=1, size(errors, dim=2)
      if (errors(1,s,1)<0._R_P) cycle
      print "(A)", solver//'-'//trim(str(.true.,s))
      do d=1, size(errors, dim=3)
        if (d==1) then
          print "(A,F8.1,A,F10.3,A,E10.3,A,E10.3,A)", &
            "  & ", time_steps(d), " & ", f*time_steps(d), " & ", &
            errors(1, s, d), " & " , errors(2, s, d), " & / & /"
        else
          print "(A,F8.1,A,F10.3,A,E10.3,A,E10.3,A,F9.2,A,F9.2)", &
            "  & ", time_steps(d), " & ", f*time_steps(d), " & ", &
            errors(1, s, d), " & " , errors(2, s, d), " & ", &
            orders(1, s, d-1), " & " , orders(2, s, d-1)
        endif
      enddo
    enddo
  elseif (present(mean_time_steps)) then
    do s=1, size(errors, dim=2)
      if (errors(1,s,1)<0._R_P) cycle
      print "(A)", solver//'-'//trim(str(.true.,s))
      do d=1, size(errors, dim=3)
        if (d==1) then
          print "(A,F8.1,A,F10.3,A,E10.3,A,E10.3,A)", &
            "  & ", mean_time_steps(s, d), " & ", f*mean_time_steps(s, d), " & ", &
            errors(1, s, d), " & " , errors(2, s, d), " & / & /"
        else
          print "(A,F8.1,A,F10.3,A,E10.3,A,E10.3,A,F9.2,A,F9.2)", &
            "  & ", mean_time_steps(s, d), " & ", f*mean_time_steps(s, d), " & ", &
            errors(1, s, d), " & " , errors(2, s, d), " & ", &
            orders(1, s, d-1), " & " , orders(2, s, d-1)
        endif
      enddo
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_analysis
endprogram integrate_oscillation
