!< Test FOODiE with the integration of Oscillation equations.
program integrate_oscillation
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Oscillation equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, FR_P, str, strz
use type_oscillation, only : oscillation
use Data_Type_Command_Line_Interface, only : Type_Command_Line_Interface
use foodie, only : adams_bashforth_integrator, &
                   euler_explicit_integrator, &
                   leapfrog_integrator, &
                   ls_runge_kutta_integrator, &
                   tvd_runge_kutta_integrator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Command_Line_Interface) :: cli                                               !< Command line interface handler.
type(oscillation)                 :: oscillator                                        !< Oscillation field.
integer,   parameter              :: space_dimension=2                                 !< Space dimensions.
real(R_P), parameter              :: initial_state(1:space_dimension)=[0._R_P, 1._R_P] !< Initial state.
real(R_P)                         :: f                                                 !< Oscillation frequency.
real(R_P)                         :: t_final                                           !< Final time.
real(R_P)                         :: Dt                                                !< Time step.
real(R_P), allocatable            :: solution(:,:)                                     !< Solution at each time step.
integer(I_P)                      :: error                                             !< Error handler.
integer(I_P)                      :: stages_steps                                      !< Number of stages/steps used.
character(99)                     :: solver                                            !< Solver used.
logical                           :: plots                                             !< Flag for activating plots saving.
logical                           :: results                                           !< Flag for activating results saving.
character(len=:), allocatable     :: output                                            !< Output files basename.
character(99)                     :: output_cli                                        !< Output files basename.
logical                           :: errors_analysis                                   !< Flag for activating errors analysis.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'oscillation',                                              &
              authors     = 'Fortran-FOSS-Programmers',                                 &
              license     = 'GNU GPLv3',                                                &
              description = 'Test FOODiE library on Oscillation equations integration', &
              examples    = ["oscillation --solver euler --results  ",                  &
                             "oscillation --solver ls-runge-kutta -r",                  &
                             "oscillation --solver adams-bashforth  ",                  &
                             "oscillation --solver all --plots -r   "])
call cli%add(switch='--solver', switch_ab='-s', help='ODE solver used', required=.false., def='all', act='store')
call cli%add(switch='--ss', help='Stages/steps used', required=.false., act='store', def='-1')
call cli%add(switch='--frequency', switch_ab='-f', help='Oscillation frequency', required=.false., def='1e-4', act='store')
call cli%add(switch='--time_step', switch_ab='-Dt', help='Integration time step', required=.false., def='100.d0', act='store')
call cli%add(switch='--t_final', switch_ab='-tf', help='Final integration time', required=.false., def='1e6', act='store')
call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.')
call cli%add(switch='--plots', switch_ab='-p', help='Save plots of results', required=.false., act='store_true', def='.false.')
call cli%add(switch='--output', help='Output files basename', required=.false., act='store', def='unset')
call cli%add(switch='--errors_analysis', help='Peform errors analysis', required=.false., act='store_true', def='.false.')
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='-s', val=solver, error=error) ; if (error/=0) stop
call cli%get(switch='--ss', val=stages_steps, error=error) ; if (error/=0) stop
call cli%get(switch='-f', val=f, error=error) ; if (error/=0) stop
call cli%get(switch='-Dt', val=Dt, error=error) ; if (error/=0) stop
call cli%get(switch='-tf', val=t_final, error=error) ; if (error/=0) stop
call cli%get(switch='-r', val=results, error=error) ; if (error/=0) stop
call cli%get(switch='-p', val=plots, error=error) ; if (error/=0) stop
call cli%get(switch='--output', val=output_cli, error=error) ; if (error/=0) stop
call cli%get(switch='--errors_analysis', val=errors_analysis, error=error) ; if (error/=0) stop
if (errors_analysis) then
  call perform_errors_analysis
else
  ! Integrate Oscillation equations
  call init
  select case(trim(adjustl(solver)))
  case('adams-bashforth')
    call test_ab
  case('euler')
    call test_euler
  case('leapfrog')
    call test_leapfrog
  case('ls-runge-kutta')
    call test_ls_rk
  case('tvd-runge-kutta')
    call test_tvd_rk
  case('all')
    call test_ab
    call test_euler
    call test_leapfrog
    call test_ls_rk
    call test_tvd_rk
  case default
    print "(A)", 'Error: unknown solver "'//trim(adjustl(solver))//'"'
    print "(A)", 'Valid solver names are:'
    print "(A)", '  + adams-bashforth'
    print "(A)", '  + euler'
    print "(A)", '  + leapfrog'
    print "(A)", '  + ls-runge-kutta'
    print "(A)", '  + tvd-runge-kutta'
    print "(A)", '  + all'
  endselect
endif
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine init()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Initialize the integration variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (mod(t_final, Dt)/=0) then
    print "(A)", 'Error: the final integration time must be an exact multiple of the time step used!'
    print "(A)", 'Final integration time: '//str(.true.,t_final)
    print "(A)", 'Time step: '//str(.true.,Dt)
    stop
  endif
  if (allocated(solution)) deallocate(solution) ; allocate(solution(0:space_dimension, 0:int(t_final/Dt)))
  if (trim(adjustl(output_cli))/='unset') then
    output = trim(adjustl(output_cli))//'-'//trim(strz(10,int(t_final/Dt)))//'-time_steps'//'-'//trim(adjustl(solver))
  else
    output = 'oscillation_integration-'//trim(strz(10,int(t_final/Dt)))//'-time_steps'//'-'//trim(adjustl(solver))
  endif
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

  subroutine save_results(title, basename)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save plots of results.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN) :: title                     !< Plot title.
  character(*), intent(IN) :: basename                  !< Output files basename.
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
    write(rawfile, '(A)')'ZONE T="FOODiE time serie"'
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

  subroutine test_ab()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit Adams-Bashforth class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(tvd_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter               :: rk_stages=5           !< Runge-Kutta stages number.
  type(oscillation)                :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(adams_bashforth_integrator) :: ab_integrator         !< Adams-Bashforth integrator.
  integer, parameter               :: ab_steps=4            !< Adams-Bashforth steps number.
  integer                          :: step                  !< Time steps counter.
  integer(I_P)                     :: s                     !< AB steps counter.
  integer(I_P)                     :: steps_range(1:2)      !< Steps used.
  character(len=:), allocatable    :: title                 !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Oscillation equations by means of Adams-Bashforth class of solvers'
  steps_range = [1, ab_steps] ; if (stages_steps>0) steps_range = [stages_steps, stages_steps]
  do s=steps_range(1), steps_range(2)
    print "(A)", ' AB-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, explicit Adams-Bashforth, t='//str(n=t_final)//' steps='//trim(str(.true., s))
    call ab_integrator%init(steps=s)
    select case(s)
    case(1, 2, 3)
      call rk_integrator%init(stages=s)
    case(4)
      call rk_integrator%init(stages=5)
    endselect
    call oscillator%init(initial_state=initial_state, f=f, steps=s)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = oscillator%output()
    step = 0
    do while(solution(0, step)<t_final)
      step = step + 1
      if (s>=step) then
        ! the time steps from 1 to s - 1 must be computed with other scheme...
        call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
      else
        call ab_integrator%integrate(U=oscillator, Dt=Dt, t=solution(0, step-s:step-1))
      endif
      solution(0, step) = step * Dt
      solution(1:space_dimension, step) = oscillator%output()
    enddo
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)))
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_ab

  subroutine test_euler()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit forward Euler ODE solver.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(euler_explicit_integrator) :: euler_integrator !< Euler integrator.
  integer(I_P)                    :: step             !< Time steps counter.
  character(len=:), allocatable   :: title            !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Oscillation equations by means of explicit Euler solver'
  title = 'Oscillation equation integration, explicit Euler, t='//str(n=t_final)
  call oscillator%init(initial_state=initial_state, f=f)
  solution(0, 0) = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    call euler_integrator%integrate(U=oscillator, Dt=Dt, t=solution(0, step))
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  call save_results(title=title, basename=output)
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_euler

  subroutine test_leapfrog()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit leapfrog class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(tvd_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter               :: rk_stages=5           !< Runge-Kutta stages number.
  type(oscillation)                :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(oscillation)                :: filter                !< Filter displacement.
  type(leapfrog_integrator)        :: lf_integrator         !< Leapfrog integrator.
  integer                          :: step                  !< Time steps counter.
  character(len=:), allocatable    :: title                 !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Oscillation equations by means of leapfrog (RAW filtered) class of solvers'
  title = 'Oscillation equation integration, leapfrog (RAW filtered), t='//str(n=t_final)
  call lf_integrator%init()
  call rk_integrator%init(stages=rk_stages)
  call oscillator%init(initial_state=initial_state, f=f, steps=2)
  solution(0, 0) = 0._R_P
  solution(1:space_dimension, 0) = oscillator%output()
  step = 0
  do while(solution(0, step)<t_final)
    step = step + 1
    if (2>=step) then
      ! the time steps from 1 to 2 must be computed with other scheme...
      call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
    else
      call lf_integrator%integrate(U=oscillator, filter=filter, Dt=Dt, t=solution(0, step))
    endif
    solution(0, step) = step * Dt
    solution(1:space_dimension, step) = oscillator%output()
  enddo
  call save_results(title=title, basename=output)
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_leapfrog

  subroutine test_ls_rk()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit low storage Runge-Kutta class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(ls_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter              :: rk_stages=5           !< Runge-Kutta stages number.
  integer, parameter              :: registers=2           !< Runge-Kutta stages number.
  type(oscillation)               :: rk_stage(1:registers) !< Runge-Kutta registers.
  integer(I_P)                    :: s                     !< RK stages counter.
  integer(I_P)                    :: step                  !< Time steps counter.
  integer(I_P)                    :: stages_range(1:2)     !< Stages used.
  character(len=:), allocatable   :: title                 !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Oscillation equations by means of low storage (2N) Runge-Kutta class of solvers'
  stages_range = [1, rk_stages] ; if (stages_steps>0) stages_range = [stages_steps, stages_steps]
  do s=stages_range(1), stages_range(2)
    if (s==2) cycle ! 2 stages not yet implemented
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, explicit low storage Runge-Kutta, t='//str(n=t_final)//&
      ' steps='//trim(str(.true., s))
    call rk_integrator%init(stages=s)
    call oscillator%init(initial_state=initial_state, f=f)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = oscillator%output()
    step = 0
    do while(solution(0, step)<t_final)
      step = step + 1
      call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
      solution(0, step) = step * Dt
      solution(1:space_dimension, step) = oscillator%output()
    enddo
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)))
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_ls_rk

  subroutine test_tvd_rk()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit TVD/SSP Runge-Kutta class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(tvd_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter               :: rk_stages=5           !< Runge-Kutta stages number.
  type(oscillation)                :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  integer(I_P)                     :: s                     !< RK stages counter.
  integer(I_P)                     :: step                  !< Time steps counter.
  integer(I_P)                     :: stages_range(1:2)     !< Stages used.
  character(len=:), allocatable    :: title                 !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Oscillation equations by means of TVD/SSP Runge-Kutta class of solvers'
  stages_range = [1, rk_stages] ; if (stages_steps>0) stages_range = [stages_steps, stages_steps]
  do s=stages_range(1), stages_range(2)
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    title = 'Oscillation equation integration, explicit TVD/SSP Runge-Kutta, t='//str(n=t_final)//' steps='//trim(str(.true., s))
    call rk_integrator%init(stages=s)
    call oscillator%init(initial_state=initial_state, f=f)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = oscillator%output()
    step = 0
    do while(solution(0, step)<t_final)
      step = step + 1
      call rk_integrator%integrate(U=oscillator, stage=rk_stage(1:s), Dt=Dt, t=solution(0, step))
      solution(0, step) = step * Dt
      solution(1:space_dimension, step) = oscillator%output()
    enddo
    call save_results(title=title, basename=output//'-'//trim(str(.true., s)))
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_tvd_rk

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
  integer(I_P), parameter          :: NDt=6                        !< Number of time steps exercised.
  real(R_P),    parameter          :: time_steps(1:NDt)=[5000._R_P, &
                                                         2500._R_P, &
                                                         1250._R_P, &
                                                         625._R_P,  &
                                                         320._R_P,  &
                                                         100._R_P] !< Time steps exercised.
  type(adams_bashforth_integrator) :: ab_integrator                !< Adams-Bashforth integrator.
  type(leapfrog_integrator)        :: lf_integrator                !< Leapfrog integrator.
  type(ls_runge_kutta_integrator)  :: ls_integrator                !< Low storage Runge-Kutta integrator.
  type(tvd_runge_kutta_integrator) :: rk_integrator                !< TVD/SSP Runge-Kutta integrator.
  integer, parameter               :: ab_steps=4                   !< Adams-Bashforth steps number.
  integer, parameter               :: rk_stages=5                  !< Runge-Kutta stages number.
  type(oscillation)                :: rk_stage(1:rk_stages)        !< Runge-Kutta stages.
  type(oscillation)                :: filter                       !< Filter displacement.
  real(R_P)                        :: ab_errors(1:space_dimension,&
                                                1:ab_steps,       &
                                                1:NDt)             !< Adams-Bashforth errors.
  real(R_P)                        :: lf_errors(1:space_dimension,&
                                                1:NDt)             !< Leapfrog errors.
  real(R_P)                        :: ls_errors(1:space_dimension,&
                                                1:rk_stages,      &
                                                1:NDt)             !< Low storage Runge-Kutta errors.
  real(R_P)                        :: rk_errors(1:space_dimension,&
                                                1:rk_stages,      &
                                                1:NDt)             !< TVD/SSP Runge-Kutta errors.
  real(R_P)                        :: ab_orders(1:space_dimension,&
                                                1:ab_steps,       &
                                                1:NDt-1)           !< Adams-Bashforth orders.
  real(R_P)                        :: lf_orders(1:space_dimension,&
                                                1:NDt-1)           !< Leapfrog orders.
  real(R_P)                        :: ls_orders(1:space_dimension,&
                                                1:rk_stages,      &
                                                1:NDt-1)           !< Low storage Runge-Kutta orders.
  real(R_P)                        :: rk_orders(1:space_dimension,&
                                                1:rk_stages,      &
                                                1:NDt-1)           !< TVD/SSP Runge-Kutta orders.
  integer                          :: step                         !< Time steps counter.
  integer(I_P)                     :: s                            !< Steps/stages counter.
  integer(I_P)                     :: d                            !< Time steps-exercised counter.
  character(len=:), allocatable    :: title                        !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  results = .true. ! activate results saving
  ! Analyze errors of Adams-Bashforth solvers
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'adams-bashforth'
    call init
    do s=1, ab_steps
      if (s==1) cycle ! 1st order scheme surely not stable for this test
      title = 'Oscillation equation integration, explicit Adams-Bashforth, t='//str(n=t_final)//' steps='//trim(str(.true., s))
      call ab_integrator%init(steps=s)
      select case(s)
      case(1, 2, 3)
        call rk_integrator%init(stages=s)
      case(4)
        call rk_integrator%init(stages=5)
      endselect
      call oscillator%init(initial_state=initial_state, f=f, steps=s)
      solution(0, 0) = 0._R_P
      solution(1:space_dimension, 0) = oscillator%output()
      step = 0
      do while(solution(0, step)<t_final)
        step = step + 1
        if (s>=step) then
          ! the time steps from 1 to s - 1 must be computed with other scheme...
          call rk_integrator%integrate(U=oscillator, stage=rk_stage(1:s), Dt=Dt, t=solution(0, step))
        else
          call ab_integrator%integrate(U=oscillator, Dt=Dt, t=solution(0, step-s:step-1))
        endif
        solution(0, step) = step * Dt
        solution(1:space_dimension, step) = oscillator%output()
      enddo
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)))
      ab_errors(:, s, d) = compute_errors()
    enddo
  enddo
  do s=1, ab_steps
    if (s==1) cycle ! 1st order scheme surely not stable for this test
    do d=1, NDt-1
      ab_orders(:, s, d) = estimate_orders(solver_error=ab_errors(:, s, d:d+1), Dt_used=time_steps(d:d+1))
    enddo
  enddo
  ! Analyze errors of leapfrog solver
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'leapfrog'
    call init
    title = 'Oscillation equation integration, leapfrog (RAW filtered), t='//str(n=t_final)
    call lf_integrator%init()
    call rk_integrator%init(stages=rk_stages)
    call oscillator%init(initial_state=initial_state, f=f, steps=2)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = oscillator%output()
    step = 0
    do while(solution(0, step)<t_final)
      step = step + 1
      if (2>=step) then
        ! the time steps from 1 to 2 must be computed with other scheme...
        call rk_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
      else
        call lf_integrator%integrate(U=oscillator, filter=filter, Dt=Dt, t=solution(0, step))
      endif
      solution(0, step) = step * Dt
      solution(1:space_dimension, step) = oscillator%output()
    enddo
    call save_results(title=title, basename=output)
    lf_errors(:, d) = compute_errors()
  enddo
  do d=1, NDt-1
    lf_orders(:, d) = estimate_orders(solver_error=lf_errors(:, d:d+1), Dt_used=time_steps(d:d+1))
  enddo
  ! Analyze errors of low storage Runge-Kutta solvers
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'ls-runge-kutta'
    call init
    do s=1, rk_stages
      if (s==1) cycle ! 1st order scheme surely not stable for this test
      if (s==2) cycle ! 2 stages not yet implemented
      if (s==3) cycle ! 3 stages not yet implemented
      if (s==4) cycle ! 4 stages not yet implemented
      title = 'Oscillation equation integration, explicit low storage Runge-Kutta, t='//str(n=t_final)//&
        ' steps='//trim(str(.true., s))
      call ls_integrator%init(stages=s)
      call oscillator%init(initial_state=initial_state, f=f)
      solution(0, 0) = 0._R_P
      solution(1:space_dimension, 0) = oscillator%output()
      step = 0
      do while(solution(0, step)<t_final)
        step = step + 1
        call ls_integrator%integrate(U=oscillator, stage=rk_stage, Dt=Dt, t=solution(0, step))
        solution(0, step) = step * Dt
        solution(1:space_dimension, step) = oscillator%output()
      enddo
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)))
      ls_errors(:, s, d) = compute_errors()
    enddo
  enddo
  do s=1, rk_stages
    if (s==1) cycle ! 1st order scheme surely not stable for this test
    if (s==2) cycle ! 2 stages not yet implemented
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    do d=1, NDt-1
      ls_orders(:, s, d) = estimate_orders(solver_error=ls_errors(:, s, d:d+1), Dt_used=time_steps(d:d+1))
    enddo
  enddo
  ! Analyze errors of TVD/SSP Runge-Kutta solvers
  do d=1, NDt ! loop over exercised time steps
    Dt = time_steps(d)
    solver = 'tvd-runge-kutta'
    call init
    do s=1, rk_stages
      if (s==1) cycle ! 1st order scheme surely not stable for this test
      if (s==4) cycle ! 4 stages not yet implemented
      title = 'Oscillation equation integration, explicit TVD/SSP Runge-Kutta, t='//str(n=t_final)//' steps='//trim(str(.true., s))
      call rk_integrator%init(stages=s)
      call oscillator%init(initial_state=initial_state, f=f)
      solution(0, 0) = 0._R_P
      solution(1:space_dimension, 0) = oscillator%output()
      step = 0
      do while(solution(0, step)<t_final)
        step = step + 1
        call rk_integrator%integrate(U=oscillator, stage=rk_stage(1:s), Dt=Dt, t=solution(0, step))
        solution(0, step) = step * Dt
        solution(1:space_dimension, step) = oscillator%output()
      enddo
      call save_results(title=title, basename=output//'-'//trim(str(.true., s)))
      rk_errors(:, s, d) = compute_errors()
    enddo
  enddo
  do s=1, rk_stages
    if (s==1) cycle ! 1st order scheme surely not stable for this test
    if (s==4) cycle ! 4 stages not yet implemented
    do d=1, NDt-1
      rk_orders(:, s, d) = estimate_orders(solver_error=rk_errors(:, s, d:d+1), Dt_used=time_steps(d:d+1))
    enddo
  enddo
  call print_analysis(ab_steps=ab_steps, rk_stages=rk_stages, NDt=NDt, time_steps=time_steps, &
                      ab_errors=ab_errors, lf_errors=lf_errors, ls_errors=ls_errors, rk_errors=rk_errors, &
                      ab_orders=ab_orders, lf_orders=lf_orders, ls_orders=ls_orders, rk_orders=rk_orders)

  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine perform_errors_analysis

  function compute_errors() result(error_L2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the L2 norm of numerical error with respect the exact solution.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P)    :: error_L2(1:space_dimension) !< L2 norm of the numerical error.
  real(R_P)    :: ex_sol(1:space_dimension)   !< Exact solution.
  integer(I_P) :: s                           !< Steps/stages counter.
  integer(I_P) :: v                           !< Variables counter.
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

  subroutine print_analysis(ab_steps, rk_stages, NDt, time_steps, &
                            ab_errors, lf_errors, ls_errors, rk_errors, ab_orders, lf_orders, ls_orders, rk_orders)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Print summary of the error analysis.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN) :: ab_steps         !< Adams-Bashforth steps number.
  integer(I_P), intent(IN) :: rk_stages        !< Runge-Kutta stages number.
  integer(I_P), intent(IN) :: NDt              !< Number of time steps exercised.
  real(R_P),    intent(IN) :: time_steps(:)    !< Time steps exercised.
  real(R_P),    intent(IN) :: ab_errors(:,:,:) !< Adams-Bashforth errors.
  real(R_P),    intent(IN) :: lf_errors(:,:)   !< Leapfrog errors.
  real(R_P),    intent(IN) :: ls_errors(:,:,:) !< Low storage Runge-Kutta errors.
  real(R_P),    intent(IN) :: rk_errors(:,:,:) !< TVD/SSP Runge-Kutta errors.
  real(R_P),    intent(IN) :: ab_orders(:,:,:) !< Adams-Bashforth orders.
  real(R_P),    intent(IN) :: lf_orders(:,:)   !< Leapfrog orders.
  real(R_P),    intent(IN) :: ls_orders(:,:,:) !< Low storage Runge-Kutta orders.
  real(R_P),    intent(IN) :: rk_orders(:,:,:) !< TVD/SSP Runge-Kutta orders.
  integer(I_P)             :: s                !< Steps/stages counter.
  integer(I_P)             :: d                !< Time steps-exercised counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", "Solver & Time Step & Error X & Error Y & Observed Order of Accuracy X & Observed Order of Accuracy Y"
  do s=1, ab_steps
    if (s==1) cycle ! 1st order scheme surely not stable for this test
    print "(A)", "Adams-Bashforth, "//trim(str(.true.,s))//" steps"
    do d=1, NDt
      if (d==1) then
        print "(A,F8.1,A,E10.3,A,E10.3,A)", "  & ", time_steps(d), " & ", &
                                         ab_errors(1, s, d), " & " , ab_errors(2, s, d), " & / & /"
      else
        print "(A,F8.1,A,E10.3,A,E10.3,A,F9.2,A,F9.2)", "  & ", time_steps(d), " & ", &
                                                        ab_errors(1, s, d), " & " , ab_errors(2, s, d), " & ", &
                                                        ab_orders(1, s, d-1), " & " , ab_orders(2, s, d-1)
      endif
    enddo
  enddo
  print "(A)", "Leapfrog 2 steps"
  do d=1, NDt
    if (d==1) then
      print "(A,F8.1,A,E10.3,A,E10.3,A)", "  & ", time_steps(d), " & ", &
                                          lf_errors(1, d), " & " , lf_errors(2, d), " & / & /"
    else
      print "(A,F8.1,A,E10.3,A,E10.3,A,F9.2,A,F9.2)", "  & ", time_steps(d), " & ", &
                                                      lf_errors(1, d), " & " , lf_errors(2, d), " & ", &
                                                      lf_orders(1, d-1), " & " , lf_orders(2, d-1)
    endif
  enddo
  do s=1, rk_stages
    if (s==1) cycle ! 1st order scheme surely not stable for this test
    if (s==2) cycle ! 2 stages not yet implemented
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", "Low Storage Runge-Kutta, "//trim(str(.true.,s))//" stages"
    do d=1, NDt
      if (d==1) then
        print "(A,F8.1,A,E10.3,A,E10.3,A)", "  & ", time_steps(d), " & ", &
                                            ls_errors(1, s, d), " & " , ls_errors(2, s, d), " & / & /"
      else
        print "(A,F8.1,A,E10.3,A,E10.3,A,F9.2,A,F9.2)", "  & ", time_steps(d), " & ", &
                                                        ls_errors(1, s, d), " & " , ls_errors(2, s, d), " & ", &
                                                        ls_orders(1, s, d-1), " & " , ls_orders(2, s, d-1)
      endif
    enddo
  enddo
  do s=1, rk_stages
    if (s==1) cycle ! 1st order scheme surely not stable for this test
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", "TVD/SSP Runge-Kutta, "//trim(str(.true.,s))//" stages"
    do d=1, NDt
      if (d==1) then
        print "(A,F8.1,A,E10.3,A,E10.3,A)", "  & ", time_steps(d), " & ", &
                                             rk_errors(1, s, d), " & " , rk_errors(2, s, d), " & / & /"
      else
        print "(A,F8.1,A,E10.3,A,E10.3,A,F9.2,A,F9.2)", "  & ", time_steps(d), " & ", &
                                                        rk_errors(1, s, d), " & " , rk_errors(2, s, d), " & ", &
                                                        rk_orders(1, s, d-1), " & " , rk_orders(2, s, d-1)
      endif
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_analysis
endprogram integrate_oscillation
