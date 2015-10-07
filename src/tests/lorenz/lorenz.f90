!< Test FOODiE with the integration of Lorenz equations.
program integrate_lorenz
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Lorenz equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, FR_P, str
use type_lorenz, only : lorenz
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
type(Type_Command_Line_Interface) :: cli                                           !< Command line interface handler.
type(lorenz)                      :: attractor                                     !< Lorenz field.
integer,   parameter              :: num_steps=2000                                !< Maximum time steps.
integer,   parameter              :: space_dimension=3                             !< Space dimensions.
real(R_P), parameter              :: sigma=10._R_P                                 !< Lorenz' \(\sigma\).
real(R_P), parameter              :: rho=28._R_P                                   !< Lorenz' \(\rho\).
real(R_P), parameter              :: beta=8._R_P/3._R_P                            !< Lorenz' \(\beta\).
real(R_P), parameter              :: dt=0.01_R_P                                   !< Time step.
real(R_P), parameter              :: initial_state(1:space_dimension)=[1., 1., 1.] !< Initial state.
real(R_P)                         :: solution(0:space_dimension, 0:num_steps)      !< Solution at each time step.
integer(I_P)                      :: error                                         !< Error handler.
character(99)                     :: solver                                        !< Solver used.
logical                           :: plots                                         !< Flag for activating plots saving.
logical                           :: results                                       !< Flag for activating results saving.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'lorenz',                                              &
              authors     = 'Fortran-FOSS-Programmers',                            &
              license     = 'GNU GPLv3',                                           &
              description = 'Test FOODiE library on Lorenz equations integration', &
              examples    = ["lorenz --solver euler --results  ",                  &
                             "lorenz --solver ls-runge-kutta -r",                  &
                             "lorenz --solver adams-bashforth  ",                  &
                             "lorenz --solver all --plots -r   "])
call cli%add(switch='--solver', switch_ab='-s', help='ODE solver used', required=.true., act='store')
call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.')
call cli%add(switch='--plots', switch_ab='-p', help='Save plots of results', required=.false., act='store_true', def='.false.')
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='-s', val=solver, error=error) ; if (error/=0) stop
call cli%get(switch='-r', val=results, error=error) ; if (error/=0) stop
call cli%get(switch='-p', val=plots, error=error) ; if (error/=0) stop
! integrate Lorenz equations
select case(trim(adjustl(solver)))
case('adams-bashforth')
  call test_ab()
case('euler')
  call test_euler()
case('leapfrog')
  call test_leapfrog()
case('ls-runge-kutta')
  call test_ls_rk()
case('tvd-runge-kutta')
  call test_tvd_rk()
case('all')
  call test_ab()
  call test_euler()
  call test_leapfrog()
  call test_ls_rk()
  call test_tvd_rk()
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
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine save_results(title, filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save plots of results.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN) :: title    !< Plot title.
  character(*), intent(IN) :: filename !< Output filename.
  integer(I_P)             :: rawfile  !< Raw file unit for saving results.
  type(pyplot)             :: plt      !< Plot file handler.
  integer(I_P)             :: i        !< Counter.
  integer(I_P)             :: s        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (results) then
    open(newunit=rawfile, file=filename//'.dat')
    write(rawfile, '(A)')'# '//title
    write(rawfile, '(A)')'# VARIABLES: "t" "x" "y" "z"'
    do s=1, num_steps
      write(rawfile, '(4('//FR_P//',1X))')(solution(i, s), i=0, 3)
    enddo
    close(rawfile)
  endif
  if (plots) then
    call plt%initialize(grid=.true., xlabel='time', title=title, legend=.true.)
    call plt%add_plot(x=solution(0, :), y=solution(1, :), label='x', linestyle='r-', linewidth=1)
    call plt%add_plot(x=solution(0, :), y=solution(2, :), label='y', linestyle='b-', linewidth=1)
    call plt%add_plot(x=solution(0, :), y=solution(3, :), label='z', linestyle='g-', linewidth=1)
    call plt%savefig(filename//'.png')

    call plt%initialize(grid=.true., title=title//'-path', legend=.true.)
    call plt%add_plot(x=solution(1, :), y=solution(2, :), label='x-y path', linestyle='r-', linewidth=1)
    call plt%add_plot(x=solution(1, :), y=solution(3, :), label='x-z path', linestyle='g-', linewidth=1)
    call plt%add_plot(x=solution(2, :), y=solution(3, :), label='y-z path', linestyle='b-', linewidth=1)
    call plt%savefig('path-'//filename//'.png')
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
  type(lorenz)                     :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(adams_bashforth_integrator) :: ab_integrator         !< Adams-Bashforth integrator.
  integer, parameter               :: ab_steps=4            !< Adams-Bashforth steps number.
  type(lorenz)                     :: previous(1:ab_steps)  !< Previous time steps solutions.
  integer(I_P)                     :: s                     !< AB steps counter.
  integer                          :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Lorenz equations by means of Adams-Bashforth class of solvers'
  do s=1, ab_steps
    print "(A)", ' AB-'//trim(str(.true.,s))
    call ab_integrator%init(steps=s)
    select case(s)
    case(1, 2, 3)
      call rk_integrator%init(stages=s)
    case(4)
      call rk_integrator%init(stages=5)
    endselect
    call attractor%init(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta, steps=s)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = attractor%output()
    do step = 1, num_steps
      if (s>=step) then
        ! the time steps from 1 to s - 1 must be computed with other scheme...
        call rk_integrator%integrate(U=attractor, stage=rk_stage, dt=dt, t=solution(0, step))
        previous(step) = attractor
      else
        call ab_integrator%integrate(U=attractor, previous=previous(1:s), dt=dt, t=solution(0, step-s:step-1))
      endif
      solution(0, step) = step * dt
      solution(1:space_dimension, step) = attractor%output()
    enddo
    call save_results(title='FOODiE test: Lorenz equation integration, explicit Adams-Bashforth '//trim(str(.true., s))//' steps', &
                      filename='lorenz_integration-ab-'//trim(str(.true., s)))
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
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Lorenz equations by means of explicit Euler solver'
  call attractor%init(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
  solution(0, 0) = 0._R_P
  solution(1:space_dimension, 0) = attractor%output()
  do step = 1, num_steps
    call euler_integrator%integrate(U=attractor, dt=dt, t=solution(0, step))
    solution(0, step) = step * dt
    solution(1:space_dimension, step) = attractor%output()
  enddo
  call save_results(title='FOODiE test: Lorenz equation integration, explicit Euler', &
                    filename='lorenz_integration-euler')
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
  type(lorenz)                     :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(lorenz)                     :: filter                !< Filter displacement.
  type(leapfrog_integrator)        :: lf_integrator         !< Leapfrog integrator.
  type(lorenz)                     :: previous(1:2)         !< Previous time steps solutions.
  integer                          :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Lorenz equations by means of leapfrog (RAW filtered) class of solvers'
  call lf_integrator%init(nu=1.0_R_P, alpha=0._R_P)
  call rk_integrator%init(stages=rk_stages)
  call attractor%init(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta, steps=2)
  solution(0, 0) = 0._R_P
  solution(1:space_dimension, 0) = attractor%output()
  do step = 1, num_steps
    if (2>=step) then
      ! the time steps from 1 to 2 must be computed with other scheme...
      call rk_integrator%integrate(U=attractor, stage=rk_stage, dt=dt, t=solution(0, step))
      previous(step) = attractor
    else
      call lf_integrator%integrate(U=attractor, previous=previous, dt=dt, t=solution(0, step), filter=filter)
    endif
    solution(0, step) = step * dt
    solution(1:space_dimension, step) = attractor%output()
  enddo
  call save_results(title='FOODiE test: Lorenz equation integration, explicit leapfrog scheme',&
                    filename='lorenz_integration-lf-RAW-filter')
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
  type(lorenz)                    :: rk_stage(1:registers) !< Runge-Kutta stages.
  integer(I_P)                    :: s                     !< RK stages counter.
  integer(I_P)                    :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Lorenz equations by means of low storage (2N) Runge-Kutta class of solvers'
  do s=1, rk_stages
    if (s==2) cycle ! 2 stages not yet implemented
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    call rk_integrator%init(stages=s)
    call attractor%init(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = attractor%output()
    do step = 1, num_steps
      call rk_integrator%integrate(U=attractor, stage=rk_stage, dt=dt, t=solution(0, step))
      solution(0, step) = step * dt
      solution(1:space_dimension, step) = attractor%output()
    enddo
    call save_results(title='FOODiE test: Lorenz equation integration, explicit low storage Runge-Kutta '//trim(str(.true., s))//&
                      ' stages', filename='lorenz_integration-lsrk-'//trim(str(.true., s)))
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
  type(lorenz)                     :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  integer(I_P)                     :: s                     !< RK stages counter.
  integer(I_P)                     :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Lorenz equations by means of TVD/SSP Runge-Kutta class of solvers'
  do s=1, rk_stages
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    call rk_integrator%init(stages=s)
    call attractor%init(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = attractor%output()
    do step = 1, num_steps
      call rk_integrator%integrate(U=attractor, stage=rk_stage(1:s), dt=dt, t=solution(0, step))
      solution(0, step) = step * dt
      solution(1:space_dimension, step) = attractor%output()
    enddo
    call save_results(title='FOODiE test: Lorenz equation integration, explicit TVD Runge-Kutta '//trim(str(.true., s))//' stages',&
                      filename='lorenz_integration-tvdrk-'//trim(str(.true., s)))
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_tvd_rk
endprogram integrate_lorenz
