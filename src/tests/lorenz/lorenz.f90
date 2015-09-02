!< Test FOODiE with the integration of Lorenz equations.
program integrate_lorenz
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Lorenz equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, str
use type_lorenz, only : lorenz
use Data_Type_Command_Line_Interface, only : Type_Command_Line_Interface
use foodie, only : euler_explicit_integrator, tvd_runge_kutta_integrator, adams_bashforth_integrator
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
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'lorenz',                                              &
              authors     = 'Fortran-FOSS-Programmers',                            &
              license     = 'GNU GPLv3',                                           &
              description = 'Test FOODiE library on Lorenz equations integration', &
              examples    = ["lorenz --solver euler          ",                    &
                             "lorenz --solver runge-kutta    ",                    &
                             "lorenz --solver adams-bashforth",                    &
                             "lorenz --solver all            "])
call cli%add(switch='--solver', switch_ab='-s', help='ODE solver used', required=.true., act='store', error=error)
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='-s', val=solver, error=error) ; if (error/=0) stop
! integrate Lorenz equations
select case(trim(adjustl(solver)))
case('euler')
  call test_euler()
case('runge-kutta')
  call test_rk()
case('adams-bashforth')
  call test_ab()
case('all')
  call test_euler()
  call test_rk()
  call test_ab()
case default
  print "(A)", 'Error: unknown solver "'//trim(adjustl(solver))//'"'
  print "(A)", 'Valid solver names are:'
  print "(A)", '  + euler'
  print "(A)", '  + runge-kutta'
  print "(A)", '  + adams-bashforth'
  print "(A)", '  + all'
endselect
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine save_plots(title, filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save plots of results.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN) :: title    !< Plot title.
  character(*), intent(IN) :: filename !< Output filename.
  type(pyplot)             :: plt      !< Plot file handler.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call plt%initialize(grid=.true., title=title, legend=.true.)
  call plt%add_plot(x=solution(1, :), y=solution(2, :), label='x-y path', linestyle='r-', linewidth=1)
  call plt%add_plot(x=solution(2, :), y=solution(3, :), label='y-z path', linestyle='b-', linewidth=1)
  call plt%savefig(filename)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_plots

  subroutine test_euler()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit forward Euler ODE solver.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(euler_explicit_integrator) :: euler_integrator !< Euler integrator.
  integer(I_P)                    :: step             !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Lorenz equations by means of explicit Euler solver'
  ! initialize field
  call attractor%init(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
  solution(0, 0) = 0._R_P
  solution(1:space_dimension, 0) = attractor%output()
  ! integrate field
  do step = 1, num_steps
    call euler_integrator%integrate(field=attractor, dt=dt)
    solution(0, step) = step * dt
    solution(1:space_dimension, step) = attractor%output()
  enddo
  call save_plots(title='FOODiE test: Lorenz equation integration, explicit Euler', &
                  filename='lorenz_integration-euler.png')
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_euler

  subroutine test_rk()
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
    ! initialize the RK integrator accordingly to the number of stages used
    rk_integrator = tvd_runge_kutta_integrator(stages=s)
    ! initialize field
    call attractor%init(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = attractor%output()
    ! integrate field
    do step = 1, num_steps
      call rk_integrator%integrate(field=attractor, stage=rk_stage(1:s), dt=dt)
      solution(0, step) = step * dt
      solution(1:space_dimension, step) = attractor%output()
    enddo
    call save_plots(title='FOODiE test: Lorenz equation integration, explicit Runge-Kutta '//trim(str(.true., s))//' stages', &
                    filename='lorenz_integration-rk-'//trim(str(.true., s))//'.png')
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_rk

  subroutine test_ab()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit Adams-Bashforth class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(euler_explicit_integrator)  :: euler_integrator !< Euler integrator.
  type(adams_bashforth_integrator) :: ab_integrator    !< Adams-Bashforth integrator.
  integer, parameter               :: ab_steps=3       !< Adams-Bashforth steps number.
  integer(I_P)                     :: s                !< AB steps counter.
  integer                          :: step             !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Lorenz equations by means of Adams-Bashforth class of solvers'
  do s=1, ab_steps
    print "(A)", ' AB-'//trim(str(.true.,s))
    ! initialize the AB integrator accordingly to the number of time steps used
    call ab_integrator%init(steps=s)
    ! initialize field
    call attractor%init(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta, steps=s)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = attractor%output()
    ! integrate field
    do step = 1, num_steps
      if (s>step) then
        ! the time steps from 1 to s - 1 must be computed with other scheme...
        call euler_integrator%integrate(field=attractor, dt=dt)
      else
        call ab_integrator%integrate(field=attractor, dt=dt)
      endif
      solution(0, step) = step * dt
      solution(1:space_dimension, step) = attractor%output()
    enddo
    call save_plots(title='FOODiE test: Lorenz equation integration, explicit Adams-Bashforth '//trim(str(.true., s))//' steps', &
                    filename='lorenz_integration-ab-'//trim(str(.true., s))//'.png')
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_ab
endprogram integrate_lorenz
