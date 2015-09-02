!< Test FOODiE with the integration of Oscillation equations.
program integrate_oscillation
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Oscillation equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, str
use type_oscillation, only : oscillation
use Data_Type_Command_Line_Interface, only : Type_Command_Line_Interface
use foodie, only : euler_explicit_integrator, tvd_runge_kutta_integrator, adams_bashforth_integrator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Command_Line_Interface) :: cli                                               !< Command line interface handler.
type(oscillation)                 :: attractor                                         !< Oscillation field.
integer,   parameter              :: num_steps=1e4                                     !< Maximum time steps.
integer,   parameter              :: space_dimension=2                                 !< Space dimensions.
real(R_P), parameter              :: f=1e-4_R_P                                        !< Frequency.
real(R_P), parameter              :: dt = 100._R_P                                     !< Time step.
real(R_P), parameter              :: initial_state(1:space_dimension)=[0._R_P, 1._R_P] !< Initial state.
real(R_P)                         :: solution(0:space_dimension, 0:num_steps)          !< Solution at each time step.
integer(I_P)                      :: error                                         !< Error handler.
character(99)                     :: solver                                        !< Solver used.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'oscillation',                                              &
              authors     = 'Fortran-FOSS-Programmers',                                 &
              license     = 'GNU GPLv3',                                                &
              description = 'Test FOODiE library on Oscillation equations integration', &
              examples    = ["oscillation --solver euler          ",                    &
                             "oscillation --solver runge-kutta    ",                    &
                             "oscillation --solver adams-bashforth",                    &
                             "oscillation --solver all            "])
call cli%add(switch='--solver', switch_ab='-s', help='ODE solver used', required=.true., act='store', error=error)
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='-s', val=solver, error=error) ; if (error/=0) stop
! integrate Oscillation equations
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
  call plt%initialize(grid=.true., xlabel='time', title=title, legend=.true.)
  call plt%add_plot(x=solution(0, :), y=solution(1, :), label='u', linestyle='r-', linewidth=1)
  call plt%add_plot(x=solution(0, :), y=solution(2, :), label='v', linestyle='b-', linewidth=1)
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
  print "(A)", 'Integrating Oscillation equations by means of explicit Euler solver'
  ! initialize field
  call attractor%init(initial_state=initial_state, f=f)
  solution(0, 0) = 0._R_P
  solution(1:space_dimension, 0) = attractor%output()
  ! integrate field
  do step = 1, num_steps
    call euler_integrator%integrate(field=attractor, dt=dt)
    solution(0, step) = step * dt
    solution(1:space_dimension, step) = attractor%output()
  enddo
  call save_plots(title='FOODiE test: Oscillation equation integrations, explicit Euler', &
                  filename='oscillation_integration-euler.png')
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
  type(oscillation)                :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  integer(I_P)                     :: s                     !< RK stages counter.
  integer(I_P)                     :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Oscillation equations by means of TVD/SSP Runge-Kutta class of solvers'
  do s=1, rk_stages
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    ! initialize the RK integrator accordingly to the number of stages used
    rk_integrator = tvd_runge_kutta_integrator(stages=s)
    ! initialize field
    call attractor%init(initial_state=initial_state, f=f)
    solution(0, 0) = 0._R_P
    solution(1:space_dimension, 0) = attractor%output()
    ! integrate field
    do step = 1, num_steps
      call rk_integrator%integrate(field=attractor, stage=rk_stage(1:s), dt=dt)
      solution(0, step) = step * dt
      solution(1:space_dimension, step) = attractor%output()
    enddo
    call save_plots(title='FOODiE test: Oscillation equation integration, explicit Runge-Kutta '//trim(str(.true., s))//' stages', &
                    filename='oscillation_integration-rk-'//trim(str(.true., s))//'.png')
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
  print "(A)", 'Integrating Oscillation equations by means of Adams-Bashforth class of solvers'
  do s=1, ab_steps
    print "(A)", ' AB-'//trim(str(.true.,s))
    ! initialize the AB integrator accordingly to the number of time steps used
    call ab_integrator%init(steps=s)
    ! initialize field
    call attractor%init(initial_state=initial_state, f=f, steps=s)
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
    call save_plots(title='FOODiE test: Oscillation equation integration, explicit '//'Adams-Bashforth '//&
                          trim(str(.true., s))//' steps', &
                    filename='oscillation_integration-ab-'//trim(str(.true., s))//'.png')
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_ab
endprogram integrate_oscillation
