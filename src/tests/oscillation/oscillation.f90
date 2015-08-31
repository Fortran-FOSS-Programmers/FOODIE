!< Test FOODiE with the integration of Oscillation equations.
program integrate_oscillation
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Oscillation equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, str
use type_oscillation, only : oscillation
use foodie, only : euler_explicit_integrator, tvd_runge_kutta_integrator, adams_bashforth_integrator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(euler_explicit_integrator)  :: euler_integrator                                  !< Euler integrator.
type(tvd_runge_kutta_integrator) :: rk_integrator                                     !< Runge-Kutta integrator.
type(adams_bashforth_integrator) :: ab_integrator                                     !< Adams-Bashforth integrator.
integer, parameter               :: rk_stages=5                                       !< Runge-Kutta stages number.
integer, parameter               :: ab_steps=3                                        !< Adams-Bashforth steps number.
type(oscillation)                :: rk_stage(1:rk_stages)                             !< Runge-Kutta stages.
type(oscillation)                :: attractor                                         !< Oscillation field.
integer                          :: step                                              !< Time steps counter.
integer,   parameter             :: num_steps=1e4                                     !< Maximum time steps.
integer,   parameter             :: space_dimension=2                                 !< Space dimensions.
real(R_P), parameter             :: f=1e-4_R_P                                        !< Frequency.
real(R_P), parameter             :: dt = 100._R_P                                     !< Time step.
real(R_P), parameter             :: initial_state(1:space_dimension)=[0._R_P, 1._R_P] !< Initial state.
real(R_P)                        :: solution(0:space_dimension, 0:num_steps)          !< Solution at each time step.
type(pyplot)                     :: plt                                               !< Plot file handler.
integer(I_P)                     :: s                                                 !< RK stages/AB steps counter.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! solve Oscillation equations by means of the Forward Explicit Euler scheme
! initialize field
attractor = oscillation(initial_state=initial_state, f=f)
solution(0, 0) = 0._R_P
solution(1:space_dimension, 0) = attractor%output()
! integrate field
do step = 1, num_steps
  call euler_integrator%integrate(field=attractor, dt=dt)
  solution(0, step) = step * dt
  solution(1:space_dimension, step) = attractor%output()
enddo
! save plot of results
call plt%initialize(grid=.true., xlabel='time', &
                    title='FOODiE test: Oscillation equation integration, explicit Euler', legend=.true.)
call plt%add_plot(x=solution(0, :), y=solution(1, :), label='u', linestyle='r-', linewidth=1)
call plt%add_plot(x=solution(0, :), y=solution(2, :), label='v', linestyle='g-', linewidth=1)
call plt%savefig('oscillation_integration-euler.png')

! solve Oscillation equations by means of the TVD/SSP RK schemes (testing all available schemes)
do s=1, rk_stages
  if (s==4) cycle ! 4 stages not yet implemented
  ! initialize the RK integrator accordingly to the number of stages used
  rk_integrator = tvd_runge_kutta_integrator(stages=s)
  ! initialize field
  attractor = oscillation(initial_state=initial_state, f=f)
  solution(0, 0) = 0._R_P
  solution(1:space_dimension, 0) = attractor%output()
  ! integrate field
  do step = 1, num_steps
    call rk_integrator%integrate(field=attractor, stage=rk_stage(1:s), dt=dt)
    solution(0, step) = step * dt
    solution(1:space_dimension, step) = attractor%output()
  enddo
  ! save plot of results
  call plt%initialize(grid=.true., xlabel='time', title='FOODiE test: Oscillation equation integration, explicit Runge-Kutta '//&
                      trim(str(.true., s))//' stages', legend=.true.)
  call plt%add_plot(x=solution(0, :), y=solution(1, :), label='u', linestyle='r-', linewidth=1)
  call plt%add_plot(x=solution(0, :), y=solution(2, :), label='v', linestyle='g-', linewidth=1)
  call plt%savefig('oscillation_integration-rk-'//trim(str(.true., s))//'.png')
enddo

! solve Oscillation equations by means of the Adams-Bashforth schemes (testing all available schemes)
do s=1, ab_steps
  ! initialize the AB integrator accordingly to the number of time steps used
  ab_integrator = adams_bashforth_integrator(steps=s)
  ! initialize field
  attractor = oscillation(initial_state=initial_state, f=f, steps=s)
  solution(0, 0) = 0._R_P
  solution(1:space_dimension, 0) = attractor%output()
  ! integrate field
  do step = 1, num_steps
    if (s>step) then
      ! the time steps from 1 to s - 1 must be computed with other scheme...
      call euler_integrator%integrate(field=attractor, dt=dt)
    else
      call ab_integrator%integrate(field=attractor, steps=s, dt=dt)
    endif
    solution(0, step) = step * dt
    solution(1:space_dimension, step) = attractor%output()
    call attractor%update_previous_steps()
  enddo
  ! save plot of results
  call plt%initialize(grid=.true., xlabel='time', title='FOODiE test: Oscillation equation integration, explicit '//&
                     'Adams-Bashforth '//trim(str(.true., s))//' steps', legend=.true.)
  call plt%add_plot(x=solution(0, :), y=solution(1, :), label='u', linestyle='r-', linewidth=1)
  call plt%add_plot(x=solution(0, :), y=solution(2, :), label='v', linestyle='g-', linewidth=1)
  call plt%savefig('oscillation_integration-ab-'//trim(str(.true., s))//'.png')
enddo
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram integrate_oscillation
