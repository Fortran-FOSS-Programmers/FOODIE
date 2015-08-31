!< Test FOODiE with the integration of Lorenz equations.
program integrate_lorenz
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Lorenz equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, str
use type_lorenz, only : lorenz
use foodie, only : euler_explicit_integrator, tvd_runge_kutta_integrator, adams_bashforth_integrator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(euler_explicit_integrator)  :: euler_integrator                                !< Euler integrator.
type(tvd_runge_kutta_integrator) :: rk_integrator                                   !< Runge-Kutta integrator.
type(adams_bashforth_integrator) :: ab_integrator                                   !< Adams-Bashforth integrator.
integer, parameter               :: rk_stages=5                                     !< Runge-Kutta stages number.
integer, parameter               :: ab_steps=3                                      !< Adams-Bashforth steps number.
type(lorenz)                     :: rk_stage(1:rk_stages)                           !< Runge-Kutta stages.
type(lorenz)                     :: attractor                                       !< Lorenz field.
integer                          :: step                                            !< Time steps counter.
integer,   parameter             :: num_steps=2000                                  !< Maximum time steps.
integer,   parameter             :: space_dimension=3                               !< Space dimensions.
real(R_P), parameter             :: sigma=10._R_P                                   !< Lorenz' \(\sigma\).
real(R_P), parameter             :: rho=28._R_P                                     !< Lorenz' \(\rho\).
real(R_P), parameter             :: beta=8._R_P/3._R_P                              !< Lorenz' \(\beta\).
real(R_P), parameter             :: dt=0.01_R_P                                     !< Time step.
real(R_P), parameter             :: initial_state(1:space_dimension)=[1., 1., 1.]   !< Initial state.
real(R_P)                        :: solution(0:space_dimension, 0:num_steps)        !< Solution at each time step.
type(pyplot)                     :: plt                                             !< Plot file handler.
integer(I_P)                     :: s                                               !< RK stages/AB steps counter.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! solve Lorenz equations by means of the Forward Explicit Euler scheme
! initialize field
attractor = lorenz(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
solution(0, 0) = 0._R_P
solution(1:space_dimension, 0) = attractor%output()
! integrate field
do step = 1, num_steps
  call euler_integrator%integrate(field=attractor, dt=dt)
  solution(0, step) = step * dt
  solution(1:space_dimension, step) = attractor%output()
enddo
! save plot of results
call plt%initialize(grid=.true., xlabel='time', title='FOODiE test: Lorenz equation integration, explicit Euler', legend=.true.)
call plt%add_plot(x=solution(0, :), y=solution(1, :), label='sigma', linestyle='g-', linewidth=1, markersize=4)
call plt%add_plot(x=solution(0, :), y=solution(2, :), label='rho', linestyle='b-x', linewidth=1, markersize=4)
call plt%add_plot(x=solution(0, :), y=solution(2, :), label='beta', linestyle='r-o', linewidth=1, markersize=4)
call plt%savefig('lorenz_integration-euler.png')

! solve Lorenz equations by means of the TVD/SSP RK schemes (testing all available schemes)
do s=1, rk_stages
  if (s==4) cycle ! 4 stages not yet implemented
  ! initialize the RK integrator accordingly to the number of stages used
  rk_integrator = tvd_runge_kutta_integrator(stages=s)
  ! initialize field
  attractor = lorenz(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
  solution(0, 0) = 0._R_P
  solution(1:space_dimension, 0) = attractor%output()
  ! integrate field
  do step = 1, num_steps
    call rk_integrator%integrate(field=attractor, stage=rk_stage(1:s), dt=dt)
    solution(0, step) = step * dt
    solution(1:space_dimension, step) = attractor%output()
  enddo
  ! save plot of results
  call plt%initialize(grid=.true., xlabel='time', title='FOODiE test: Lorenz equation integration, explicit Runge-Kutta '//&
                      trim(str(.true., s))//' stages', legend=.true.)
  call plt%add_plot(x=solution(0, :), y=solution(1, :), label='sigma', linestyle='g-', linewidth=1, markersize=4)
  call plt%add_plot(x=solution(0, :), y=solution(2, :), label='rho', linestyle='b-x', linewidth=1, markersize=4)
  call plt%add_plot(x=solution(0, :), y=solution(2, :), label='beta', linestyle='r-o', linewidth=1, markersize=4)
  call plt%savefig('lorenz_integration-rk-'//trim(str(.true., s))//'.png')
enddo

! solve Lorenz equations by means of the Adams-Bashforth schemes (testing all available schemes)
do s=1, ab_steps
  ! initialize the AB integrator accordingly to the number of time steps used
  ab_integrator = adams_bashforth_integrator(steps=s)
  ! initialize field
  attractor = lorenz(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta, steps=s)
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
  call plt%initialize(grid=.true., xlabel='time', title='FOODiE test: Lorenz equation integration, explicit Adams-Bashforth '//&
                      trim(str(.true., s))//' steps', legend=.true.)
  call plt%add_plot(x=solution(0, :), y=solution(1, :), label='sigma', linestyle='g-', linewidth=1, markersize=4)
  call plt%add_plot(x=solution(0, :), y=solution(2, :), label='rho', linestyle='b-x', linewidth=1, markersize=4)
  call plt%add_plot(x=solution(0, :), y=solution(2, :), label='beta', linestyle='r-o', linewidth=1, markersize=4)
  call plt%savefig('lorenz_integration-ab-'//trim(str(.true., s))//'.png')
enddo
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram integrate_lorenz
