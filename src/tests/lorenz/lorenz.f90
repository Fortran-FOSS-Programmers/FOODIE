!< Test FOODiE with the integration of Lorenz equations.
program integrate_lorenz
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Lorenz equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P
use type_lorenz, only : lorenz
use foodie, only : euler_explicit_integrator, tvd_runge_kutta_integrator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(euler_explicit_integrator)  :: euler_integrator                                !< Euler integrator.
type(tvd_runge_kutta_integrator) :: rk_integrator                                   !< Runge-Kutta integrator.
integer, parameter               :: rk_stages=3                                     !< Runge-Kutta stages number.
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
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! solve Lorenz equations by means of the Forward Explicit Euler scheme
! initialize field
attractor = lorenz(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
solution(0, 0) = 0.
solution(1:3, 0) = attractor%output()
! integrate field
do step = 1, num_steps
  call euler_integrator%integrate(field=attractor, dt=dt)
  solution(0,   step) = step * dt
  solution(1:3, step) = attractor%output()
enddo
! save plot of results
call plt%initialize(grid=.true., xlabel='time', title='FOODiE test: Lorenz equation integration, explicit Euler', legend=.true.)
call plt%add_plot(x=solution(0, :), y=solution(1, :), label='sigma', linestyle='g-', linewidth=1, markersize=4)
call plt%add_plot(x=solution(0, :), y=solution(2, :), label='rho', linestyle='b-x', linewidth=1, markersize=4)
call plt%add_plot(x=solution(0, :), y=solution(2, :), label='beta', linestyle='r-o', linewidth=1, markersize=4)
call plt%savefig('lorenz_integration-euler.png')

! solve Lorenz equations by means of the SSPRK(3,3) scheme
! initialize the RK integrator accordingly to the number of stages used
rk_integrator = tvd_runge_kutta_integrator(stages=rk_stages)
! initialize field
attractor = lorenz(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
solution(0, 0) = 0.
solution(1:3, 0) = attractor%output()
! integrate field
do step = 1, num_steps
  call rk_integrator%integrate(field=attractor, stage=rk_stage, dt=dt)
  solution(0,   step) = step * dt
  solution(1:3, step) = attractor%output()
enddo
! save plot of results
call plt%initialize(grid=.true., xlabel='time', title='FOODiE test: Lorenz equation integration, explicit Runge-Kutta 1', &
                    legend=.true.)
call plt%add_plot(x=solution(0, :), y=solution(1, :), label='sigma', linestyle='g-', linewidth=1, markersize=4)
call plt%add_plot(x=solution(0, :), y=solution(2, :), label='rho', linestyle='b-x', linewidth=1, markersize=4)
call plt%add_plot(x=solution(0, :), y=solution(2, :), label='beta', linestyle='r-o', linewidth=1, markersize=4)
call plt%savefig('lorenz_integration-rk-3.png')
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram integrate_lorenz
