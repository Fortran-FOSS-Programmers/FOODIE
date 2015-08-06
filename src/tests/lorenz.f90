!< Test FOODiE with the integration of Lorenz equations.
program integrate_lorenz
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Lorenz equations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P
use type_lorenz, only : lorenz
use foodie, only : euler_explicit
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(lorenz)           :: attractor
integer                :: step
integer,parameter      :: num_steps=2000
integer,parameter      :: space_dimension=3
real(R_P),   parameter :: sigma=10.
real(R_P),   parameter :: rho=28.
real(R_P),   parameter :: beta=8./3.
real(R_P),   parameter :: dt=0.01
real(R_P),   parameter :: initial_state(1:space_dimension)=[1., 1., 1.]
real(R_P)              :: solution(0:space_dimension, 0:num_steps)
type(pyplot)           :: plt
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
attractor = lorenz(initial_state=initial_state, sigma=sigma, rho=rho, beta=beta)
solution(0, 0) = 0.
solution(1:3, 0) = attractor%output()
do step = 1, num_steps
  call euler_explicit(field=attractor, dt=dt)
  solution(0,   step) = step * dt
  solution(1:3, step) = attractor%output()
enddo
call plt%initialize(grid=.true., xlabel='time', title='FOODiE test: Lorenz equation integration, explicit Euler', legend=.true.)
call plt%add_plot(x=solution(0, :), y=solution(1, :), label='sigma', linestyle='g-', linewidth=1, markersize=4)
call plt%add_plot(x=solution(0, :), y=solution(2, :), label='rho', linestyle='b-x', linewidth=1, markersize=4)
call plt%add_plot(x=solution(0, :), y=solution(2, :), label='beta', linestyle='r-o', linewidth=1, markersize=4)
call plt%savefig('lorenz_integration.png')
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram integrate_lorenz
