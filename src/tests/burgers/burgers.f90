!< Test FOODiE with the integration of Burgers equation.
program integrate_burgers
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Burgers equation.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, str
use type_burgers, only : burgers
use foodie, only : euler_explicit_integrator, tvd_runge_kutta_integrator, adams_bashforth_integrator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(euler_explicit_integrator)  :: euler_integrator            !< Euler integrator.
type(tvd_runge_kutta_integrator) :: rk_integrator               !< Runge-Kutta integrator.
type(adams_bashforth_integrator) :: ab_integrator               !< Adams-Bashforth integrator.
integer, parameter               :: rk_stages=5                 !< Runge-Kutta stages number.
integer, parameter               :: ab_steps=3                                      !< Adams-Bashforth steps number.
type(burgers)                    :: rk_stage(1:rk_stages)       !< Runge-Kutta stages.
real(R_P), parameter             :: pi = 4._R_P * atan(1._R_P)  !< Pi greek.
type(burgers)                    :: domain                      !< Burgers field domain.
integer                          :: step                        !< Time steps counter.
real(R_P)                        :: dt                          !< Time step.
real(R_P)                        :: t=0.                        !< Initial time.
real(R_P)                        :: t_final=0.6                 !< Final time.
integer,   parameter             :: Ni=100                      !< Number of grid nodes.
real(R_P), parameter             :: h=2._R_P * pi / Ni          !< Space step discretization.
real(R_P), parameter             :: nu=1._R_P                   !< Viscosity.
real(R_P)                        :: initial_state(1:Ni)         !< Initial state.
real(R_P)                        :: x(1:Ni)                     !< Nodes values.
real(R_P), allocatable           :: final_state(:)              !< Final state.
type(pyplot)                     :: plt                         !< Plot file handler.
integer(I_P)                     :: s                           !< RK stages/AB steps counter.
integer(I_P)                     :: i                           !< Space counter.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! initial state
do i=1, Ni
  x(i) = h * (i - 1)
  initial_state(i) = 10._R_P * sin(x(i))
enddo

! solve Burgers equations by means of the Forward Explicit Euler scheme
! initialize field
call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu)
! integrate field
dt = domain%dt(CFL=0.1_R_P)
do while(t<t_final)
  call euler_integrator%integrate(field=domain, dt=dt)
  t = t + dt
enddo
final_state = domain%output()
! save plot of results
call plt%initialize(grid=.true., xlabel='x', title='FOODiE test: Burgers equation integration at t=0.6, explicit Euler')
call plt%add_plot(x=x, y=final_state, label='U', linestyle='g-', linewidth=1, markersize=4)
call plt%savefig('burgers_integration-euler.png')

! solve Burgers equations by means of the TVD/SSP RK schemes (testing all available schemes)
do s=1, rk_stages
  if (s==4) cycle ! 4 stages not yet implemented
  ! initialize the RK integrator accordingly to the number of stages used
  rk_integrator = tvd_runge_kutta_integrator(stages=s)
  ! initialize field
  call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu)
  ! integrate field
  dt = domain%dt(CFL=0.1_R_P)
  do while(t<t_final)
    call rk_integrator%integrate(field=domain, stage=rk_stage(1:s), dt=dt)
    t = t + dt
  enddo
  final_state = domain%output()
  ! save plot of results
  call plt%initialize(grid=.true., xlabel='x', title='FOODiE test: Burgers equation integration at t=0.6, explicit Runge-Kutta '//&
                      trim(str(.true., s))//' stages')
  call plt%add_plot(x=x, y=final_state, label='U', linestyle='g-', linewidth=1, markersize=4)
  call plt%savefig('burgers_integration-rk-'//trim(str(.true., s))//'.png')
enddo

! solve Burgers equations by means of the Adams-Bashforth schemes (testing all available schemes)
do s=1, ab_steps
  ! initialize the AB integrator accordingly to the number of time steps used
  call ab_integrator%init(steps=s)
  ! initialize field
  call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu, steps=s)
  ! integrate field
  dt = domain%dt(CFL=0.1_R_P)
  step = 1
  do while(t<t_final)
    if (s>step) then
      ! the time steps from 1 to s - 1 must be computed with other scheme...
      call euler_integrator%integrate(field=domain, dt=dt)
    else
      call ab_integrator%integrate(field=domain, steps=s, dt=dt)
    endif
    t = t + dt
    step = step + 1
    call domain%update_previous_steps()
  enddo
  final_state = domain%output()
  ! save plot of results
  call plt%initialize(grid=.true., xlabel='x', title='FOODiE test: Burgers equation integration at t=0.6, explicit '//&
                      'Adams-Bashforth '//trim(str(.true., s))//' steps')
  call plt%add_plot(x=x, y=final_state, label='U', linestyle='g-', linewidth=1, markersize=4)
  call plt%savefig('burgers_integration-ab-'//trim(str(.true., s))//'.png')
enddo
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram integrate_burgers
