!< Test FOODiE with the integration of Burgers equation.
program integrate_burgers
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Burgers equation.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, str
use type_burgers, only : burgers
use Data_Type_Command_Line_Interface, only : Type_Command_Line_Interface
use foodie, only : euler_explicit_integrator, tvd_runge_kutta_integrator, adams_bashforth_integrator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Command_Line_Interface) :: cli                 !< Command line interface handler.
type(burgers)                     :: domain              !< Burgers field domain.
real(R_P),    parameter           :: CFL=0.1_R_P         !< CFL value.
real(R_P),    parameter           :: t_final=0.6_R_P     !< Final time.
real(R_P),    parameter           :: nu=1._R_P           !< Viscosity.
integer(I_P), parameter           :: Ni=100              !< Number of grid nodes.
real(R_P)                         :: h                   !< Space step discretization.
real(R_P)                         :: initial_state(1:Ni) !< Initial state.
real(R_P)                         :: x(1:Ni)             !< Nodes values.
  real(R_P), allocatable          :: final_state(:)      !< Final state.
integer(I_P)                      :: error               !< Error handler.
character(99)                     :: solver              !< Solver used.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'burgers',                                             &
              authors     = 'Fortran-FOSS-Programmers',                            &
              license     = 'GNU GPLv3',                                           &
              description = 'Test FOODiE library on Burgers equation integration', &
              examples    = ["burgers --solver euler          ",                   &
                             "burgers --solver runge-kutta    ",                   &
                             "burgers --solver adams-bashforth",                   &
                             "burgers --solver all            "])
call cli%add(switch='--solver', switch_ab='-s', help='ODE solver used', required=.true., act='store', error=error)
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='-s', val=solver, error=error) ; if (error/=0) stop
! create Burgers field initial state
call init()
! integrate Burgers equation
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
  subroutine init()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Initialize the field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), parameter :: pi=4._R_P * atan(1._R_P)  !< Pi greek.
  integer(I_P)         :: i                         !< Space counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  h = 2._R_P * pi / Ni
  do i=1, Ni
    x(i) = h * (i - 1)
    initial_state(i) = 10._R_P * sin(x(i))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  subroutine save_plots(title, filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save plots of results.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN) :: title    !< Plot title.
  character(*), intent(IN) :: filename !< Output filename.
  type(pyplot)             :: plt      !< Plot file handler.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call plt%initialize(grid=.true., xlabel='x', title=title)
  call plt%add_plot(x=x, y=final_state, label='U', linestyle='b-', linewidth=1)
  call plt%savefig(filename)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_plots

  subroutine test_euler()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit forward Euler ODE solver.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(euler_explicit_integrator) :: euler_integrator !< Euler integrator.
  real(R_P)                       :: dt               !< Time step.
  real(R_P)                       :: t                !< Time.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Burgers equation by means of explicit Euler solver'
  ! initialize field
  call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu)
  ! integrate field
  dt = domain%dt(CFL=CFL)
  t = 0._R_P
  do while(t<t_final)
    call euler_integrator%integrate(field=domain, dt=dt)
    t = t + dt
  enddo
  final_state = domain%output()
  call save_plots(title='FOODiE test: Burgers equation integration at t=0.6, explicit Euler', &
                  filename='burgers_integration-euler.png')
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
  type(burgers)                    :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  real(R_P)                        :: dt                    !< Time step.
  real(R_P)                        :: t                     !< Time.
  integer(I_P)                     :: s                     !< RK stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Burgers equation by means of TVD/SSP Runge-Kutta class of solvers'
  do s=1, rk_stages
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    ! initialize the RK integrator accordingly to the number of stages used
    rk_integrator = tvd_runge_kutta_integrator(stages=s)
    ! initialize field
    call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu)
    ! integrate field
    dt = domain%dt(CFL=CFL)
    t = 0._R_P
    do while(t<t_final)
      call rk_integrator%integrate(field=domain, stage=rk_stage(1:s), dt=dt)
      t = t + dt
    enddo
    final_state = domain%output()
    call save_plots(title='FOODiE test: Burgers equation integration at t=0.6, explicit Runge-Kutta '//&
                          trim(str(.true., s))//' stages', &
                    filename='burgers_integration-rk-'//trim(str(.true., s))//'.png')
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
  integer                          :: step             !< Time steps counter.
  real(R_P)                        :: dt               !< Time step.
  real(R_P)                        :: t                !< Time.
  integer(I_P)                     :: s                !< AB steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Burgers equation by means of Adams-Bashforth class of solvers'
  do s=1, ab_steps
    print "(A)", ' AB-'//trim(str(.true.,s))
    ! initialize the AB integrator accordingly to the number of time steps used
    call ab_integrator%init(steps=s)
    ! initialize field
    call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu, steps=s)
    ! integrate field
    dt = domain%dt(CFL=CFL)
    t = 0._R_P
    step = 1
    do while(t<t_final)
      if (s>step) then
        ! the time steps from 1 to s - 1 must be computed with other scheme...
        call euler_integrator%integrate(field=domain, dt=dt)
      else
        call ab_integrator%integrate(field=domain, dt=dt)
      endif
      t = t + dt
      step = step + 1
    enddo
    final_state = domain%output()
    call save_plots(title='FOODiE test: Burgers equation integration at t=0.6, explicit '//&
                          'Adams-Bashforth '//trim(str(.true., s))//' steps', &
                    filename='burgers_integration-ab-'//trim(str(.true., s))//'.png')
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_ab
endprogram integrate_burgers
