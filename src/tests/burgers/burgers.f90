!< Test FOODiE with the integration of Burgers equation.
program integrate_burgers
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Burgers equation.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, FR_P, str
use type_burgers, only : burgers
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
type(Type_Command_Line_Interface) :: cli                 !< Command line interface handler.
type(burgers)                     :: domain              !< Burgers field domain.
real(R_P),    parameter           :: CFL=0.1_R_P         !< CFL value.
real(R_P),    parameter           :: t_final=0.6_R_P     !< Final time.
real(R_P),    parameter           :: nu=1._R_P           !< Viscosity.
integer(I_P), parameter           :: Ni=100              !< Number of grid nodes.
real(R_P)                         :: h                   !< Space step discretization.
real(R_P)                         :: initial_state(1:Ni) !< Initial state.
real(R_P)                         :: x(1:Ni)             !< Nodes values.
real(R_P), allocatable            :: final_state(:)      !< Final state.
integer(I_P)                      :: error               !< Error handler.
character(99)                     :: solver              !< Solver used.
logical                           :: plots               !< Flag for activating plots saving.
logical                           :: results             !< Flag for activating results saving.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'burgers',                                             &
              authors     = 'Fortran-FOSS-Programmers',                            &
              license     = 'GNU GPLv3',                                           &
              description = 'Test FOODiE library on Burgers equation integration', &
              examples    = ["burgers --solver euler --results  ",                 &
                             "burgers --solver ls-runge-kutta -r",                 &
                             "burgers --solver adams-bashforth  ",                 &
                             "burgers --solver all --plots -r   "])
call cli%add(switch='--solver', switch_ab='-s', help='ODE solver used', required=.true., act='store')
call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.')
call cli%add(switch='--plots', switch_ab='-p', help='Save plots of results', required=.false., act='store_true', def='.false.')
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='-s', val=solver, error=error) ; if (error/=0) stop
call cli%get(switch='-r', val=results, error=error) ; if (error/=0) stop
call cli%get(switch='-p', val=plots, error=error) ; if (error/=0) stop
! create Burgers field initial state
call init()
! integrate Burgers equation
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

  subroutine save_results(title, filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save plots of results.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN) :: title    !< Plot title.
  character(*), intent(IN) :: filename !< Output filename.
  integer(I_P)             :: rawfile  !< Raw file unit for saving results.
  type(pyplot)             :: plt      !< Plot file handler.
  integer(I_P)             :: i        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (results) then
    open(newunit=rawfile, file=filename//'.dat')
    write(rawfile, '(A)')'# '//title
    write(rawfile, '(A)')'# VARIABLES: "x" "U"'
    do i=1, Ni
      write(rawfile, '(2('//FR_P//',1X))')x(i), final_state(i)
    enddo
    close(rawfile)
  endif
  if (plots) then
    call plt%initialize(grid=.true., xlabel='x', title=title)
    call plt%add_plot(x=x, y=final_state, label='U', linestyle='b-', linewidth=1)
    call plt%savefig(filename//'.png')
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
  type(burgers)                    :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(adams_bashforth_integrator) :: ab_integrator         !< Adams-Bashforth integrator.
  integer, parameter               :: ab_steps=3            !< Adams-Bashforth steps number.
  integer                          :: step                  !< Time steps counter.
  real(R_P)                        :: dt                    !< Time step.
  real(R_P)                        :: t(1:ab_steps)         !< Times.
  integer(I_P)                     :: s                     !< AB steps counter.
  integer(I_P)                     :: ss                    !< AB substeps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Burgers equation by means of Adams-Bashforth class of solvers'
  do s=1, ab_steps
    print "(A)", ' AB-'//trim(str(.true.,s))
    call ab_integrator%init(steps=s)
    call rk_integrator%init(stages=s)
    call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu, steps=s)
    dt = domain%dt(CFL=CFL)
    t = 0._R_P
    step = 1
    do while(t(s)<t_final)
      if (s>=step) then
        ! the time steps from 1 to s - 1 must be computed with other scheme...
        call rk_integrator%integrate(field=domain, stage=rk_stage(1:s), dt=dt, t=t(s))
        if (step>1) then
          t(step) = t(step-1) + dt
        else
          t(step) = dt
        endif
      else
        call ab_integrator%integrate(field=domain, dt=dt, t=t)
        do ss=1, s-1
          t(ss) = t(ss + 1)
        enddo
        t(s) = t(s) + dt
      endif
      step = step + 1
    enddo
    final_state = domain%output()
    call save_results(title='FOODiE test: Burgers equation integration, t='//str(n=t_final)//' explicit '//&
                            'Adams-Bashforth '//trim(str(.true., s))//' steps', &
                      filename='burgers_integration-ab-'//trim(str(.true., s)))
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
  real(R_P)                       :: dt               !< Time step.
  real(R_P)                       :: t                !< Time.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Burgers equation by means of explicit Euler solver'
  call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu)
  dt = domain%dt(CFL=CFL)
  t = 0._R_P
  do while(t<t_final)
    call euler_integrator%integrate(field=domain, dt=dt, t=t)
    t = t + dt
  enddo
  final_state = domain%output()
  call save_results(title='FOODiE test: Burgers equation integration, t='//str(n=t_final)//' explicit Euler', &
                    filename='burgers_integration-euler')
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
  type(burgers)                    :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(burgers)                    :: filter                !< Filter displacement.
  type(leapfrog_integrator)        :: lf_integrator         !< Leapfrog integrator.
  real(R_P)                        :: dt                    !< Time step.
  real(R_P)                        :: t                     !< Time.
  integer                          :: step                  !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Burgers equation by means of leapfrog (RAW filtered) class of solvers'
  call lf_integrator%init(nu=1.0_R_P, alpha=0._R_P)
  call rk_integrator%init(stages=rk_stages)
  call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu, steps=2)
  dt = domain%dt(CFL=CFL)
  t = 0._R_P
  step = 1
  do while(t<t_final)
    if (2>=step) then
      ! the time steps from 1 to 2 must be computed with other scheme...
      call rk_integrator%integrate(field=domain, stage=rk_stage, dt=dt, t=t)
    else
      call lf_integrator%integrate(field=domain, filter=filter, dt=dt, t=t)
    endif
    t = t + dt
    step = step + 1
  enddo
  final_state = domain%output()
  call save_results(title='FOODiE test: Burgers equation integration, t='//str(n=t_final)//' explicit leapfrog scheme',&
                    filename='burgers_integration-lf-RAW-filter')
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
  type(burgers)                   :: rk_stage(1:registers) !< Runge-Kutta stages.
  real(R_P)                       :: dt                    !< Time step.
  real(R_P)                       :: t                     !< Time.
  integer(I_P)                    :: s                     !< RK stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating Burgers equation by means of low storage (2N) Runge-Kutta class of solvers'
  do s=1, rk_stages
    if (s==2) cycle ! 2 stages not yet implemented
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    call rk_integrator%init(stages=s)
    call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu)
    dt = domain%dt(CFL=CFL)
    t = 0._R_P
    do while(t<t_final)
      call rk_integrator%integrate(field=domain, stage=rk_stage, dt=dt, t=t)
      t = t + dt
    enddo
    final_state = domain%output()
    call save_results(title='FOODiE test: Burgers equation integration, t='//str(n=t_final)//' explicit low storage Runge-Kutta '//&
                            trim(str(.true., s))//' stages', &
                      filename='burgers_integration-lsrk-'//trim(str(.true., s)))
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
    call rk_integrator%init(stages=s)
    call domain%init(initial_state=initial_state, Ni=Ni, h=h, nu=nu)
    dt = domain%dt(CFL=CFL)
    t = 0._R_P
    do while(t<t_final)
      call rk_integrator%integrate(field=domain, stage=rk_stage(1:s), dt=dt, t=t)
      t = t + dt
    enddo
    final_state = domain%output()
    call save_results(title='FOODiE test: Burgers equation integration, t='//str(n=t_final)//' explicit TVD Runge-Kutta '//&
                            trim(str(.true., s))//' stages', &
                      filename='burgers_integration-tvdrk-'//trim(str(.true., s)))
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_tvd_rk
endprogram integrate_burgers
