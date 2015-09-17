!< Test FOODiE with the integration of Euler 1D PDEs system.
program integrate_euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Euler 1D PDEs system.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, FR_P, str
use type_euler_1D, only : euler_1D
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
type(Type_Command_Line_Interface) :: cli                      !< Command line interface handler.
type(euler_1D)                    :: domain                   !< Domain of Euler equations.
real(R_P),    parameter           :: CFL=0.7_R_P              !< CFL value.
real(R_P),    parameter           :: t_final=0.15_R_P         !< Final time.
integer(I_P), parameter           :: Ni=100                   !< Number of grid cells.
integer(I_P), parameter           :: Ns=1                     !< Number of differnt initial gas species.
integer(I_P), parameter           :: Nc=Ns+2                  !< Number of conservative variables.
integer(I_P), parameter           :: Np=Ns+4                  !< Number of primitive variables.
real(R_P)                         :: Dx=1._R_P/Ni             !< Space step discretization.
character(3)                      :: BC_L                     !< Left boundary condition type.
character(3)                      :: BC_R                     !< Right boundary condition type.
real(R_P)                         :: cp0(1:Ns)                !< Specific heat at constant pressure.
real(R_P)                         :: cv0(1:Ns)                !< Specific heat at constant volume.
real(R_P)                         :: initial_state(1:Np,1:Ni) !< Initial state of primitive variables.
real(R_P)                         :: x(1:Ni)                  !< Cell center x-abscissa values.
real(R_P), allocatable            :: final_state(:,:)         !< Final state.
integer(I_P)                      :: error                    !< Error handler.
character(99)                     :: solver                   !< Solver used.
logical                           :: plots                    !< Flag for activating plots saving.
logical                           :: results                  !< Flag for activating results saving.
logical                           :: time_serie               !< Flag for activating time serie-results saving.
logical                           :: verbose                  !< Flag for activating more verbose output.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'euler-1D',                                              &
              authors     = 'Fortran-FOSS-Programmers',                              &
              license     = 'GNU GPLv3',                                             &
              description = 'Test FOODiE library on 1D Euler equations integration', &
              examples    = ["euler-1D --solver euler --results",                    &
                             "euler-1D --solver runge-kutta -r ",                    &
                             "euler-1D --solver adams-bashforth",                    &
                             "euler-1D --solver all --plots -r "])
call cli%add(switch='--solver', switch_ab='-s', help='ODE solver used', required=.true., act='store', error=error)
call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--plots', switch_ab='-p', help='Save plots of results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--tserie', switch_ab='-t', help='Save time-serie-results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.', error=error)
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='-s', val=solver, error=error) ; if (error/=0) stop
call cli%get(switch='-r', val=results, error=error) ; if (error/=0) stop
call cli%get(switch='-p', val=plots, error=error) ; if (error/=0) stop
call cli%get(switch='-t', val=time_serie, error=error) ; if (error/=0) stop
call cli%get(switch='--verbose', val=verbose, error=error) ; if (error/=0) stop
! create Euler field initial state
call init()
! integrate Euler equation
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
  integer(I_P) :: i !< Space counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Sod's problem
  BC_L = 'TRA'
  BC_R = 'TRA'
  cp0(1) = 1040._R_P
  cv0(1) = 743._R_P
  do i=1, Ni/2
    x(i) = Dx * i - 0.5_R_P * Dx
    initial_state(:, i) = [1._R_P, & ! rho(s)
                           0._R_P, & ! u
                           1._R_P, & ! p
                           1._R_P, & ! sum(rho(s))
                           cp0/cv0]  ! gamma = cp/cv
  enddo
  do i=Ni/2 + 1, Ni
    x(i) = Dx * i - 0.5_R_P * Dx
    initial_state(:, i) = [0.125_R_P, & ! rho(s)
                           0._R_P,    & ! u
                           0.1_R_P,   & ! p
                           0.125_R_P, & ! sum(rho(s))
                           cp0/cv0]     ! gamma = cp/cv
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  subroutine save_results(title, filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save results.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN) :: title    !< Plot title.
  character(*), intent(IN) :: filename !< Output filename.
  integer(I_P)             :: rawfile  !< Raw file unit for saving results.
  type(pyplot)             :: plt      !< Plot file handler.
  integer(I_P)             :: i        !< Counter.
  integer(I_P)             :: v        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (results.or.plots) final_state = domain%output()
  if (results) then
    open(newunit=rawfile, file=filename//'.dat')
    write(rawfile, '(A)')'# '//title
    write(rawfile, '(A)')'# VARIABLES: "x" "rho(1)" "rho(2)"... "rho(Ns)" "u" "p" "rho" "gamma"'
    do i=1, Ni
      write(rawfile, '('//trim(str(.true.,Np+1))//'('//FR_P//',1X))')x(i), (final_state(v, i), v=1, Np)
    enddo
    close(rawfile)
  endif
  if (plots) then
    call plt%initialize(grid=.true., xlabel='x', title=title)
    do v=1, Ns
      call plt%add_plot(x=x, y=final_state(v, :), label='rho('//trim(str(.true.,v))//')', linestyle='b-', linewidth=1)
    enddo
    call plt%add_plot(x=x, y=final_state(Ns+1, :), label='u', linestyle='r-', linewidth=1)
    call plt%add_plot(x=x, y=final_state(Ns+2, :), label='p', linestyle='g-', linewidth=1)
    call plt%add_plot(x=x, y=final_state(Ns+3, :), label='rho', linestyle='o-', linewidth=1)
    call plt%add_plot(x=x, y=final_state(Ns+4, :), label='gamma', linestyle='c-', linewidth=1)
    call plt%savefig(filename//'.png')
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_results

  subroutine save_time_serie(title, filename, finish, t)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save time-serie results.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN), optional :: title    !< Plot title.
  character(*), intent(IN), optional :: filename !< Output filename.
  logical,      intent(IN), optional :: finish   !< Flag for triggering the file closing.
  real(R_P),    intent(IN)           :: t        !< Current integration time.
  integer(I_P), save                 :: tsfile   !< File unit for saving time serie results.
  integer(I_P)                       :: i        !< Counter.
  integer(I_P)                       :: v        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (time_serie) then
    final_state = domain%output()
    if (present(filename).and.present(title)) then
      open(newunit=tsfile, file=filename)
      write(tsfile, '(A)')'# '//title
    endif
    write(tsfile, '(A)')'# VARIABLES: "x" "rho(1)" "rho(2)"... "rho(Ns)" "u" "p" "rho" "gamma"'
    write(tsfile, '(A)')'# Time: '//str(n=t)
    do i=1, Ni
      write(tsfile, '('//trim(str(.true.,Np+1))//'('//FR_P//',1X))')x(i), (final_state(v, i), v=1, Np)
    enddo
    if (present(finish)) then
      if (finish) close(tsfile)
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_time_serie

  subroutine test_ab()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit Adams-Bashforth class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(tvd_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter               :: rk_stages=5           !< Runge-Kutta stages number.
  type(euler_1D)                   :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(adams_bashforth_integrator) :: ab_integrator         !< Adams-Bashforth integrator.
  integer, parameter               :: ab_steps=3            !< Adams-Bashforth steps number.
  integer                          :: step                  !< Time steps counter.
  real(R_P)                        :: dt                    !< Time step.
  real(R_P)                        :: t                     !< Time.
  integer(I_P)                     :: s                     !< AB steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of Adams-Bashforth class of solvers'
  do s=1, ab_steps
    print "(A)", ' AB-'//trim(str(.true.,s))
    call ab_integrator%init(steps=s)
    call rk_integrator%init(stages=s)
    select case(s)
    case(1)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, steps=s, ord=1)
    case(2, 3)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, steps=s, ord=3)
    case(5)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, steps=s, ord=7)
    endselect
    t = 0._R_P
    call save_time_serie(title='FOODiE test: 1D Euler equations integration, explicit Adams-Bashforth, t='//str(n=t_final)// &
                               trim(str(.true., s))//' steps', &
                         filename='euler_1D_integration-ab-'//trim(str(.true., s))//'-time_serie.dat', &
                         t=t)
    step = 1
    do while(t<t_final)
      if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
      dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=0.1_R_P*CFL)
      if (s>step) then
        ! the time steps from 1 to s - 1 must be computed with other scheme...
        call rk_integrator%integrate(field=domain, stage=rk_stage(1:s), dt=dt)
      else
        call ab_integrator%integrate(field=domain, dt=dt)
      endif
      t = t + dt
      step = step + 1
      call save_time_serie(t=t)
    enddo
    call save_time_serie(t=t, finish=.true.)
    call save_results(title='FOODiE test: 1D Euler equations integration, explicit Adams-Bashforth, t='//str(n=t_final)// &
                            trim(str(.true., s))//' steps', &
                      filename='euler_1D_integration-ab-'//trim(str(.true., s)))
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
  print "(A)", 'Integrating 1D Euler equations by means of explicit Euler solver'
  call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0)
  t = 0._R_P
  call save_time_serie(title='FOODiE test: 1D Euler equations integration, explicit Euler, t='//str(n=t_final), &
                       filename='euler_1D_integration-euler-time_serie.dat', &
                       t=t)
  do while(t<t_final)
    if (verbose) print "(A)", '  Time step: '//str(n=dt)//', Time: '//str(n=t)
    dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=CFL)
    call euler_integrator%integrate(field=domain, dt=dt)
    t = t + dt
    call save_time_serie(t=t)
  enddo
  call save_time_serie(t=t, finish=.true.)
  call save_results(title='FOODiE test: 1D Euler equations integration, explicit Euler, t='//str(n=t_final), &
                    filename='euler_1D_integration-euler')
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
  type(euler_1D)                   :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(euler_1D)                   :: filter                !< Filter displacement.
  type(leapfrog_integrator)        :: lf_integrator         !< Leapfrog integrator.
  integer                          :: step                  !< Time steps counter.
  real(R_P)                        :: dt                    !< Time step.
  real(R_P)                        :: t                     !< Time.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of leapfrog (RAW filtered) class of solvers'
  call lf_integrator%init(nu=1.0_R_P, alpha=0._R_P)
  call rk_integrator%init(stages=rk_stages)
  call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, steps=2, ord=3)
  t = 0._R_P
  call save_time_serie(title='FOODiE test: 1D Euler equations integration, explicit leapfrog, t='//str(n=t_final), &
                       filename='euler_1D_integration-lf-RAW-filter-time_serie.dat', &
                       t=t)
  step = 1
  do while(t<t_final)
    if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
    dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=0.1_R_P*CFL)
    if (2>step) then
      ! the time steps from 1 to s - 1 must be computed with other scheme...
      call rk_integrator%integrate(field=domain, stage=rk_stage, dt=dt)
    else
      call lf_integrator%integrate(field=domain, filter=filter, dt=dt)
    endif
    t = t + dt
    step = step + 1
    call save_time_serie(t=t)
  enddo
  call save_time_serie(t=t, finish=.true.)
  call save_results(title='FOODiE test: 1D Euler equations integration, explicit leapfrog, t='//str(n=t_final), &
                    filename='euler_1D_integration-lf-RAW-filter')
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
  type(euler_1D)                  :: rk_stage(1:registers) !< Runge-Kutta stages.
  real(R_P)                       :: dt                    !< Time step.
  real(R_P)                       :: t                     !< Time.
  integer(I_P)                    :: s                     !< RK stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of low storage (2N) Runge-Kutta class of solvers'
  do s=1, rk_stages
    if (s==2) cycle ! 2 stages not yet implemented
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    call rk_integrator%init(stages=s)
    select case(s)
    case(1)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, ord=1)
    case(2, 3)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, ord=3)
    case(5)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, ord=7)
    endselect
    t = 0._R_P
    call save_time_serie(title='FOODiE test: 1D Euler equations integration, explicit low storage Runge-Kutta, t='//&
                               str(n=t_final)//trim(str(.true., s))//' stages', &
                         filename='euler_1D_integration-lsrk-'//trim(str(.true., s))//'-time_serie.dat', &
                         t=t)
    do while(t<t_final)
      if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
      dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=CFL)
      call rk_integrator%integrate(field=domain, stage=rk_stage, dt=dt)
      t = t + dt
      call save_time_serie(t=t)
    enddo
    call save_time_serie(t=t, finish=.true.)
    call save_results(title='FOODiE test: 1D Euler equations integration, explicit low storage Runge-Kutta, t='//str(n=t_final)// &
                            trim(str(.true., s))//' stages', &
                      filename='euler_1D_integration-lsrk-'//trim(str(.true., s)))
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
  type(euler_1D)                   :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  real(R_P)                        :: dt                    !< Time step.
  real(R_P)                        :: t                     !< Time.
  integer(I_P)                     :: s                     !< RK stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of TVD/SSP Runge-Kutta class of solvers'
  do s=1, rk_stages
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(.true.,s))
    call rk_integrator%init(stages=s)
    select case(s)
    case(1)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, ord=1)
    case(2, 3)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, ord=3)
    case(5)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, ord=7)
    endselect
    t = 0._R_P
    call save_time_serie(title='FOODiE test: 1D Euler equations integration, explicit TVD Runge-Kutta, t='//str(n=t_final)// &
                               trim(str(.true., s))//' stages', &
                         filename='euler_1D_integration-tvdrk-'//trim(str(.true., s))//'-time_serie.dat', &
                         t=t)
    do while(t<t_final)
      if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
      dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=CFL)
      call rk_integrator%integrate(field=domain, stage=rk_stage(1:s), dt=dt)
      t = t + dt
      call save_time_serie(t=t)
    enddo
    call save_time_serie(t=t, finish=.true.)
    call save_results(title='FOODiE test: 1D Euler equations integration, explicit TVD Runge-Kutta, t='//str(n=t_final)// &
                            trim(str(.true., s))//' stages', &
                      filename='euler_1D_integration-tvdrk-'//trim(str(.true., s)))
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_tvd_rk
endprogram integrate_euler_1D
