!< Test FOODIE with the integration of Euler 1D PDEs system.
program integrate_euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODIE with the integration of Euler 1D PDEs system.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use flap, only : command_line_interface
use foodie, only : integrator_adams_bashforth, &
                   euler_explicit_integrator, &
                   leapfrog_integrator, &
                   ls_runge_kutta_integrator, &
                   tvd_runge_kutta_integrator
use penf, only : R_P, I_P, FR_P, str
use pyplot_module, only :  pyplot
use type_euler_1D, only : euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(command_line_interface)  :: cli                !< Command line interface handler.
type(euler_1D)                :: domain             !< Domain of Euler equations.
integer(I_P)                  :: Ns                 !< Number of differnt initial gas species.
integer(I_P)                  :: Nc                 !< Number of conservative variables, Nc=Ns+2.
integer(I_P)                  :: Np                 !< Number of primitive variables, Np=Ns+4.
integer(I_P)                  :: Ni                 !< Number of grid cells.
real(R_P)                     :: Dx                 !< Space step discretization.
real(R_P)                     :: CFL                !< CFL value.
real(R_P)                     :: t_final            !< Final time.
character(3)                  :: BC_L               !< Left boundary condition type.
character(3)                  :: BC_R               !< Right boundary condition type.
real(R_P), allocatable        :: cp0(:)             !< Specific heat at constant pressure [1:Ns].
real(R_P), allocatable        :: cv0(:)             !< Specific heat at constant volume [1:Ns].
real(R_P), allocatable        :: initial_state(:,:) !< Initial state of primitive variables [1:Np,1:Ni].
real(R_P), allocatable        :: xcenter(:)         !< Cell center x-abscissa values, [1:Ni].
real(R_P), allocatable        :: xnode(:)           !< Cell node x-abscissa values, [0:Ni].
real(R_P), allocatable        :: av_xnode(:)        !< Average-grid cell node x-abscissa values, [0:Ni].
real(R_P), allocatable        :: final_state(:,:)   !< Final state.
real(R_P), allocatable        :: av_state(:,:)      !< Average-grid final state.
character(len=:), allocatable :: variables          !< Variables names list.
character(len=:), allocatable :: output             !< Output files basename.
integer(I_P)                  :: error              !< Error handler.
integer(I_P)                  :: stages_steps       !< Number of stages/steps used.
character(99)                 :: solver             !< Solver used.
character(99)                 :: problem            !< Problem solved.
character(99)                 :: output_cli         !< Output files basename.
logical                       :: plots              !< Flag for activating plots saving.
logical                       :: results            !< Flag for activating results saving.
logical                       :: time_serie         !< Flag for activating time serie-results saving.
logical                       :: verbose            !< Flag for activating more verbose output.
integer(I_P)                  :: av_Ni              !< Average the solution over an average-grid.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'euler-1D',                                              &
              authors     = 'Fortran-FOSS-Programmers',                              &
              license     = 'GNU GPLv3',                                             &
              description = 'Test FOODIE library on 1D Euler equations integration', &
              examples    = ["euler-1D --solver euler --results  ",                  &
                             "euler-1D --solver ls-runge-kutta -r",                  &
                             "euler-1D --solver adams-bashforth  ",                  &
                             "euler-1D --solver all --plots -r   "])
call cli%add(switch='--solver', switch_ab='-s', help='ODE solver', required=.true., act='store')
call cli%add(switch='--problem', help='Problem solved', required=.false., def='sod', act='store', choices='sod,smooth')
call cli%add(switch='--Ni', help='Number finite volumes', required=.false., act='store', def='100')
call cli%add(switch='--av_Ni', help='Number finite volumes over average the solution', required=.false., act='store', def='-1')
call cli%add(switch='--ss', help='Stages/steps used', required=.false., act='store', def='-1')
call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.')
call cli%add(switch='--plots', switch_ab='-p', help='Save plots of results', required=.false., act='store_true', def='.false.')
call cli%add(switch='--tserie', switch_ab='-t', help='Save time-serie-results', required=.false., act='store_true', def='.false.')
call cli%add(switch='--output', help='Output files basename', required=.false., act='store', def='unset')
call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.')
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='--solver', val=solver, error=error) ; if (error/=0) stop
call cli%get(switch='--problem', val=problem, error=error) ; if (error/=0) stop
call cli%get(switch='--Ni', val=Ni, error=error) ; if (error/=0) stop
call cli%get(switch='--av_Ni', val=av_Ni, error=error) ; if (error/=0) stop
call cli%get(switch='--ss', val=stages_steps, error=error) ; if (error/=0) stop
call cli%get(switch='--results', val=results, error=error) ; if (error/=0) stop
call cli%get(switch='--plots', val=plots, error=error) ; if (error/=0) stop
call cli%get(switch='--tserie', val=time_serie, error=error) ; if (error/=0) stop
call cli%get(switch='--output', val=output_cli, error=error) ; if (error/=0) stop
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
  real(R_P), parameter :: pi=4._R_P * atan(1._R_P) !< Pi greek.
  integer(I_P)         :: i                        !< Space counter.
  integer(I_P)         :: s                        !< Species counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(xcenter(1:Ni))
  allocate(xnode(0:Ni))
  Dx=1._R_P/Ni
  do i=1, Ni
    xcenter(i) = Dx * i - 0.5_R_P * Dx
  enddo
  do i=0, Ni
    xnode(i) = Dx * i
  enddo
  if ((av_Ni>0).and.(av_Ni/=Ni)) then
    allocate(av_xnode(0:av_Ni))
    do i=0, av_Ni
      av_xnode(i) = 1._R_P/av_Ni * i
    enddo
  endif
  select case(trim(adjustl(problem)))
  case('sod')
    print "(A)", 'Solving "'//trim(adjustl(problem))//'" problem'
    t_final = 0.2_R_P
    CFL = 0.7_R_P
    Ns = 1
    Nc = Ns + 2
    Np = Ns + 4
    allocate(initial_state(1:Np, 1:Ni))
    allocate(cp0(1:Ns))
    allocate(cv0(1:Ns))
    variables = 'VARIABLES="x"'
    do s=1, Ns
      variables = variables//' "rho('//trim(str(s,.true.))//')"'
    enddo
    variables = variables//'  "u" "p" "rho" "gamma"'
    BC_L = 'TRA'
    BC_R = 'TRA'
    cp0(1) = 1040._R_P
    cv0(1) = 743._R_P
    do i=1, Ni/2
      initial_state(:, i) = [1._R_P, & ! rho(s)
                             0._R_P, & ! u
                             1._R_P, & ! p
                             1._R_P, & ! sum(rho(s))
                             cp0/cv0]  ! gamma = cp/cv
    enddo
    do i=Ni/2 + 1, Ni
      initial_state(:, i) = [0.125_R_P, & ! rho(s)
                             0._R_P,    & ! u
                             0.1_R_P,   & ! p
                             0.125_R_P, & ! sum(rho(s))
                             cp0/cv0]     ! gamma = cp/cv
    enddo
  case('smooth')
    print "(A)", 'Solving "'//trim(adjustl(problem))//'" problem'
    t_final = 0.1_R_P
    CFL = 0.7_R_P
    Ns = 1
    Nc = Ns + 2
    Np = Ns + 4
    allocate(initial_state(1:Np, 1:Ni))
    allocate(cp0(1:Ns))
    allocate(cv0(1:Ns))
    variables = 'VARIABLES="x"'
    do s=1, Ns
      variables = variables//' "rho('//trim(str(s,.true.))//')"'
    enddo
    variables = variables//'  "u" "p" "rho" "gamma"'
    BC_L = 'TRA'
    BC_R = 'TRA'
    cp0(1) = 1040._R_P
    cv0(1) = 743._R_P
    do i=1, Ni
      initial_state(:, i) = [1._R_P + 4._R_P / 5._R_P  * sin(         pi * xcenter(i) * 0.5_R_P) + &
                                      1._R_P / 10._R_P * sin(5._R_P * pi * xcenter(i) * 0.5_R_P),  & ! rho(s)
                             0.5_R_P * (xcenter(i) - 0.5_R_P)**4,                                  & ! u
                             10._R_P + 2._R_P * xcenter(i)**4,                                     & ! p
                             1._R_P + 4._R_P / 5._R_P  * sin(         pi * xcenter(i) * 0.5_R_P) + &
                                      1._R_P / 10._R_P * sin(5._R_P * pi * xcenter(i) * 0.5_R_P),  & ! sum(rho(s))
                             cp0/cv0]     ! gamma = cp/cv
    enddo
  case default
    print "(A)", 'Error: unknown problem "'//trim(adjustl(problem))//'"'
    print "(A)", 'Valid problem names are:'
    print "(A)", '  + sod'
    print "(A)", '  + smooth'
    stop
  endselect

  if (trim(adjustl(output_cli))/='unset') then
    output = trim(adjustl(output_cli))//'-'//trim(adjustl(solver))
  else
    output = 'euler_1D_integration-'//trim(adjustl(solver))
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  subroutine average_solution()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Average the solution over an average grid.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P):: i, ii, i1, i2 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if ((av_Ni>0).and.(av_Ni/=Ni)) then
    if (allocated(av_state)) deallocate(av_state) ; allocate(av_state(1:Np, 1:av_Ni))
    do i=1, av_Ni
      i1 = minloc(array=xcenter, dim=1, mask=(xcenter>=av_xnode(i-1)))
      i2 = maxloc(array=xcenter, dim=1, mask=(xcenter<=av_xnode(i)))
      av_state(:, i) = 0._R_P
      do ii=i1, i2
        av_state(:, i) = av_state(:, i) + final_state(:, ii)
      enddo
      av_state(:, i) = av_state(:, i) / (i2-i1+1)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine average_solution

  subroutine save_results(title, basename)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Save results.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN) :: title    !< Plot title.
  character(*), intent(IN) :: basename !< Output basename.
  integer(I_P)             :: rawfile  !< Raw file unit for saving results.
  type(pyplot)             :: plt      !< Plot file handler.
  integer(I_P)             :: v        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (results.or.plots) final_state = domain%output()
  if (results) then
    open(newunit=rawfile, file=basename//'-'//trim(str(Ni,.true.))//'_cells.dat')
    write(rawfile, '(A)')'TITLE="'//title//'"'
    write(rawfile, '(A)')variables
    write(rawfile, '(A)')'ZONE T="FOODIE: '//trim(str(Ni,.true.))//' cells", I='//trim(str(Ni+1,.true.))//&
                         ', J=1, K=1, DATAPACKING=BLOCK, VARLOCATION=([1]=NODAL,[2-'//trim(str(Np+1,.true.)) //']=CELLCENTERED)'
    write(rawfile, '('//trim(str(Ni+1,.true.))//'('//FR_P//',1X))')xnode
    do v=1, Np
      write(rawfile, '('//trim(str(Ni,.true.))//'('//FR_P//',1X))')final_state(v, :)
    enddo
    if ((av_Ni>0).and.(av_Ni/=Ni)) then
      print "(A)", ' Average solution from Ni: '//trim(str(Ni,.true.))//' to av_Ni: '//trim(str(av_Ni,.true.))
      call average_solution
      write(rawfile, '(A)')'ZONE T="FOODIE: '//trim(str(Ni,.true.))//' cells averaged over '//trim(str(av_Ni,.true.))//&
                           ' cells", I='//trim(str(av_Ni+1,.true.))//&
                           ', J=1, K=1, DATAPACKING=BLOCK, VARLOCATION=([1]=NODAL,[2-'//trim(str(Np+1,.true.)) //']=CELLCENTERED)'
      write(rawfile, '('//trim(str(av_Ni+1,.true.))//'('//FR_P//',1X))')av_xnode
      do v=1, Np
        write(rawfile, '('//trim(str(av_Ni,.true.))//'('//FR_P//',1X))')av_state(v, :)
      enddo
    endif
    close(rawfile)
  endif
  if (plots) then
    call plt%initialize(grid=.true., xlabel='x', title=title)
    do v=1, Ns
      call plt%add_plot(x=xcenter, y=final_state(v, :), label='rho('//trim(str(v,.true.))//')', linestyle='b-', linewidth=1)
    enddo
    call plt%add_plot(x=xcenter, y=final_state(Ns+1, :), label='u', linestyle='r-', linewidth=1)
    call plt%add_plot(x=xcenter, y=final_state(Ns+2, :), label='p', linestyle='g-', linewidth=1)
    call plt%add_plot(x=xcenter, y=final_state(Ns+3, :), label='rho', linestyle='o-', linewidth=1)
    call plt%add_plot(x=xcenter, y=final_state(Ns+4, :), label='gamma', linestyle='c-', linewidth=1)
    call plt%savefig(basename//'-'//trim(str(Ni,.true.))//'_cells.png')
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
  integer(I_P)                       :: v        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (time_serie) then
    final_state = domain%output()
    if (present(filename).and.present(title)) then
      open(newunit=tsfile, file=filename)
      write(tsfile, '(A)')'TITLE="'//title//'"'
    endif
    write(tsfile, '(A)')variables//' "t"'
    write(tsfile, '(A)')'ZONE T="'//str(n=t)//'", I='//trim(str(Ni+1,.true.))//&
                        ', J=1, K=1, DATAPACKING=BLOCK, VARLOCATION=([1]=NODAL,[2-'//trim(str(Np+2,.true.)) //']=CELLCENTERED)'
    write(tsfile, '('//trim(str(Ni+1,.true.))//'('//FR_P//',1X))')xnode
    do v=1, Np
      write(tsfile, '('//trim(str(Ni,.true.))//'('//FR_P//',1X))')final_state(v, :)
    enddo
    write(tsfile, '('//trim(str(Ni,.true.))//'('//FR_P//',1X))')(t, v=1,Ni)
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
  type(integrator_adams_bashforth) :: ab_integrator         !< Adams-Bashforth integrator.
  integer, parameter               :: ab_steps=4            !< Adams-Bashforth steps number.
  type(euler_1D)                   :: previous(1:ab_steps)  !< Previous time steps solutions.
  integer                          :: step                  !< Time steps counter.
  real(R_P)                        :: dt                    !< Time step.
  real(R_P)                        :: t(1:ab_steps)         !< Times.
  integer(I_P)                     :: s                     !< AB steps counter.
  integer(I_P)                     :: ss                    !< AB substeps counter.
  integer(I_P)                     :: steps_range(1:2)      !< Steps used.
  character(len=:), allocatable    :: title                 !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of Adams-Bashforth class of solvers'
  steps_range = [1, ab_steps] ; if (stages_steps>0) steps_range = [stages_steps, stages_steps]
  do s=steps_range(1), steps_range(2)
    print "(A)", ' AB-'//trim(str(s,.true.))
    title = '1D Euler equations integration, explicit Adams-Bashforth, t='//str(n=t_final)//trim(str( s,.true.))//' steps'
    call ab_integrator%init(steps=s)
    select case(s)
    case(1)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, steps=s, ord=1)
      call rk_integrator%init(stages=s)
    case(2, 3)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, steps=s, ord=3)
      call rk_integrator%init(stages=s)
    case(4, 5)
      call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, steps=s, ord=7)
      call rk_integrator%init(stages=5)
    endselect
    t = 0._R_P
    call save_time_serie(title=title, filename=output//'-'//trim(str( s,.true.))//'-time_serie.dat', t=t(s))
    step = 1
    do while(t(s)<t_final)
      if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
      dt = domain%dt(Nmax=0, Tmax=t_final, t=t(s), CFL=0.1_R_P*CFL)
      if (s>=step) then
        ! the time steps from 1 to s - 1 must be computed with other scheme...
        call rk_integrator%integrate(U=domain, stage=rk_stage, dt=dt, t=t(s))
        previous(step) = domain
        if (step>1) then
          t(step) = t(step-1) + dt
        else
          t(step) = dt
        endif
      else
        call ab_integrator%integrate(U=domain, previous=previous(1:s), dt=dt, t=t)
        do ss=1, s-1
          t(ss) = t(ss + 1)
        enddo
        t(s) = t(s) + dt
      endif
      step = step + 1
      call save_time_serie(t=t(s))
    enddo
    call save_time_serie(t=t(s), finish=.true.)
    call save_results(title=title, basename=output//'-'//trim(str( s,.true.)))
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
  character(len=:), allocatable   :: title            !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of explicit Euler solver'
  title = '1D Euler equations integration, explicit Euler, t='//str(n=t_final)
  call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0)
  t = 0._R_P
  call save_time_serie(title=title, filename=output//'-time_serie.dat', t=t)
  do while(t<t_final)
    if (verbose) print "(A)", '  Time step: '//str(n=dt)//', Time: '//str(n=t)
    dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=CFL)
    call euler_integrator%integrate(U=domain, dt=dt, t=t)
    t = t + dt
    call save_time_serie(t=t)
  enddo
  call save_time_serie(t=t, finish=.true.)
  call save_results(title=title, basename=output)
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_euler

  subroutine test_leapfrog()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit leapfrog class of ODE solvers.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(tvd_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
  integer, parameter               :: rk_stages=2           !< Runge-Kutta stages number.
  type(euler_1D)                   :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
  type(euler_1D)                   :: filter                !< Filter displacement.
  type(leapfrog_integrator)        :: lf_integrator         !< Leapfrog integrator.
  type(euler_1D)                   :: previous(1:2)         !< Previous time steps solutions.
  integer                          :: step                  !< Time steps counter.
  real(R_P)                        :: dt                    !< Time step.
  real(R_P)                        :: t                     !< Time.
  character(len=:), allocatable    :: title                 !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of leapfrog (RAW filtered) class of solvers'
  title = '1D Euler equations integration, explicit leapfrog (RAW filtered), t='//str(n=t_final)
  call lf_integrator%init(nu=1.0_R_P, alpha=0._R_P)
  call rk_integrator%init(stages=rk_stages)
  call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, steps=2, ord=3)
  t = 0._R_P
  call save_time_serie(title=title, filename=output//'-time_serie.dat', t=t)
  step = 1
  do while(t<t_final)
    if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
    dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=0.1_R_P*CFL)
    if (2>=step) then
      ! the time steps from 1 to s - 1 must be computed with other scheme...
      call rk_integrator%integrate(U=domain, stage=rk_stage, dt=dt, t=t)
      previous(step) = domain
    else
      call lf_integrator%integrate(U=domain, previous=previous, dt=dt, t=t, filter=filter)
    endif
    t = t + dt
    step = step + 1
    call save_time_serie(t=t)
  enddo
  call save_time_serie(t=t, finish=.true.)
  call save_results(title=title, basename=output)
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
  integer(I_P)                    :: stages_range(1:2)     !< Stages used.
  character(len=:), allocatable   :: title                 !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of low storage (2N) Runge-Kutta class of solvers'
  stages_range = [1, rk_stages] ; if (stages_steps>0) stages_range = [stages_steps, stages_steps]
  do s=stages_range(1), stages_range(2)
    if (s==2) cycle ! 2 stages not yet implemented
    if (s==3) cycle ! 3 stages not yet implemented
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(s,.true.))
    title = '1D Euler equations integration, explicit low storage Runge-Kutta, t='//str(n=t_final)//trim(str( s,.true.))//' stages'
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
    call save_time_serie(title=title, filename=output//'-'//trim(str( s,.true.))//'-time_serie.dat', t=t)
    do while(t<t_final)
      if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
      dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=CFL)
      call rk_integrator%integrate(U=domain, stage=rk_stage, dt=dt, t=t)
      t = t + dt
      call save_time_serie(t=t)
    enddo
    call save_time_serie(t=t, finish=.true.)
    call save_results(title=title, basename=output//'-'//trim(str( s,.true.)))
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
  integer(I_P)                     :: stages_range(1:2)     !< Stages used.
  character(len=:), allocatable    :: title                 !< Output files title.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of TVD/SSP Runge-Kutta class of solvers'
  stages_range = [1, rk_stages] ; if (stages_steps>0) stages_range = [stages_steps, stages_steps]
  do s=stages_range(1), stages_range(2)
    if (s==4) cycle ! 4 stages not yet implemented
    print "(A)", ' RK-'//trim(str(s,.true.))
    title = '1D Euler equations integration, explicit TVD/SSP Runge-Kutta, t='//str(n=t_final)//trim(str( s,.true.))//' stages'
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
    call save_time_serie(title=title, filename=output//'-'//trim(str( s,.true.))//'-time_serie.dat', &
                         t=t)
    do while(t<t_final)
      if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
      dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=CFL)
      call rk_integrator%integrate(U=domain, stage=rk_stage(1:s), dt=dt, t=t)
      t = t + dt
      call save_time_serie(t=t)
    enddo
    call save_time_serie(t=t, finish=.true.)
    call save_results(title=title, basename=output//'-'//trim(str( s,.true.)))
  enddo
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_tvd_rk
endprogram integrate_euler_1D
