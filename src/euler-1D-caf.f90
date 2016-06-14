!< Test FOODIE with the integration of Euler 1D PDEs system.
program integrate_euler_1D_caf
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODIE with the integration of Euler 1D PDEs system.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use flap, only : command_line_interface
use foodie, only : tvd_runge_kutta_integrator
use IR_Precision, only : R_P, I_P, FR_P, str, strz
use pyplot_module, only : pyplot
use type_euler_1D_caf, only : euler_1D_caf
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(command_line_interface)     :: cli                   !< Command line interface handler.
type(tvd_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
integer, parameter               :: rk_stages=5           !< Runge-Kutta stages number.
type(euler_1D_caf)               :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
integer, parameter               :: ord=7                 !< Space reconstruciton order,
real(R_P)                        :: t                     !< Time.
real(R_P),    parameter          :: CFL=0.7_R_P           !< CFL value.
integer(I_P), parameter          :: Ns=1                  !< Number of differnt initial gas species.
integer(I_P), parameter          :: Nc=Ns+2               !< Number of conservative variables.
integer(I_P), parameter          :: Np=Ns+4               !< Number of primitive variables.
character(len=:), allocatable    :: BC_L                  !< Left boundary condition type.
character(len=:), allocatable    :: BC_R                  !< Right boundary condition type.
real(R_P)                        :: Dx                    !< Space step discretization.
real(R_P)                        :: cp0(1:Ns)             !< Specific heat at constant pressure.
real(R_P)                        :: cv0(1:Ns)             !< Specific heat at constant volume.
real(R_P), allocatable           :: initial_state(:,:)    !< Initial state of primitive variables.
real(R_P), allocatable           :: x(:)                  !< Cell center x-abscissa values.
real(R_P), allocatable           :: final_state(:,:)      !< Final state.
integer(I_P)                     :: error                 !< Error handler.
integer(I_P)                     :: profiling(1:2)        !< Tic-toc profiling counters.
integer(I_P)                     :: count_rate            !< Counting rate of system clock.
real(R_P)                        :: system_clocks         !< Profiling result.
integer(I_P)                     :: steps                 !< Time steps counter.
type(euler_1D_caf)               :: domain                !< Domain of Euler equations.
! coarrays-related variables
integer(I_P)                     :: Ni_image              !< Space dimension of local image.
#ifdef CAF
integer(I_P)                     :: Ni[*]                 !< Number of grid cells.
integer(I_P)                     :: steps_max[*]          !< Maximum number of time steps.
logical                          :: results[*]            !< Flag for activating results saving.
logical                          :: plots[*]              !< Flag for activating plots saving.
logical                          :: time_serie[*]         !< Flag for activating time serie-results saving.
logical                          :: verbose[*]            !< Flag for activating more verbose output.
real(R_P), allocatable           :: Dt(:)[:]              !< Time step.
#else
integer(I_P)                     :: Ni                    !< Number of grid cells.
integer(I_P)                     :: steps_max             !< Maximum number of time steps.
logical                          :: results               !< Flag for activating results saving.
logical                          :: plots                 !< Flag for activating plots saving.
logical                          :: time_serie            !< Flag for activating time serie-results saving.
logical                          :: verbose               !< Flag for activating more verbose output.
real(R_P), allocatable           :: Dt(:)                 !< Time step.
#endif
integer(I_P)                     :: me                    !< ID of this_image()
integer(I_P)                     :: we                    !< Number of CAF images used.
character(len=:), allocatable    :: id                    !< My ID.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef CAF
me = this_image()
we = num_images()
#else
me = 1
we = 1
#endif
id = trim(strz(3, me))//'> '
! setting Command Line Interface
call cli%init(progname    = 'euler-1D-caf',                                                       &
              authors     = 'Fortran-FOSS-Programmers',                                           &
              license     = 'GNU GPLv3',                                                          &
              description = 'Test FOODIE library on 1D Euler equations integration, CAF enabled', &
              examples    = ["euler-1D-caf --results  ",                                          &
                             "euler-1D-caf -r -t -v -p",                                          &
                             "euler-1D-caf            ",                                          &
                             "euler-1D-caf --plots -r "])
call cli%add(switch='--Ni', help='Number finite volumes used', required=.false., act='store', def='100', error=error)
call cli%add(switch='--steps', help='Number time steps performed', required=.false., act='store', def='30', error=error)
call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--plots', switch_ab='-p', help='Save plots of results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--tserie', switch_ab='-t', help='Save time-serie-results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.', error=error)
! parsing Command Line Interface
if (me==1) then ! only master image do CLI stuffs
  call cli%parse(error=error)
  call cli%get(switch='--Ni', val=Ni, error=error) ; if (error/=0) stop
  call cli%get(switch='--steps', val=steps_max, error=error) ; if (error/=0) stop
  call cli%get(switch='-r', val=results, error=error) ; if (error/=0) stop
  call cli%get(switch='-p', val=plots, error=error) ; if (error/=0) stop
  call cli%get(switch='-t', val=time_serie, error=error) ; if (error/=0) stop
  call cli%get(switch='--verbose', val=verbose, error=error) ; if (error/=0) stop
endif
call init()
system_clocks = 0._R_P
do steps=1, steps_max
  if (verbose) print "(A)", id//'    Time step: '//str(n=Dt(me))//', Time: '//str(n=t)
  Dt(me) = domain%dt(Nmax=steps_max, Tmax=0._R_P, t=t, CFL=CFL)
  call synchronize
  call system_clock(profiling(1), count_rate)
  call rk_integrator%integrate(U=domain, stage=rk_stage, Dt=Dt(me), t=t)
  call system_clock(profiling(2), count_rate)
  system_clocks = system_clocks + real(profiling(2) - profiling(1), kind=R_P)/count_rate
  t = t + Dt(me)
  call save_time_serie(t=t)
enddo
system_clocks = system_clocks / steps_max
if (verbose) print "(A)", id//'    Time step: '//str(n=Dt(me))//', Time: '//str(n=t)
call save_time_serie(t=t, finish=.true.)
call save_results(title='FOODIE test: 1D Euler equations integration, explicit TVD Runge-Kutta'// &
                        trim(str(.true., rk_stages))//' stages', &
                  filename='euler_1D_caf_integration-tvdrk-'//trim(str(.true., rk_stages))//'-image-'//trim(strz(3, me)))

print "(A,I5,A,F23.15)", id, we, ' ', system_clocks
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine init()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Initialize the simulation.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P)                  :: i        !< Space counter.
  real(R_P)                     :: x_L      !< Left abscissa of local image.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef CAF
  ! CAF images comunications
  sync all
  if (me/=1) then
    Ni = Ni[1]
    steps_max = steps_max[1]
    results = results[1]
    plots = plots[1]
    time_serie = time_serie[1]
    verbose = verbose[1]
  endif
#endif
  ! init simulation
  if (mod(Ni, we)/=0) error stop 'error: the number of cells Ni must be a multiple of the number of CAF images used!'
  Ni_image = Ni / we
  allocate(x(1:Ni_image))
  allocate(initial_state(1:Np, 1:Ni_image))
  Dx = 1._R_P / Ni
  ! Sod's problem
  cp0(1) = 1040._R_P
  cv0(1) = 743._R_P
  if (we==1) then
      BC_L = 'TRA'
      BC_R = 'TRA'
  else
    if (me==1) then
      BC_L = 'TRA'
      BC_R = 'CON-'//trim(strz(2, me+1))
    elseif (me==we) then
      BC_L = 'CON-'//trim(strz(2, me-1))
      BC_R = 'TRA'
    else
      BC_L = 'CON-'//trim(strz(2, me-1))
      BC_R = 'CON-'//trim(strz(2, me+1))
    endif
  endif
  if (me>1) then
    x_L = Ni_image * Dx * (me - 1)
  else
    x_L = 0._R_P
  endif
  do i=1, Ni_image
    x(i) = x_L + Dx * i - 0.5_R_P * Dx
    if (x(i)<=0.5_R_P) then
      initial_state(:, i) = [1._R_P, & ! rho(s)
                             0._R_P, & ! u
                             1._R_P, & ! p
                             1._R_P, & ! sum(rho(s))
                             cp0/cv0]  ! gamma = cp/cv
    else
      initial_state(:, i) = [0.125_R_P, & ! rho(s)
                             0._R_P,    & ! u
                             0.1_R_P,   & ! p
                             0.125_R_P, & ! sum(rho(s))
                             cp0/cv0]     ! gamma = cp/cv
    endif
  enddo
  if (verbose) then
    print '(A)', id//'image '//trim(str(.true., me))//' of '//trim(str(.true., we))
    print '(A)', id//'Number of total cells: '//trim(str(.true., Ni))
    print '(A)', id//'Number of time steps: '//trim(str(.true., steps_max))
    print '(A)', id//'Save final results: '//trim(str(results))
    print '(A)', id//'Save plots of results: '//trim(str(plots))
    print '(A)', id//'Save time serie of results: '//trim(str(time_serie))
    print '(A)', id//'Left BC: '//BC_L
    print '(A)', id//'Right BC: '//BC_R
    print '(A)', id//'Space resolution: '//trim(str(.true., Dx))
    print '(A)', id//'X(1) X(N): '//trim(str(.true., x(1)))//' '//trim(str(.true., x(Ni_image)))
    print '(A)', id//'Density value: '//trim(str(n=initial_state(1, 1)))//' '//trim(str(n=initial_state(1, Ni_image)))
  endif
  ! initialize integrator and domain
  call rk_integrator%init(stages=rk_stages)
  call domain%init(Ni=Ni_image, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, &
                   me=me, we=we, ord=ord)
#ifdef CAF
  allocate(Dt(1:we)[*])
#else
  allocate(Dt(1:we))
#endif
  ! initialize time serie file
  call save_time_serie(title='FOODIE test: 1D Euler equations integration, explicit TVD Runge-Kutta'// &
                             trim(str(.true., rk_stages))//' stages', &
                       filename='euler_1D_caf_integration-tvdrk-'//&
                                trim(str(.true., rk_stages))//'-image-'//&
                                trim(strz(3, me))//&
                                '-time_serie.dat', &
                       t=t)
  ! initialize time variables
  t = 0._R_P
  Dt = 0._R_P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  subroutine synchronize()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Synchronize CAF images.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P) :: i !< Images counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef CAF
  if (we>1) then
    sync all
    ! reduction on minumum value of Dt
    do i=1, we
      if (i/=me) Dt(i) = Dt(i)[i]
      Dt(me) = min(Dt(me), Dt(i))
    enddo
  endif
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine synchronize

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
    write(rawfile, '(A)')'VARIABLES="x" "rho(1)" "u" "p" "rho" "gamma"'
    write(rawfile, '(A)')'ZONE T="'//str(n=t)//'"'
    do i=1, Ni_image
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
    write(tsfile, '(A)')'VARIABLES="x" "rho(1)" "u" "p" "rho" "gamma"'
    write(tsfile, '(A)')'ZONE T="'//str(n=t)//'"'
    do i=1, Ni_image
      write(tsfile, '('//trim(str(.true.,Np+1))//'('//FR_P//',1X))')x(i), (final_state(v, i), v=1, Np)
    enddo
    if (present(finish)) then
      if (finish) close(tsfile)
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_time_serie
endprogram integrate_euler_1D_caf
