!< Test without FOODIE with the integration of Euler 1D PDEs system by OpenMP paradigm.
program integrate_euler_1D_openmp_no_foodie
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test without FOODIE with the integration of Euler 1D PDEs system by OpenMP paradigm.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, FR_P, str
use type_euler_1D_openmp_no_foodie, only : euler_1D_omp_nf, tvd_runge_kutta_integrator
use Data_Type_Command_Line_Interface, only : Type_Command_Line_Interface
use pyplot_module, only :  pyplot
#ifdef OPENMP
use OMP_LIB
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Command_Line_Interface) :: cli                   !< Command line interface handler.
type(tvd_runge_kutta_integrator)  :: rk_integrator         !< Runge-Kutta integrator.
integer, parameter                :: rk_stages=5           !< Runge-Kutta stages number.
type(euler_1D_omp_nf)             :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
real(R_P)                         :: dt                    !< Time step.
real(R_P)                         :: t                     !< Time.
type(euler_1D_omp_nf)             :: domain                !< Domain of Euler equations.
real(R_P),    parameter           :: CFL=0.7_R_P           !< CFL value.
integer(I_P), parameter           :: Ns=1                  !< Number of differnt initial gas species.
integer(I_P), parameter           :: Nc=Ns+2               !< Number of conservative variables.
integer(I_P), parameter           :: Np=Ns+4               !< Number of primitive variables.
character(3)                      :: BC_L                  !< Left boundary condition type.
character(3)                      :: BC_R                  !< Right boundary condition type.
integer(I_P)                      :: Ni                    !< Number of grid cells.
real(R_P)                         :: Dx                    !< Space step discretization.
real(R_P)                         :: cp0(1:Ns)             !< Specific heat at constant pressure.
real(R_P)                         :: cv0(1:Ns)             !< Specific heat at constant volume.
real(R_P), allocatable            :: initial_state(:,:)    !< Initial state of primitive variables.
real(R_P), allocatable            :: x(:)                  !< Cell center x-abscissa values.
real(R_P), allocatable            :: final_state(:,:)      !< Final state.
integer(I_P)                      :: error                 !< Error handler.
integer(I_P)                      :: steps_max             !< Maximum number of time steps.
integer(I_P)                      :: omp_threads           !< Number of OpenMP threads.
logical                           :: plots                 !< Flag for activating plots saving.
logical                           :: results               !< Flag for activating results saving.
logical                           :: time_serie            !< Flag for activating time serie-results saving.
logical                           :: verbose               !< Flag for activating more verbose output.
integer(I_P)                      :: profiling(1:2)        !< Tic-toc profiling counters.
integer(I_P)                      :: count_rate            !< Counting rate of system clock.
real(R_P)                         :: system_clocks         !< Profiling result.
integer(I_P)                      :: steps                 !< Time steps counter.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! setting Command Line Interface
call cli%init(progname    = 'euler-1D-openmp-no-foodie',                                          &
              authors     = 'Fortran-FOSS-Programmers',                                           &
              license     = 'GNU GPLv3',                                                          &
              description = 'Test 1D Euler equations integration without FOODIE, OpenMP enabled', &
              examples    = ["euler-1D-openmp-no-foodie --results  ",                             &
                             "euler-1D-openmp-no-foodie -r -t -v -p",                             &
                             "euler-1D-openmp-no-foodie            ",                             &
                             "euler-1D-openmp-no-foodie --plots -r "])
call cli%add(switch='--Ni', help='Number finite volumes used', required=.false., act='store', def='100', error=error)
call cli%add(switch='--steps', help='Number time steps performed', required=.false., act='store', def='30', error=error)
call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--plots', switch_ab='-p', help='Save plots of results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--tserie', switch_ab='-t', help='Save time-serie-results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.', error=error)
call cli%add(switch='--omp_threads', help='Number of OpenMP threads used', required=.false., act='store', def='1', error=error)
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='--omp_threads', val=omp_threads, error=error) ; if (error/=0) stop
call cli%get(switch='--Ni', val=Ni, error=error) ; if (error/=0) stop
call cli%get(switch='--steps', val=steps_max, error=error) ; if (error/=0) stop
call cli%get(switch='-r', val=results, error=error) ; if (error/=0) stop
call cli%get(switch='-p', val=plots, error=error) ; if (error/=0) stop
call cli%get(switch='-t', val=time_serie, error=error) ; if (error/=0) stop
call cli%get(switch='--verbose', val=verbose, error=error) ; if (error/=0) stop
#ifdef OPENMP
call omp_set_dynamic(.false.)
call omp_set_num_threads(omp_threads)
! check OpenMP parallel environment correctness
!$OMP PARALLEL      &
!$OMP DEFAULT(none) &
!$OMP SHARED(omp_threads)
omp_threads = omp_get_num_threads()
!$OMP END PARALLEL
#endif
call init()
call rk_integrator%init(stages=rk_stages)
call domain%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, ord=7)
t = 0._R_P
call save_time_serie(title='FOODIE test: 1D Euler equations integration, explicit TVD Runge-Kutta'// &
                           trim(str(.true., rk_stages))//' stages', &
                     filename='euler_1D_omp_nf_integration-tvdrk-'//trim(str(.true., rk_stages))//'-time_serie.dat', &
                     t=t)
system_clocks = 0._R_P
do steps=1, steps_max
  if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
  dt = domain%dt(Nmax=steps_max, Tmax=0._R_P, t=t, CFL=CFL)
  call system_clock(profiling(1), count_rate)
  call rk_integrator%integrate(U=domain, stage=rk_stage, dt=dt, t=t)
  call system_clock(profiling(2), count_rate)
  system_clocks = system_clocks + real(profiling(2) - profiling(1), kind=R_P)/count_rate
  t = t + dt
  call save_time_serie(t=t)
enddo
system_clocks = system_clocks / steps_max
if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
call save_time_serie(t=t, finish=.true.)
call save_results(title='FOODIE test: 1D Euler equations integration, explicit TVD Runge-Kutta'// &
                        trim(str(.true., rk_stages))//' stages', &
                  filename='euler_1D_omp_nf_integration-tvdrk-'//trim(str(.true., rk_stages)))
print "(I5,A,F23.15)", omp_threads, ' ', system_clocks
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
  allocate(x(1:Ni))
  allocate(initial_state(1:Np, 1:Ni))
  Dx=1._R_P/Ni
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
endprogram integrate_euler_1D_openmp_no_foodie
