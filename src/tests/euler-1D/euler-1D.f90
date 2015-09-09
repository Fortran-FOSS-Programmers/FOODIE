!< Test FOODiE with the integration of Euler 1D PDEs system.
program integrate_euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Euler 1D PDEs system.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, FR_P, str
use type_euler_1D, only : euler_1D
use Data_Type_Command_Line_Interface, only : Type_Command_Line_Interface
use foodie, only : euler_explicit_integrator, tvd_runge_kutta_integrator, adams_bashforth_integrator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Command_Line_Interface) :: cli                      !< Command line interface handler.
type(euler_1D)                    :: domain                   !< Domain of Euler equations.
real(R_P),    parameter           :: CFL=0.1_R_P              !< CFL value.
real(R_P),    parameter           :: t_final=0.15_R_P         !< Final time.
integer(I_P), parameter           :: Ni=100                   !< Number of grid cells.
integer(I_P), parameter           :: Ng=3                     !< Number of ghost cells.
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
logical                           :: plots               !< Flag for activating plots saving.
logical                           :: results             !< Flag for activating results saving.
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
! parsing Command Line Interface
call cli%parse(error=error)
call cli%get(switch='-s', val=solver, error=error) ; if (error/=0) stop
call cli%get(switch='-r', val=results, error=error) ; if (error/=0) stop
call cli%get(switch='-p', val=plots, error=error) ; if (error/=0) stop
! create Euler field initial state
call init()
! integrate Euler equation
select case(trim(adjustl(solver)))
case('euler')
  call test_euler()
case('runge-kutta')
  ! call test_rk()
case('adams-bashforth')
  ! call test_ab()
case('all')
  call test_euler()
  ! call test_rk()
  ! call test_ab()
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

  subroutine test_euler()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Test explicit forward Euler ODE solver.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(euler_explicit_integrator) :: euler_integrator !< Euler integrator.
  real(R_P)                       :: dt               !< Time step.
  real(R_P)                       :: t                !< Time.
  integer(I_P)                    :: s                !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  print "(A)", 'Integrating 1D Euler equations by means of explicit Euler solver'
  ! initialize field
  call domain%init(Ni=Ni, Ng=Ng, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0)
  ! integrate field
  t = 0._R_P
  s = 0
  do while(t<t_final)
    print*, s
    dt = domain%dt(Nmax=0, Tmax=t_final, t=t, CFL=CFL)
    call euler_integrator%integrate(field=domain, dt=dt)
    final_state = domain%output()
    call save_results(title='FOODiE test: 1D Euler equations integration at t=0.6, explicit Euler', &
                      filename='euler_integration-euler')
    t = t + dt
    s = s + 1
  enddo
  final_state = domain%output()
  call save_results(title='FOODiE test: 1D Euler equations integration at t=0.6, explicit Euler', &
                    filename='euler_integration-euler')
  print "(A)", 'Finish!'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_euler
endprogram integrate_euler_1D
