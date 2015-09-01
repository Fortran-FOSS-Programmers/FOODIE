!< Test FOODiE with the integration of Euler 1D PDEs system.
program integrate_euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test FOODiE with the integration of Euler 1D PDEs system.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! use IR_Precision, only : R_P, I_P, str
! use type_euler_1D, only : euler_1D
! use foodie, only : euler_explicit_integrator, tvd_runge_kutta_integrator, adams_bashforth_integrator
! use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
! type(euler_explicit_integrator)  :: euler_integrator            !< Euler integrator.
! type(tvd_runge_kutta_integrator) :: rk_integrator               !< Runge-Kutta integrator.
! type(adams_bashforth_integrator) :: ab_integrator               !< Adams-Bashforth integrator.
! integer, parameter               :: rk_stages=5                 !< Runge-Kutta stages number.
! integer, parameter               :: ab_steps=3                                      !< Adams-Bashforth steps number.
! type(burgers)                    :: rk_stage(1:rk_stages)       !< Runge-Kutta stages.
! real(R_P), parameter             :: pi = 4._R_P * atan(1._R_P)  !< Pi greek.
! type(burgers)                    :: domain                      !< Burgers field domain.
! integer                          :: step                        !< Time steps counter.
! real(R_P)                        :: dt                          !< Time step.
! real(R_P)                        :: t=0.                        !< Initial time.
! real(R_P)                        :: t_final=0.6                 !< Final time.
! integer,   parameter             :: Ni=100                      !< Number of grid nodes.
! real(R_P), parameter             :: h=2._R_P * pi / Ni          !< Space step discretization.
! real(R_P), parameter             :: nu=1._R_P                   !< Viscosity.
! real(R_P)                        :: initial_state(1:Ni)         !< Initial state.
! real(R_P)                        :: x(1:Ni)                     !< Nodes values.
! real(R_P), allocatable           :: final_state(:)              !< Final state.
! type(pyplot)                     :: plt                         !< Plot file handler.
! integer(I_P)                     :: s                           !< RK stages/AB steps counter.
! integer(I_P)                     :: i                           !< Space counter.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram integrate_euler_1D
