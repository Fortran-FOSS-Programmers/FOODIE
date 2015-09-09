module foodie
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODiE, Fortran Object oriented Ordinary Differential Equations integration library.
!<
!< FOODiE is a KISS library for solving systems of Ordinary Differential Equation (ODE) into the Initial Values Problems (IVP)
!< contest. The mathematical formulation of the problem is:
!<
!< $$U_t = R(t,U)$$
!< $$U_0 = F(0)$$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function and *F* is the (vectorial) initial conditions function.
!<
!< FOODiE is aimed to the time-like integration of the above system of ODE. To this aim, different numerical schemes are provided:
!<
!<+ *forward explicit Euler* scheme, a 1st order accurate scheme;
!<+ *TVD or SSP Runge-Kutta* class schemes:
!<    + TVD(1,1): 1 stage, 1st order scheme, namely the forward explicit Euler one;
!<    + SSP(2,2): 2 stages, 2nd order scheme;
!<    + SSP(3,3): 3 stages, 3rd order scheme;
!<    + SSP(5,4): 5 stages, 4th order scheme;
!<
!<### Usage
!<
!< FOODiE schemes must be applied to only subclass extensions of the abstract class *type_integrand*.
!<
!< To use FOODiE you must:
!<
!<#### extend type_integrand abstract class provided by FOODiE implementing your concrete integrand field
!<
!< For example for the Lorenz' ODE system
!<
!<```fortran
!< type, extends(integrand) :: lorenz
!<   !< Lorenz equations field.
!<   !<
!<   !< It is a FOODiE integrand class.
!<   private
!<   real(R_P), dimension(:), allocatable :: state        !< Solution vector.
!<   real(R_P)                            :: sigma=0._R_P !< Lorenz \(\sigma\).
!<   real(R_P)                            :: rho=0._R_P   !< Lorenz \(\rho\).
!<   real(R_P)                            :: beta=0._R_P  !< Lorenz \(\beta\).
!<   contains
!<     procedure, pass(self), public :: output                                          !< Extract Lorenz field.
!<     procedure, pass(self), public :: t => dLorenz_dt                                 !< Time derivate, resiuduals function.
!<     procedure, pass(lhs),  public :: integrand_multiply_real => lorenz_multiply_real !< lorenz * real operator.
!<     procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_lorenz !< Real * Lorenz operator.
!<     procedure, pass(lhs),  public :: add => add_lorenz                               !< Lorenz + Lorenz oprator.
!<     procedure, pass(lhs),  public :: assign_integrand => lorenz_assign_lorenz        !< Lorenz = Lorenz.
!<     procedure, pass(lhs),  public :: assign_real => lorenz_assign_real               !< Lorenz = real.
!< endtype lorenz
!<```
!<
!<#### use one of the provided FOODiE integrator
!<
!< For example using the forward explicit Euler scheme to the above Lorenz' ODE system
!<
!<```fortran
!< use foodie, only : euler_explicit_integrator
!< use type_lorenz, only : lorenz
!< type(euler_explicit_integrator) :: euler_integrator
!< type(lorenz)                    :: attractor
!< real                            :: dt=0.01
!< do step = 1, num_steps
!<   call euler_integrator%integrate(field=attractor, dt=dt)
!< enddo
!<```
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use type_integrand, only : integrand
use foodie_integrator_euler_explicit, only : euler_explicit_integrator
use foodie_integrator_tvd_runge_kutta, only : tvd_runge_kutta_integrator
use foodie_integrator_adams_bashforth, only : adams_bashforth_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: integrand
public :: euler_explicit_integrator
public :: tvd_runge_kutta_integrator
public :: adams_bashforth_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule foodie
