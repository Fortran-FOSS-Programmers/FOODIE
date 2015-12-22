!< FOODIE, Fortran Object oriented Ordinary Differential Equations integration library.
module foodie
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE, Fortran Object oriented Ordinary Differential Equations integration library.
!<
!< FOODIE is a KISS library for solving systems of Ordinary Differential Equation (ODE) into the Initial Values Problems (IVP)
!< contest. The mathematical formulation of the problem is:
!<
!< $$U_t = R(t,U)$$
!< $$U_0 = F(0)$$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function and *F* is the (vectorial) initial conditions function.
!<
!< FOODIE is aimed to the time-like integration of the above system of ODE. To this aim, different numerical schemes are provided:
!<
!<+ *explicit Adams-Bashforth* class of schemes:
!<    + 1 step, namely the explicit forward Euler scheme, 1st order accurate;
!<    + 2 steps, 2nd order accurate;
!<    + 3 steps, 3rd order accurate;
!<    + 4 steps, 4th order accurate;
!<+ *forward explicit Euler* scheme, a 1st order accurate;
!<+ *explicit Leapfrog*:
!<    + Unfiltered, 2nd order accurate, (mostly) unstable;
!<    + Robert-Asselin filtered, 1st order accurate;
!<    + Robert-Asselin-Williams filter, 2nd order accurate;
!<+ *explicit low storage Runge-Kutta 2N* class schemes:
!<    + LS(1,1): 1 stage, 1st order accurate, namely the forward explicit Euler one;
!<    + LS(5,4): 5 stages, 4th order accurate;
!<+ *explicit TVD or SSP Runge-Kutta* class schemes:
!<    + TVD(1,1): 1 stage, 1st order accurate, namely the forward explicit Euler one;
!<    + SSP(2,2): 2 stages, 2nd order accurate;
!<    + SSP(3,3): 3 stages, 3rd order accurate;
!<    + SSP(5,4): 5 stages, 4th order accurate;
!<+ *explicit embedded Runge-Kutta* class schemes:
!<    + DP(7,4): 7 stages, 4th order accurate, Dormand and Prince scheme;
!<+ *implicit Adams-Moulton* class of schemes:
!<    + 0 step, namely the implicit backward Euler scheme, 1st order accurate;
!<    + 1 step, 2nd order accurate;
!<    + 2 steps, 3rd order accurate;
!<    + 3 steps, 4th order accurate;
!<+ *predictor-corrector Adams-Bashforth-Moulton* class of schemes:
!<    + P=AB(1)-C=AM(0) step, namely the explicit/implicit forward/backward Euler scheme, 1st order accurate;
!<    + P=AB(2)-C=AM(1) step, 2nd order accurate;
!<    + P=AB(3)-C=AM(2) steps, 3rd order accurate;
!<    + P=AB(4)-C=AM(3) steps, 4th order accurate;
!<
!<### Usage
!<
!< FOODIE schemes must be applied to only subclass extensions of the abstract class *integrand*.
!<
!< To use FOODIE you must:
!<
!<#### extend integrand abstract class provided by FOODIE implementing your concrete integrand field
!<
!< For example for the Lorenz' ODE system
!<
!<```fortran
!< type, extends(integrand) :: lorenz
!<   !< Lorenz equations field.
!<   !<
!<   !< It is a FOODIE integrand class.
!<   private
!<   real(R_P), dimension(:), allocatable :: state        !< Solution vector.
!<   real(R_P)                            :: sigma=0._R_P !< Lorenz \(\sigma\).
!<   real(R_P)                            :: rho=0._R_P   !< Lorenz \(\rho\).
!<   real(R_P)                            :: beta=0._R_P  !< Lorenz \(\beta\).
!<   contains
!<     procedure, pass(self), public :: t => dLorenz_dt                                 !< Time derivate, resiuduals function.
!<     procedure, pass(lhs),  public :: integrand_multiply_real => lorenz_multiply_real !< lorenz * real operator.
!<     procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_lorenz !< Real * Lorenz operator.
!<     procedure, pass(lhs),  public :: add => add_lorenz                               !< Lorenz + Lorenz oprator.
!<     procedure, pass(lhs),  public :: assign_integrand => lorenz_assign_lorenz        !< Lorenz = Lorenz.
!<     procedure, pass(lhs),  public :: assign_real => lorenz_assign_real               !< Lorenz = real.
!<     ...
!< endtype lorenz
!<```
!<
!<#### use one of the provided FOODIE integrator
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
use foodie_adt_integrand, only : integrand
use foodie_integrator_adams_bashforth, only : adams_bashforth_integrator
use foodie_integrator_adams_bashforth_moulton, only : adams_bashforth_moulton_integrator
use foodie_integrator_adams_moulton, only : adams_moulton_integrator
use foodie_integrator_backward_differentiation_formula, only : bdf_integrator
use foodie_integrator_emd_runge_kutta, only : emd_runge_kutta_integrator
use foodie_integrator_euler_explicit, only : euler_explicit_integrator
use foodie_integrator_leapfrog, only : leapfrog_integrator
use foodie_integrator_low_storage_runge_kutta, only : ls_runge_kutta_integrator
use foodie_integrator_tvd_runge_kutta, only : tvd_runge_kutta_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: integrand
public :: adams_bashforth_integrator
public :: adams_bashforth_moulton_integrator
public :: adams_moulton_integrator
public :: bdf_integrator
public :: emd_runge_kutta_integrator
public :: euler_explicit_integrator
public :: leapfrog_integrator
public :: ls_runge_kutta_integrator
public :: tvd_runge_kutta_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule foodie
