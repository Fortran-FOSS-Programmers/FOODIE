!< FOODIE integrator: provide an explicit Euler scheme, it being 1st order accurate.
module foodie_integrator_euler_explicit
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE integrator: provide an explicit Euler scheme, it being 1st order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the Forward Explicit Euler scheme implemented is:
!<
!< $$ U^{n+1} = U^n +\Delta t R(t, U^n) $$
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie_kinds, only : R_P
use foodie_adt_integrand, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: euler_explicit_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: euler_explicit_integrator
  !< FOODIE integrator: provide an explicit Euler scheme, it being 1st order accurate.
  !<
  !< @note The integrator can be used directly without any initialization.
  contains
    procedure, nopass, public :: integrate !< Integrate integrand field.
endtype euler_explicit_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine integrate(U, Dt, t)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with explicit Euler scheme, 1st order.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(integrand),    intent(INOUT) :: U  !< Field to be integrated.
  real(R_P),           intent(IN)    :: Dt !< Time step.
  real(R_P), optional, intent(IN)    :: t  !< Time.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  U = U + U%t(t=t) * Dt
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate
endmodule foodie_integrator_euler_explicit
