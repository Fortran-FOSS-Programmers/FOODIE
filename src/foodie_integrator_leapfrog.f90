!< FOODIE integrator: provide an explicit class of leapfrog multi-step schemes, 2nd order accurate.
module foodie_integrator_leapfrog
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE integrator: provide an explicit class of leapfrog multi-step schemes, 2nd order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the leapfrog class scheme implemented (see [3]) is:
!<
!< $$ U^{n+2} = U^{n} + 2\Delta t \cdot R(t^{n+1}, U^{n+1}) $$
!<
!< Optionally, the Robert-Asselin-Williams (RAW) filter (see [3]) is applied to the computed integration steps:
!< $$ \Delta = \frac{\nu}{2}(U^{n} - 2 U^{n+1} + U^{n+2}) $$
!< $$ U^{n+1} = U^{n+1} + \Delta * \alpha $$
!< $$ U^{n+2} = U^{n+2} + \Delta * (\alpha-1) $$
!< Note that for \(\alpha=1\) the filter reverts back to the standard Robert-Asselin scheme.
!< The filter coefficients should be taken as \(\nu \in (0,1]\) and \(\alpha \in (0.5,1]\). The default values are
!<
!<  + \(\nu=0.01\)
!<  + \(\alpha=0.53\)
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit. The filter coefficients \(\nu,\,\alpha \) define the actual scheme.
!<
!<#### Bibliography
!<
!< [1] *The integration of a low order spectral form of the primitive meteorological equations*, Robert, A. J., J. Meteor. Soc.
!< Japan,vol. 44, pages 237--245, 1966.
!<
!< [2] *Frequency filter for time integrations*, Asselin, R., Monthly Weather Review, vol. 100, pages 487--490, 1972.
!<
!< [3] *The RAW filter: An improvement to the Robert–Asselin filter in semi-implicit integrations*, Williams, P.D., Monthly
!< Weather Review, vol. 139(6), pages 1996--2007, June 2011.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie_kinds, only : R_P, I_P
use foodie_adt_integrand, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: leapfrog_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: leapfrog_integrator
  !< FOODIE integrator: provide an explicit class of leapfrog multi-step schemes, 2nd order accurate.
  !<
  !< @note The integrator could be used without initialialization (initialize the time filter coefficients) if the defulat values
  !< are suitable for the problem.
  private
  real(R_P) :: nu=0.01_R_P    !< Robert-Asselin filter coefficient.
  real(R_P) :: alpha=0.53_R_P !< Robert-Asselin-Williams filter coefficient.
  contains
    procedure, pass(self), public :: init      !< Initialize (create) the integrator.
    procedure, pass(self), public :: integrate !< Integrate integrand field.
endtype leapfrog_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  elemental subroutine init(self, nu, alpha)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual leapfrog integrator: initialize the filter coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(leapfrog_integrator), intent(INOUT) :: self  !< LF integrator.
  real(R_P), optional,        intent(IN)    :: nu    !< Williams-Robert-Asselin filter coefficient.
  real(R_P), optional,        intent(IN)    :: alpha !< Robert-Asselin filter coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%nu = 0.01_R_P
  self%alpha = 0.53_R_P
  if (present(nu)) self%nu = nu
  if (present(alpha)) self%alpha = alpha
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  subroutine integrate(self, U, previous, Dt, t, filter)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with leapfrog class scheme.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(leapfrog_integrator), intent(IN)    :: self          !< LF integrator.
  class(integrand),           intent(INOUT) :: U             !< Field to be integrated.
  class(integrand),           intent(INOUT) :: previous(1:2) !< Previous time steps solutions of integrand field.
  real(R_P),                  intent(in)    :: Dt            !< Time step.
  real(R_P),                  intent(IN)    :: t             !< Time.
  class(integrand), optional, intent(INOUT) :: filter        !< Filter field displacement.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  U = previous(1) + previous(2)%t(t=t) * (Dt * 2._R_P)
  if (present(filter)) then
    filter = (previous(1) - previous(2) * 2._R_P + U) * self%nu * 0.5_R_P
    previous(2) = previous(2) + filter * self%alpha
    U = U + filter * (self%alpha - 1._R_P)
  endif
  previous(1) = previous(2)
  previous(2) = U
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate
endmodule foodie_integrator_leapfrog
