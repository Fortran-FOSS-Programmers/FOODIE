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
use foodie_adt_integrand, only : integrand
use foodie_kinds, only : I_P, R_P
use foodie_utils, only : is_admissible
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: euler_explicit_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
character(len=1), parameter :: supported_stages_steps='1' !< List of supported stages/steps number. Valid format is `1-2,4,9-23...`.
integer(I_P),     parameter :: min_ss=1                   !< Minimum number of stages/steps supported.
integer(I_P),     parameter :: max_ss=1                   !< Maximum number of stages/steps supported.

type :: euler_explicit_integrator
  !< FOODIE integrator: provide an explicit Euler scheme, it being 1st order accurate.
  !<
  !< @note The integrator can be used directly without any initialization.
  private
  contains
    private
    procedure, nopass, public :: integrate        !< Integrate integrand field.
    procedure, nopass, public :: min_stages_steps !< Return the minimum number of stages/steps supported.
    procedure, nopass, public :: max_stages_steps !< Return the maximum number of stages/steps supported.
    procedure, nopass, public :: is_supported     !< Check if the queried number of stages/steps is supported or not.
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

  pure function min_stages_steps()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the minimum number of stages/steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P) :: min_stages_steps !< Minimum number of stages/steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  min_stages_steps = min_ss
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction min_stages_steps

  pure function max_stages_steps()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the maximum number of stages/steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P) :: max_stages_steps !< Maximum number of stages/steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  max_stages_steps = max_ss
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction max_stages_steps

  elemental function is_supported(stages_steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Check if the queried number of stages/steps is supported or not.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN) :: stages_steps !< Number of stages/steps used.
  logical                  :: is_supported !< Is true is the stages number is in *supported_stages_steps*.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_supported = is_admissible(n=stages_steps, adm_range=trim(supported_stages_steps))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_supported

endmodule foodie_integrator_euler_explicit
