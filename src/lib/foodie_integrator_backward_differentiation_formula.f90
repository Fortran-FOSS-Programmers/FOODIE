!< FOODIE integrator: provide an implicit class of Backward Differentiation Formula schemes, from 1st to 6th order accurate.
module foodie_integrator_backward_differentiation_formula
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE integrator: provide an implicit class of Backward Differentiation Formula schemes, from 1st to 6th order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the Backward Differentiation Formula class scheme implemented is:
!<
!< $$ U^{n+N_s} + \sum_{s=1}^{N_s} \alpha_s U^{n+N_s-s} = \Delta t \left[ \beta R(t^{n+N_s}, U^{n+N_s}) \right] $$
!<
!< where \(N_s\) is the number of previous steps considered.
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are implicit. The coefficients \(\alpha_s\) and \(\beta\) define the actual scheme, that is selected accordingly
!< to the number of *steps* used.
!<
!< Currently, the following schemes are available:
!<
!< | ` Step ` | ` beta `   | ` alpha 1 `  | ` alpha 2 ` | ` alpha 3 `  | ` alpha 4 ` | ` alpha 5 ` | ` alpha 6 ` |
!< |----------|------------|--------------|-------------|--------------|-------------|-------------|-------------|
!< | ` 1 `    | `  1 `     | ` -1 `       |             |              |             |             |             |
!< | ` 2 `    | ` 2/3 `    | ` -4/3 `     | ` 1/3 `     |              |             |             |             |
!< | ` 3 `    | ` 6/11 `   | ` -18/11 `   | ` 9/11  `   | ` -2/11 `    |             |             |             |
!< | ` 4 `    | ` 12/25 `  | ` -48/25 `   | ` 36/25 `   | ` -16/25 `   | ` 3/25 `    |             |             |
!< | ` 5 `    | ` 60/137 ` | ` -300/137 ` | ` 300/137 ` | ` -200/137 ` | ` 75/137 `  | ` -12/137 ` |             |
!< | ` 6 `    | ` 60/147 ` | ` -360/147 ` | ` 450/147 ` | ` -400/147 ` | ` 225/147 ` | ` -72/147 ` | ` 10/147 `  |
!<
!<#### Bibliography
!<
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie_adt_integrand, only : integrand
use foodie_kinds, only : I_P, R_P
use foodie_utils, only : is_admissible
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: back_df_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
character(len=99), parameter :: supported_steps='1-6' !< List of supported steps number. Valid format is `1-2,4,9-23...`.
integer(I_P),      parameter :: min_ss=1              !< Minimum number of steps supported.
integer(I_P),      parameter :: max_ss=6              !< Maximum number of steps supported.

type :: back_df_integrator
  !< FOODIE integrator: provide an implicit class of Backward-Differentiation-Formula multi-step schemes, from 1st to 6th order
  !< accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *alpha* and *beta* coefficients) before used.
  !<
  !< ### List of errors status
  !<+ error=0 => no error;
  !<+ error=1 => bad (unsupported) number of required time steps;
  private
  integer(I_P)           :: steps=0   !< Number of time steps.
  real(R_P), allocatable :: a(:)      !< \(\alpha\) coefficients.
  real(R_P)              :: b=0.0_R_P !< \(\beta\) coefficient.
  integer(I_P)           :: error=0   !< Error status flag: trap occurrences of errors.
  contains
    private
    procedure, pass(self), public :: destroy         !< Destroy the integrator.
    procedure, pass(self), public :: init            !< Initialize (create) the integrator.
    procedure, pass(self), public :: integrate       !< Integrate integrand field.
    procedure, pass(self), public :: update_previous !< Cyclic update previous time steps.
    procedure, nopass,     public :: min_steps       !< Return the minimum number of steps supported.
    procedure, nopass,     public :: max_steps       !< Return the maximum number of steps supported.
    procedure, nopass,     public :: is_supported    !< Check if the queried number of steps is supported or not.
    final                         :: finalize        !< Finalize object.
endtype back_df_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  elemental subroutine init(self, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual BDF integrator: initialize the *alpha* and *beta* coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(back_df_integrator), intent(INOUT) :: self  !< BDF integrator.
  integer(I_P),              intent(IN)    :: steps !< Number of time steps used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = steps
  if (allocated(self%a)) deallocate(self%a) ; allocate(self%a(1:steps)) ; self%a = 0.0_R_P
  select case(steps)
  case(1)
    self%a(1) = -1.0_R_P
    self%b = 1.0_R_P
  case(2)
    self%a(1) = 1.0_R_P/3.0_R_P
    self%a(2) = -4.0_R_P/3.0_R_P
    self%b = 2.0_R_P/3.0_R_P
  case(3)
    self%a(1) = -2.0_R_P/11.0_R_P
    self%a(2) = 9.0_R_P/11.0_R_P
    self%a(3) = -18.0_R_P/11.0_R_P
    self%b = 6.0_R_P/11.0_R_P
  case(4)
    self%a(1) = 3.0_R_P/25.0_R_P
    self%a(2) = -16.0_R_P/25.0_R_P
    self%a(3) = 36.0_R_P/25.0_R_P
    self%a(4) = -48.0_R_P/25.0_R_P
    self%b = 12.0_R_P/25.0_R_P
  case(5)
    self%a(1) = -12.0_R_P/137.0_R_P
    self%a(2) = 75.0_R_P/137.0_R_P
    self%a(3) = -200.0_R_P/137.0_R_P
    self%a(4) = 300.0_R_P/137.0_R_P
    self%a(5) = -300.0_R_P/137.0_R_P
    self%b = 60.0_R_P/137.0_R_P
  case(6)
    self%a(1) = 10.0_R_P/147.0_R_P
    self%a(2) = -72.0_R_P/147.0_R_P
    self%a(3) = 225.0_R_P/147.0_R_P
    self%a(4) = -400.0_R_P/147.0_R_P
    self%a(5) = 450.0_R_P/147.0_R_P
    self%a(6) = -360.0_R_P/147.0_R_P
    self%b = 60.0_R_P/147.0_R_P
  case default
    ! bad (unsupported) number of required time steps
    self%error = 1
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy the integrator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(back_df_integrator), intent(INOUT) :: self !< BDF integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = 1
  if (allocated(self%a)) deallocate(self%a)
  self%b = 0.0_R_P
  self%error = 0
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine integrate(self, U, previous, Dt, t, iterations, autoupdate)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with BDF class scheme.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(back_df_integrator),  intent(IN)    :: self         !< Actual BDF integrator.
  class(integrand),           intent(INOUT) :: U            !< Field to be integrated.
  class(integrand),           intent(INOUT) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                  intent(IN)    :: Dt           !< Time steps.
  real(R_P),                  intent(IN)    :: t(:)         !< Times.
  integer(I_P), optional,     intent(IN)    :: iterations   !< Fixed point iterations.
  logical,      optional,     intent(IN)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  integer(I_P)                              :: iterations_  !< Fixed point iterations.
  logical                                   :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  class(integrand), allocatable             :: delta        !< Delta RHS for fixed point iterations.
  integer(I_P)                              :: s            !< Steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  iterations_ = 1 ; if (present(iterations)) iterations_ = iterations
  allocate(delta, source=previous(self%steps) * (-self%a(self%steps)))
  do s=1, self%steps - 1
    delta = delta + previous(s) * (-self%a(s))
  enddo
  do s=1, iterations
    U = delta + U%t(t=t(self%steps) + Dt) * (Dt * self%b)
  enddo
  if (autoupdate_) call self%update_previous(U=U, previous=previous)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate

  pure function min_steps()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the minimum number of steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P) :: min_steps !< Minimum number of steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  min_steps = min_ss
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction min_steps

  pure function max_steps()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the maximum number of steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P) :: max_steps !< Maximum number of steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  max_steps = max_ss
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction max_steps

  elemental function is_supported(steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Check if the queried number of steps is supported or not.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN) :: steps        !< Number of time steps used.
  logical                  :: is_supported !< Is true is the steps number is in *supported_steps*.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_supported = is_admissible(n=steps, adm_range=trim(supported_steps))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_supported

  subroutine update_previous(self, U, previous)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cyclic update previous time steps.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(back_df_integrator), intent(IN)    :: self         !< Actual BDF integrator.
  class(integrand),          intent(IN)    :: U            !< Field to be integrated.
  class(integrand),          intent(INOUT) :: previous(1:) !< Previous time steps solutions of integrand field.
  integer(I_P)                             :: s            !< Steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do s=1, self%steps - 1
    previous(s) = previous(s + 1)
  enddo
  previous(self%steps) = U
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_previous

  ! private methods
  elemental subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Finalize object.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(back_df_integrator), intent(INOUT) :: self !< BDF integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule foodie_integrator_backward_differentiation_formula
