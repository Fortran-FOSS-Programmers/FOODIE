!< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property, from
!< 2nd to 3rd order accurate.

module foodie_integrator_lmm_ssp
!< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property, from
!< 2nd to 3rd order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the LMM-SSP class scheme implemented is:
!<
!< $$ U^n = \sum_{s=0}^{N_s-1}{\left[a_s U^{n-N_s+s} + \Delta t b_s \cdot R(t^{n-N_s+s}, U^{n-N_s+s}) \right]} $$
!<
!<where \(N_s\) is the number of previous steps considered.
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit. The coefficients *a,b* define the actual scheme, that is selected accordingly to the number of
!< **steps** used. Currently, the schemes provided have steps number in *[3, 5]*. The formal order of accuracy varies
!< consistently in *[2nd, 3rd]* order.
!<
!<#### Bibliography
!< [1] *Strong Stability Preserving Runge-Kutta and Multistep Time Discretizations*, S. Gottlieb, D. Ketcheson, C.W. Shu,
!< 2011, 978-981-4289-26-9, doi:10.1142/7498, World Scientific Publishing Co. Pte. Ltd.

use foodie_adt_integrand, only : integrand
use foodie_error_codes, only : ERROR_BAD_STEPS_NUMBER
use foodie_kinds, only : I_P, R_P
use foodie_integrator_object, only : integrator_object
use foodie_utils, only : is_admissible

implicit none
private
public :: integrator_lmm_ssp

character(len=99), parameter :: supported_steps='3-5' !< List of supported steps number. Valid format is `1-2,4,9-23...`.
integer(I_P),      parameter :: min_ss=3              !< Minimum number of steps supported.
integer(I_P),      parameter :: max_ss=5              !< Maximum number of steps supported.

type, extends(integrator_object) :: integrator_lmm_ssp
  !< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property,
  !< from 2nd to 3rd order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *a,b* coefficients) before used.
  private
  integer(I_P)           :: steps=0 !< Number of time steps.
  real(R_P), allocatable :: a(:)    !< *a* coefficients.
  real(R_P), allocatable :: b(:)    !< *b* coefficients.
  contains
    ! deferred methods
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    ! public methods
    procedure, pass(self) :: destroy         !< Destroy the integrator.
    procedure, pass(self) :: init            !< Initialize (create) the integrator.
    procedure, pass(self) :: integrate       !< Integrate integrand field.
    procedure, nopass     :: is_supported    !< Check if the queried number of steps is supported or not.
    procedure, nopass     :: min_steps       !< Return the minimum number of steps supported.
    procedure, nopass     :: max_steps       !< Return the maximum number of steps supported.
    procedure, pass(self) :: update_previous !< Cyclic update previous time steps.
endtype integrator_lmm_ssp

contains
  ! deferred methods
  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_lmm_ssp), intent(in)           :: self             !< Integrator.
  character(*),              intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                   :: desc             !< Description.
  character(len=:), allocatable                   :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                     :: NL=new_line('a') !< New line character.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = ''
  desc = desc//prefix_//'Strong Stability preserving Linear-Multistep-Methods class'//NL
  desc = desc//prefix_//'  Supported steps numbers: ['//trim(adjustl(supported_steps))//']'
  endfunction description

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_lmm_ssp), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),  intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  select type(rhs)
  class is (integrator_lmm_ssp)
                          lhs%steps = rhs%steps
    if (allocated(rhs%a)) lhs%a     = rhs%a
    if (allocated(rhs%b)) lhs%b     = rhs%b
  endselect
  endsubroutine integr_assign_integr

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_lmm_ssp), intent(inout) :: self !< Integrator.

  call self%destroy_abstract
  self%steps = 0
  if (allocated(self%a)) deallocate(self%a)
  if (allocated(self%b)) deallocate(self%b)
  endsubroutine destroy

  subroutine init(self, steps)
  !< Create the actual LMM-SSP integrator: initialize the *a,b* coefficients.
  !<
  !< @note If the integrator is initialized with a bad (unsupported) number of required time steps the initialization fails and
  !< the integrator error status is updated consistently for external-provided errors handling.
  class(integrator_lmm_ssp), intent(inout) :: self  !< Integrator.
  integer(I_P),              intent(in)    :: steps !< Number of time steps used.

  if (self%is_supported(steps)) then
    call self%destroy
    self%steps = steps
    allocate(self%a(1:steps)) ; self%a = 0.0_R_P
    allocate(self%b(1:steps)) ; self%b = 0.0_R_P
    select case(steps)
    case(3) ! LMM-SSP(3,2)
      self%a(1) = 1.0_R_P/4.0_R_P
      self%a(2) = 0.0_R_P
      self%a(3) = 3.0_R_P/4.0_R_P

      self%b(1) = 0.0_R_P
      self%b(2) = 0.0_R_P
      self%b(3) = 3.0_R_P/2.0_R_P
    case(4) ! LMM-SSP(4,3)
      self%a(1) = 11.0_R_P/27.0_R_P
      self%a(2) = 0.0_R_P
      self%a(3) = 0.0_R_P
      self%a(4) = 16.0_R_P/27.0_R_P

      self%b(1) = 12.0_R_P/27.0_R_P
      self%b(2) = 0.0_R_P
      self%b(3) = 0.0_R_P
      self%b(4) = 16.0_R_P/9.0_R_P
    case(5) ! LMM-SSP(5,3)
      self%a(1) = 7.0_R_P/32.0_R_P
      self%a(2) = 0.0_R_P
      self%a(3) = 0.0_R_P
      self%a(4) = 0.0_R_P
      self%a(5) = 25.0_R_P/32.0_R_P

      self%b(1) = 5.0_R_P/16.0_R_P
      self%b(2) = 0.0_R_P
      self%b(3) = 0.0_R_P
      self%b(4) = 0.0_R_P
      self%b(5) = 25.0_R_P/16.0_R_P
    endselect
  else
    call self%trigger_error(error=ERROR_BAD_STEPS_NUMBER,                           &
                            error_message='bad (unsupported) number of time steps', &
                            is_severe=.true.)
  endif
  endsubroutine init

  subroutine integrate(self, U, previous, Dt, t, autoupdate)
  !< Integrate field with LMM-SSP class scheme.
  class(integrator_lmm_ssp), intent(in)    :: self         !< Integrator.
  class(integrand),          intent(inout) :: U            !< Field to be integrated.
  class(integrand),          intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                 intent(in)    :: Dt           !< Time steps.
  real(R_P),                 intent(in)    :: t(:)         !< Times.
  logical, optional,         intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  logical                                  :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  integer(I_P)                             :: s            !< Steps counter.

  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  U = U * 0._R_P
  do s=1, self%steps
    U = U + previous(s) * self%a(s) + previous(s)%t(t=t(s)) * (Dt * self%b(s))
  enddo
  if (autoupdate_) call self%update_previous(U=U, previous=previous)
  endsubroutine integrate

  elemental function is_supported(steps)
  !< Check if the queried number of steps is supported or not.
  integer(I_P), intent(in) :: steps        !< Number of time steps used.
  logical                  :: is_supported !< Is true if the steps number is in *supported_steps*.

  is_supported = is_admissible(n=steps, adm_range=trim(supported_steps))
  endfunction is_supported

  pure function min_steps()
  !< Return the minimum number of steps supported.
  integer(I_P) :: min_steps !< Minimum number of steps supported.

  min_steps = min_ss
  endfunction min_steps

  pure function max_steps()
  !< Return the maximum number of steps supported.
  integer(I_P) :: max_steps !< Maximum number of steps supported.

  max_steps = max_ss
  endfunction max_steps

  subroutine update_previous(self, U, previous)
  !< Cyclic update previous time steps.
  class(integrator_lmm_ssp), intent(in)    :: self         !< Integrator.
  class(integrand),          intent(in)    :: U            !< Field to be integrated.
  class(integrand),          intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  integer(I_P)                             :: s            !< Steps counter.

  do s=1, self%steps - 1
    previous(s) = previous(s + 1)
  enddo
  previous(self%steps) = U
  endsubroutine update_previous
endmodule foodie_integrator_lmm_ssp
