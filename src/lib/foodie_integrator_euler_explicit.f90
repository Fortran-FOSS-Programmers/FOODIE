!< FOODIE integrator: provide explicit Euler scheme, it being 1st order accurate.

module foodie_integrator_euler_explicit
!< FOODIE integrator: provide explicit Euler scheme, it being 1st order accurate.
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

use foodie_adt_integrand, only : integrand
use foodie_kinds, only : I_P, R_P
use foodie_integrator_object, only : integrator_object
use foodie_utils, only : is_admissible

implicit none
private
public :: integrator_euler_explicit

character(len=1), parameter :: supported_stages_steps='1' !< List of supported stages/steps number. Valid format is `1-2,4,9-23...`.
integer(I_P),     parameter :: min_ss=1                   !< Minimum number of stages/steps supported.
integer(I_P),     parameter :: max_ss=1                   !< Maximum number of stages/steps supported.

type, extends(integrator_object) :: integrator_euler_explicit
  !< FOODIE integrator: provide explicit Euler scheme, it being 1st order accurate.
  !<
  !< @note The integrator can be used directly without any initialization.
  contains
    ! deferred methods
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    ! public methods
    procedure, pass(self) :: destroy          !< Destroy the integrator.
    procedure, nopass     :: integrate        !< Integrate integrand field.
    procedure, nopass     :: is_supported     !< Check if the queried number of stages/steps is supported or not.
    procedure, nopass     :: min_stages_steps !< Return the minimum number of stages/steps supported.
    procedure, nopass     :: max_stages_steps !< Return the maximum number of stages/steps supported.
endtype integrator_euler_explicit
contains
  ! deferred methods
  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_euler_explicit), intent(in)           :: self    !< Integrator.
  character(*),                     intent(in), optional :: prefix  !< Prefixing string.
  character(len=:), allocatable                          :: desc    !< Description.
  character(len=:), allocatable                          :: prefix_ !< Prefixing string, local variable.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = prefix_//'Euler, Explicit (1 step/stage) 1st order scheme'
  endfunction description

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_euler_explicit), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),         intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  endsubroutine integr_assign_integr

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_euler_explicit), intent(inout) :: self !< Integrator.

  call self%destroy_abstract
  endsubroutine destroy

  subroutine integrate(U, Dt, t)
  !< Integrate field with explicit Euler scheme, 1st order.
  class(integrand),    intent(inout) :: U  !< Field to be integrated.
  real(R_P),           intent(in)    :: Dt !< Time step.
  real(R_P), optional, intent(in)    :: t  !< Time.

  U = U + U%t(t=t) * Dt
  endsubroutine integrate

  elemental function is_supported(stages_steps)
  !< Check if the queried number of stages/steps is supported or not.
  integer(I_P), intent(in) :: stages_steps !< Number of stages/steps used.
  logical                  :: is_supported !< Is true is the stages number is in *supported_stages_steps*.

  is_supported = is_admissible(n=stages_steps, adm_range=trim(supported_stages_steps))
  endfunction is_supported

  pure function min_stages_steps()
  !< Return the minimum number of stages/steps supported.
  integer(I_P) :: min_stages_steps !< Minimum number of stages/steps supported.

  min_stages_steps = min_ss
  endfunction min_stages_steps

  pure function max_stages_steps()
  !< Return the maximum number of stages/steps supported.
  integer(I_P) :: max_stages_steps !< Maximum number of stages/steps supported.

  max_stages_steps = max_ss
  endfunction max_stages_steps
endmodule foodie_integrator_euler_explicit
