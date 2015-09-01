!< Define Oscillation field that is a concrete extension of the abstract integrand type.
module type_oscillation
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define Oscillation field that is a concrete extension of the abstract integrand type.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P
use foodie, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: oscillation
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(integrand) :: oscillation
  !< Oscillation equations field.
  !<
  !< It is a FOODiE integrand class.
  private
  real(R_P), dimension(:,:), allocatable :: state      !< Solution vector, [1:state_dims,1:time_steps_stored].
  real(R_P)                              :: f = 0._R_P !< Oscillation frequency (Hz).
  contains
    ! auxiliary methods
    procedure, pass(self), public :: init   !< Init field.
    procedure, pass(self), public :: output !< Extract Oscillation field.
    ! type_integrand deferred methods
    procedure, pass(self), public :: t => dOscillation_dt                                             !< Time derivative, residuals.
    procedure, pass(self), public :: update_previous_steps                                            !< Update previous time steps.
    procedure, pass(lhs),  public :: integrand_multiply_integrand => oscillation_multiply_oscillation !< Oscillation * oscillation.
    procedure, pass(lhs),  public :: integrand_multiply_real => oscillation_multiply_real             !< Oscillation * real.
    procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_oscillation             !< Real * Oscillation.
    procedure, pass(lhs),  public :: add => add_oscillation                                           !< Oscillation + Oscillation.
    procedure, pass(lhs),  public :: assign_integrand => oscillation_assign_oscillation               !< Oscillation = Oscillation.
    procedure, pass(lhs),  public :: assign_real => oscillation_assign_real                           !< Oscillation = real.
endtype oscillation
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine init(self, initial_state, f, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Construct an initialized Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation),      intent(INOUT) :: self          !< Oscillation field.
  real(R_P), dimension(:), intent(IN)    :: initial_state !< Intial state of the Oscillation field vector.
  real(R_P),               intent(IN)    :: f             !< Frequency.
  integer, optional,       intent(IN)    :: steps         !< Time steps stored.
  integer                                :: dsteps        !< Time steps stored, dummy variable.
  integer                                :: s             !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  dsteps = 1 ; if (present(steps)) dsteps = steps
  if (allocated(self%state)) deallocate(self%state) ; allocate(self%state(1:size(initial_state), 1:dsteps))
  do s=1, dsteps
    self%state(:, s) = initial_state
  enddo
  self%f = f
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Oscillation field state.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN)       :: self  !< Oscillation field.
  real(R_P), dimension(:), allocatable :: state !< Oscillation state vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  state = self%state(:, ubound(self%state, dim=2))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  pure function dOscillation_dt(self, n) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation),     intent(IN) :: self      !< Oscillation field.
  integer(I_P), optional, intent(IN) :: n         !< Time level.
  class(integrand),  allocatable     :: dState_dt !< Oscillation field time derivative.
  type(oscillation), allocatable     :: delta     !< Delta state used as temporary variables.
  integer                            :: dn        !< Time level, dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! preparing temporary delta
  allocate(delta)
  allocate(delta%state(1:size(self%state, dim=1), 1:size(self%state, dim=2)))
  ! Oscillation equations
  dn = size(self%state, dim=2) ; if (present(n)) dn = n
  delta%state(1, dn) = -self%f * self%state(2, dn)
  delta%state(2, dn) =  self%f * self%state(1, dn)
  ! hold Oscillation parameters constant over time
  delta%f = 0.
  call move_alloc (delta, dState_dt)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dOscillation_dt

  pure subroutine update_previous_steps(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Update previous time steps.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(INOUT) :: self !< Oscillation field.
  integer                           :: s    !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do s=1, ubound(self%state, dim=2) - 1
    self%state(:, s) = self%state(:, s + 1)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_previous_steps

  pure function oscillation_multiply_oscillation(lhs, rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a oscillation field by another one.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN) :: lhs           !< Left hand side.
  class(integrand),   intent(IN) :: rhs           !< Right hand side.
  class(integrand),  allocatable :: product       !< Product.
  type(oscillation), allocatable :: local_product !< Temporary product.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(local_product)
  select type(rhs)
  class is (oscillation)
    local_product%state = lhs%state * rhs%state
    local_product%f = lhs%f * rhs%f
  endselect
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction oscillation_multiply_oscillation

  pure function oscillation_multiply_real(lhs, rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a Oscillation field by a real scalar.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN) :: lhs           !< Left hand side.
  real(R_P),          intent(IN) :: rhs           !< Right hand side.
  class(integrand),  allocatable :: product       !< Product.
  type(oscillation), allocatable :: local_product !< Temporary product.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(local_product)
  local_product%state = lhs%state * rhs
  local_product%f = lhs%f * rhs
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction oscillation_multiply_real

  pure function real_multiply_oscillation(lhs, rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a real scalar by a Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P),          intent(IN) :: lhs           !< Left hand side.
  class(oscillation), intent(IN) :: rhs           !< Right hand side.
  class(integrand),  allocatable :: product       !< Product.
  type(oscillation), allocatable :: local_product !< Temporary product.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(local_product)
  local_product%state = rhs%state * lhs
  local_product%f     = rhs%f * lhs
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_multiply_oscillation

  pure function add_oscillation(lhs, rhs) result(sum)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Oscillation fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation),    intent(IN)  :: lhs       !< Left hand side.
  class(integrand), intent(IN)  :: rhs       !< Right hand side.
  class(integrand), allocatable :: sum       !< Sum.
  type(oscillation),     allocatable :: local_sum !< Temporary sum.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate (oscillation :: local_sum)
  select type(rhs)
  class is (oscillation)
    local_sum%state = lhs%state + rhs%state
    local_sum%f = lhs%f + rhs%f
  endselect
  call move_alloc(local_sum, sum)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction add_Oscillation

  pure subroutine oscillation_assign_oscillation(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one Oscillation field to another.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(INOUT) :: lhs !< Left hand side.
  class(integrand),   intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(rhs)
  class is (oscillation)
    if (allocated(rhs%state)) lhs%state = rhs%state
    lhs%f = rhs%f
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine oscillation_assign_oscillation

  pure subroutine oscillation_assign_real(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one real to a Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(INOUT) :: lhs !< Left hand side.
  real(R_P),          intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(lhs%state)) lhs%state = rhs
  lhs%f = rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine oscillation_assign_real
endmodule type_oscillation
