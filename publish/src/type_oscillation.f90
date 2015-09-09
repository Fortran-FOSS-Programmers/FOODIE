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
  integer(I_P)                           :: dims=0     !< Space dimensions.
  integer(I_P)                           :: steps=0    !< Number of time steps stored.
  real(R_P), dimension(:),   allocatable :: state      !< Solution vector, [1:dims].
  real(R_P), dimension(:,:), allocatable :: previous   !< Previous steps solution vector, [1:dims,1:steps].
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
  integer(I_P), optional,  intent(IN)    :: steps         !< Time steps stored.
  integer(I_P)                           :: s             !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%dims = size(initial_state)
  self%steps = 0 ; if (present(steps)) self%steps = steps
  if (allocated(self%state)) deallocate(self%state) ; allocate(self%state(1:self%dims))
  if (self%steps>0) then
    if (allocated(self%previous)) deallocate(self%previous) ; allocate(self%previous(1:self%dims, 1:self%steps))
  endif
  self%state = initial_state
  if (self%steps>0) then
    do s=1, self%steps
      self%previous(:, s) = initial_state
    enddo
  endif
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
  state = self%state
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
  integer(I_P)                       :: dn        !< Time level, dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! preparing temporary delta
  allocate(delta)
  delta%dims = self%dims
  delta%steps = self%steps
  delta%state = self%state
  if (allocated(self%previous)) delta%previous = self%previous
  delta%f = self%f
  ! Oscillation equations
  if (self%steps>=2) then ! self%previous should be used
    dn = self%steps ; if (present(n)) dn = n
    delta%state(1) = -self%f * self%previous(2, dn)
    delta%state(2) =  self%f * self%previous(1, dn)
  else ! self%previous should not be used, use directly self%state
    delta%state(1) = -self%f * self%state(2)
    delta%state(2) =  self%f * self%state(1)
  endif
  call move_alloc(delta, dState_dt)
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
  if (self%steps>0) then
    do s=1, self%steps - 1
      self%previous(:, s) = self%previous(:, s + 1)
    enddo
    self%previous(:, self%steps) = self%state
  endif
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
    local_product%dims = lhs%dims
    local_product%steps = lhs%steps
    local_product%state = lhs%state * rhs%state
    if (allocated(lhs%previous)) local_product%previous = lhs%previous
    local_product%f = lhs%f
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
  local_product%dims = lhs%dims
  local_product%steps = lhs%steps
  local_product%state = lhs%state * rhs
  if (allocated(lhs%previous)) local_product%previous = lhs%previous
  local_product%f = lhs%f
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
  local_product%dims = rhs%dims
  local_product%steps = rhs%steps
  local_product%state = rhs%state * lhs
  if (allocated(rhs%previous)) local_product%previous = rhs%previous
  local_product%f     = rhs%f
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_multiply_oscillation

  pure function add_oscillation(lhs, rhs) result(sum)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Oscillation fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN) :: lhs       !< Left hand side.
  class(integrand),   intent(IN) :: rhs       !< Right hand side.
  class(integrand),  allocatable :: sum       !< Sum.
  type(oscillation), allocatable :: local_sum !< Temporary sum.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate (oscillation :: local_sum)
  select type(rhs)
  class is (oscillation)
    local_sum%dims = lhs%dims
    local_sum%steps = lhs%steps
    local_sum%state = lhs%state + rhs%state
    if (allocated(lhs%previous)) local_sum%previous = lhs%previous
    local_sum%f = lhs%f
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
    lhs%dims = rhs%dims
    lhs%steps = rhs%steps
    if (allocated(rhs%state)) lhs%state = rhs%state
    if (allocated(rhs%previous)) lhs%previous = rhs%previous
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
  if (allocated(lhs%previous)) lhs%previous = rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine oscillation_assign_real
endmodule type_oscillation
