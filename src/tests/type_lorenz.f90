!< Define Lorenz field that is a concrete extension of the abstract integrand type.
module type_lorenz
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define Lorenz field that is a concrete extension of the abstract integrand type.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P
use foodie, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: lorenz
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(integrand) :: lorenz
  private
  real(R_P), dimension(:), allocatable :: state            !< Solution vector.
  real(R_P)                            :: sigma, rho, beta !< Lorenz parameters.
contains
   procedure, pass(self), public :: t        => dLorenz_dt
   procedure, pass(lhs),  public :: add      => add_lorenz
   procedure, pass(lhs),  public :: multiply => multiply_lorenz
   procedure, pass(lhs),  public :: assign   => assign_lorenz
   procedure, pass(self), public :: output
end type
interface lorenz
  module procedure constructor_lorenz
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  type(lorenz) function constructor_lorenz(initial_state, sigma, rho, beta)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Construct an initialized Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), dimension(:), intent(in)  :: initial_state    !< Intial state of Lorenz field vector.
  real(R_P),               intent(in)  :: sigma, rho, beta !< Lorenz parameters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  constructor_lorenz%state = initial_state
  constructor_lorenz%sigma = sigma
  constructor_lorenz%rho = rho
  constructor_lorenz%beta = beta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction constructor_lorenz

  function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Lorenz field state.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(in)            :: self  !< Lorenz field.
  real(R_P), dimension(:), allocatable :: state !< Lorenz state vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  state = self%state
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  function dLorenz_dt(self) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(in)     :: self      !< Lorenz field.
  class(integrand), allocatable :: dState_dt !< Lorenz field time derivative.
  type(lorenz),     allocatable :: delta     !< Delta state used as temporary variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! preparing temporary delta
  allocate(delta)
  allocate(delta%state(size(self%state)))
  ! Lorenz equations
  delta%state(1) = self%sigma * (self%state(2) - self%state(1))
  delta%state(2) = self%state(1) * (self%rho - self%state(3)) - self%state(2)
  delta%state(3) = self%state(1) * self%state(2) - self%beta * self%state(3)
  ! hold Lorenz parameters constant over time
  delta%sigma = 0.
  delta%rho = 0.
  delta%beta = 0.
  call move_alloc (delta, dState_dt)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dLorenz_dt

  function add_Lorenz(lhs, rhs) result(sum)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Lorenz fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),    intent(in)  :: lhs
  class(integrand), intent(in)  :: rhs
  class(integrand), allocatable :: sum
  type(lorenz),     allocatable :: local_sum
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate (lorenz :: local_sum)
  select type(rhs)
    class is (lorenz)
      local_sum%state = lhs%state + rhs%state
      local_sum%sigma = lhs%sigma + rhs%sigma
      local_sum%rho   = lhs%rho   + rhs%rho
      local_sum%beta  = lhs%beta  + rhs%beta
    class default
      stop 'add_Lorenz: rhs argument type not supported'
  endselect
  call move_alloc(local_sum, sum)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction add_Lorenz

  function multiply_Lorenz(lhs,rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a Lorenz field by a real scalar
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(in)     :: lhs
  real(R_P),     intent(in)     :: rhs
  class(integrand), allocatable :: product
  type(lorenz),     allocatable :: local_product
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate (local_product)
  local_product%state = lhs%state * rhs
  local_product%sigma = lhs%sigma * rhs
  local_product%rho   = lhs%rho   * rhs
  local_product%beta  = lhs%beta  * rhs
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction multiply_Lorenz

  subroutine assign_lorenz(lhs,rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one Lorenz field to another.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),    intent(inout) :: lhs
  class(integrand), intent(in)    :: rhs
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(rhs)
    class is (lorenz)
      lhs%state = rhs%state
      lhs%sigma = rhs%sigma
      lhs%rho   = rhs%rho
      lhs%beta  = rhs%beta
    class default
      stop 'assign_lorenz: rhs argument type not supported'
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine
endmodule type_lorenz
