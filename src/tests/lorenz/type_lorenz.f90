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
  !< Lorenz equations field.
  !<
  !< It is a FOODiE integrand class.
  private
  real(R_P), dimension(:), allocatable :: state        !< Solution vector.
  real(R_P)                            :: sigma=0._R_P !< Lorenz \(\sigma\).
  real(R_P)                            :: rho=0._R_P   !< Lorenz \(\rho\).
  real(R_P)                            :: beta=0._R_P  !< Lorenz \(\beta\).
  contains
    procedure, pass(self), public :: output                                          !< Extract Lorenz field.
    procedure, pass(self), public :: t => dLorenz_dt                                 !< Time derivate, resiuduals function.
    procedure, pass(lhs),  public :: integrand_multiply_real => lorenz_multiply_real !< lorenz * real operator.
    procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_lorenz !< Real * Lorenz operator.
    procedure, pass(lhs),  public :: add => add_lorenz                               !< Lorenz + Lorenz oprator.
    procedure, pass(lhs),  public :: assign_integrand => lorenz_assign_lorenz        !< Lorenz = Lorenz.
    procedure, pass(lhs),  public :: assign_real => lorenz_assign_real               !< Lorenz = real.
endtype lorenz
interface lorenz
  !< Overload lorenz name adding the constructor function.
  module procedure constructor_lorenz
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  function constructor_lorenz(initial_state, sigma, rho, beta) result(concrete)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Construct an initialized Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), dimension(:), intent(IN)  :: initial_state !< Intial state of Lorenz field vector.
  real(R_P),               intent(IN)  :: sigma         !< Lorenz  \(\sigma\).
  real(R_P),               intent(IN)  :: rho           !< Lorenz  \(\rho\).
  real(R_P),               intent(IN)  :: beta          !< Lorenz  \(\beta\).
  type(lorenz)                         :: concrete      !< Concrete instance of Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  concrete%state = initial_state
  concrete%sigma = sigma
  concrete%rho = rho
  concrete%beta = beta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction constructor_lorenz

  function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Lorenz field state.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(IN)            :: self  !< Lorenz field.
  real(R_P), dimension(:), allocatable :: state !< Lorenz state vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  state = self%state
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  pure function dLorenz_dt(self) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(IN)     :: self      !< Lorenz field.
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

  pure function lorenz_multiply_real(lhs, rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a Lorenz field by a real scalar.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(IN)     :: lhs           !< Left hand side.
  real(R_P),     intent(IN)     :: rhs           !< Right hand side.
  class(integrand), allocatable :: product       !< Product.
  type(lorenz),     allocatable :: local_product !< Temporary produtc.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(local_product)
  local_product%state = lhs%state * rhs
  local_product%sigma = lhs%sigma * rhs
  local_product%rho   = lhs%rho   * rhs
  local_product%beta  = lhs%beta  * rhs
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction lorenz_multiply_real

  pure function real_multiply_lorenz(lhs, rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a real scalar by a Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P),     intent(IN)     :: lhs           !< Left hand side.
  class(lorenz), intent(IN)     :: rhs           !< Right hand side.
  class(integrand), allocatable :: product       !< Product.
  type(lorenz),     allocatable :: local_product !< Temporary produtc.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(local_product)
  local_product%state = rhs%state * lhs
  local_product%sigma = rhs%sigma * lhs
  local_product%rho   = rhs%rho   * lhs
  local_product%beta  = rhs%beta  * lhs
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_multiply_lorenz

  pure function add_lorenz(lhs, rhs) result(sum)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Lorenz fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),    intent(IN)  :: lhs       !< Left hand side.
  class(integrand), intent(IN)  :: rhs       !< Right hand side.
  class(integrand), allocatable :: sum       !< Sum.
  type(lorenz),     allocatable :: local_sum !< Temporary sum.
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
      ! stop 'add_Lorenz: rhs argument type not supported'
  endselect
  call move_alloc(local_sum, sum)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction add_Lorenz

  pure subroutine lorenz_assign_lorenz(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one Lorenz field to another.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),    intent(INOUT) :: lhs !< Left hand side.
  class(integrand), intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(rhs)
    class is (lorenz)
      if (allocated(rhs%state)) lhs%state = rhs%state
      lhs%sigma = rhs%sigma
      lhs%rho = rhs%rho
      lhs%beta = rhs%beta
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine lorenz_assign_lorenz

  pure subroutine lorenz_assign_real(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one real to a Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(INOUT) :: lhs !< Left hand side.
  real(R_P),     intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(lhs%state)) lhs%state = rhs
  lhs%sigma = rhs
  lhs%rho = rhs
  lhs%beta = rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine lorenz_assign_real
endmodule type_lorenz
