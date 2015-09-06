!< Define Lorenz field that is a concrete extension of the abstract integrand type.
module type_lorenz
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define Lorenz field that is a concrete extension of the abstract integrand type.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P
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
  integer(I_P)                           :: dims=0       !< Space dimensions.
  integer(I_P)                           :: steps=0      !< Number of time steps stored.
  real(R_P), dimension(:,:), allocatable :: state        !< Solution vector, [1:dims,1:steps].
  real(R_P)                              :: sigma=0._R_P !< Lorenz \(\sigma\).
  real(R_P)                              :: rho=0._R_P   !< Lorenz \(\rho\).
  real(R_P)                              :: beta=0._R_P  !< Lorenz \(\beta\).
  contains
    ! auxiliary methods
    procedure, pass(self), public :: init   !< Init field.
    procedure, pass(self), public :: output !< Extract Lorenz field.
    ! type_integrand deferred methods
    procedure, pass(self), public :: t => dLorenz_dt                                        !< Time derivate, resiuduals function.
    procedure, pass(self), public :: update_previous_steps                                  !< Update previous time steps.
    procedure, pass(lhs),  public :: integrand_multiply_integrand => lorenz_multiply_lorenz !< Lorenz * lorenz operator.
    procedure, pass(lhs),  public :: integrand_multiply_real => lorenz_multiply_real        !< Lorenz * real operator.
    procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_lorenz        !< Real * Lorenz operator.
    procedure, pass(lhs),  public :: add => add_lorenz                                      !< Lorenz + Lorenz oprator.
    procedure, pass(lhs),  public :: assign_integrand => lorenz_assign_lorenz               !< Lorenz = Lorenz.
    procedure, pass(lhs),  public :: assign_real => lorenz_assign_real                      !< Lorenz = real.
endtype lorenz
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine init(self, initial_state, sigma, rho, beta, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Construct an initialized Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),           intent(INOUT) :: self          !< Lorenz field.
  real(R_P), dimension(:), intent(IN)    :: initial_state !< Intial state of Lorenz field vector.
  real(R_P),               intent(IN)    :: sigma         !< Lorenz  \(\sigma\).
  real(R_P),               intent(IN)    :: rho           !< Lorenz  \(\rho\).
  real(R_P),               intent(IN)    :: beta          !< Lorenz  \(\beta\).
  integer(I_P), optional,  intent(IN)    :: steps         !< Time steps stored.
  integer(I_P)                           :: s             !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%dims = size(initial_state)
  self%steps = 1 ; if (present(steps)) self%steps = steps
  if (allocated(self%state)) deallocate(self%state) ; allocate(self%state(1:self%dims, 1:self%steps))
  do s=1, self%steps
    self%state(:, s) = initial_state
  enddo
  self%sigma = sigma
  self%rho = rho
  self%beta = beta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Lorenz field state.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(IN)            :: self  !< Lorenz field.
  real(R_P), dimension(:), allocatable :: state !< Lorenz state vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  state = self%state(:, self%steps)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  pure function dLorenz_dt(self, n) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),          intent(IN) :: self      !< Lorenz field.
  integer(I_P), optional, intent(IN) :: n         !< Time level.
  class(integrand), allocatable      :: dState_dt !< Lorenz field time derivative.
  type(lorenz),     allocatable      :: delta     !< Delta state used as temporary variables.
  integer                            :: dn        !< Time level, dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! preparing temporary delta
  allocate(delta)
  delta%dims = self%dims
  delta%steps = self%steps
  allocate(delta%state(1:self%dims, 1:self%steps))
  delta%state = 0_R_P
  delta%sigma = self%sigma
  delta%rho = self%rho
  delta%beta = self%beta
  ! Lorenz equations
  dn = self%steps ; if (present(n)) dn = n
  delta%state(1, dn) = self%sigma * (self%state(2, dn) - self%state(1, dn))
  delta%state(2, dn) = self%state(1, dn) * (self%rho - self%state(3, dn)) - self%state(2, dn)
  delta%state(3, dn) = self%state(1, dn) * self%state(2, dn) - self%beta * self%state(3, dn)
  call move_alloc(delta, dState_dt)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dLorenz_dt

  pure subroutine update_previous_steps(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Update previous time steps.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(INOUT) :: self !< Lorenz field.
  integer                      :: s    !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do s=1, self%steps - 1
    self%state(:, s) = self%state(:, s + 1)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_previous_steps

  pure function lorenz_multiply_lorenz(lhs, rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a lorenz field by another one.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),    intent(IN)  :: lhs           !< Left hand side.
  class(integrand), intent(IN)  :: rhs           !< Right hand side.
  class(integrand), allocatable :: product       !< Product.
  type(lorenz),     allocatable :: local_product !< Temporary produtc.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(local_product)
  select type(rhs)
  class is (lorenz)
    local_product%dims = lhs%dims
    local_product%steps = lhs%steps
    local_product%state = lhs%state * rhs%state
    local_product%sigma = lhs%sigma
    local_product%rho = lhs%rho
    local_product%beta = lhs%beta
  endselect
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction lorenz_multiply_lorenz

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
  local_product%dims = lhs%dims
  local_product%steps = lhs%steps
  local_product%state = lhs%state * rhs
  local_product%sigma = lhs%sigma
  local_product%rho   = lhs%rho
  local_product%beta  = lhs%beta
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
  local_product%dims = rhs%dims
  local_product%steps = rhs%steps
  local_product%state = rhs%state * lhs
  local_product%sigma = rhs%sigma
  local_product%rho = rhs%rho
  local_product%beta = rhs%beta
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
    local_sum%dims = lhs%dims
    local_sum%steps = lhs%steps
    local_sum%state = lhs%state + rhs%state
    local_sum%sigma = lhs%sigma
    local_sum%rho = lhs%rho
    local_sum%beta = lhs%beta
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
    lhs%dims = rhs%dims
    lhs%steps = rhs%steps
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine lorenz_assign_real
endmodule type_lorenz
