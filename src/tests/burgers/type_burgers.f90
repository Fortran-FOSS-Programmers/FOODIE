!< Define Burgers field that is a concrete extension of the abstract integrand type.
module type_burgers
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define Burgers field that is a concrete extension of the abstract integrand type.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P
use foodie, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: burgers
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(integrand) :: burgers
  !< Burgers equations field.
  !<
  !< It is a FOODiE integrand class.
  private
  integer(I_P)                         :: Ni=0      !< Number of grid nodes.
  real(R_P)                            :: h=0._R_P  !< Space step discretization.
  real(R_P)                            :: nu=0._R_P !< Viscosity.
  real(R_P), dimension(:), allocatable :: state     !< Solution vector, whole physical domain, [1:Ni].
  contains
    ! public methods
    procedure, pass(self), public :: output                                                   !< Extract Burgers field.
    procedure, pass(self), public :: dt => compute_dt                                         !< Compute the current time step.
    procedure, pass(self), public :: t => dBurgers_dt                                         !< Time derivate, residuals function.
    procedure, pass(lhs),  public :: integrand_multiply_integrand => burgers_multiply_burgers !< burgers * burgers operator.
    procedure, pass(lhs),  public :: integrand_multiply_real => burgers_multiply_real         !< burgers * real operator.
    procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_burgers         !< Real * Burgers operator.
    procedure, pass(lhs),  public :: add => add_burgers                                       !< Burgers + Burgers oprator.
    procedure, pass(lhs),  public :: sub => sub_burgers                                       !< Burgers - Burgers oprator.
    procedure, pass(lhs),  public :: assign_integrand => burgers_assign_burgers               !< Burgers = Burgers.
    procedure, pass(lhs),  public :: assign_real => burgers_assign_real                       !< Burgers = real.
    ! operators overloading
    generic, public :: operator(*) => integrand_multiply_integrand, &
                                      integrand_multiply_real, &
                                      real_multiply_integrand
    generic, public :: operator(-) => sub
    ! private methods
    procedure, pass(self), private :: x  => dBurgers_dx   !< 1st derivative.
    procedure, pass(self), private :: xx => d2Burgers_dx2 !< 2nd derivative.
endtype burgers
interface burgers
  !< Overload burgers name adding the constructor function.
  module procedure constructor_burgers
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
function constructor_burgers(initial_state, Ni, h, nu) result(concrete)
!---------------------------------------------------------------------------------------------------------------------------------
!< Construct an initialized Burgers field.
!---------------------------------------------------------------------------------------------------------------------------------
real(R_P), dimension(:), intent(IN)  :: initial_state !< Intial state of Burgers field domain.
integer(I_P),            intent(IN)  :: Ni            !< Number of grid nodes.
real(R_P),               intent(IN)  :: h             !< Space step discretization.
real(R_P),               intent(IN)  :: nu            !< Viscosity.
type(burgers)                        :: concrete      !< Concrete instance of Burgers field.
!---------------------------------------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------------------------------------
concrete%state = initial_state
concrete%Ni = Ni
concrete%h = h
concrete%nu = nu
return
!---------------------------------------------------------------------------------------------------------------------------------
endfunction constructor_burgers

  function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Burgers field state.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN)           :: self  !< Burgers field.
  real(R_P), dimension(:), allocatable :: state !< Burgers state vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  state = self%state
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  pure function compute_dt(self, CFL) result(dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the current time step.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN) :: self !< Burgers field.
  real(R_P),      intent(IN) :: CFL  !< Courant-Friedricks-Lewi stability coefficient.
  real(R_P)                  :: dt   !< Current time step.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  dt = CFL * self%h**2 / self%nu
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_dt

  pure function dBurgers_dt(self) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Burgers field, residuals function.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN)    :: self      !< Burgers field.
  class(integrand), allocatable :: dState_dt !< Burgers field time derivative.
  type(burgers),    allocatable :: delta     !< Delta state used as temporary variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! preparing temporary variables
  allocate(delta, source=self)
  ! Burgers residuals
  delta = self%xx() * self%nu
  delta = delta - self * self%x()
  call move_alloc (delta, dState_dt)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dBurgers_dt

  pure function burgers_multiply_burgers(lhs, rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a Burgers field by another one.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN)    :: lhs           !< Left hand side.
  type(burgers),  intent(IN)    :: rhs           !< Right hand side.
  class(integrand), allocatable :: product       !< Product.
  type(burgers),    allocatable :: local_product !< Temporary produtc.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(local_product, source=lhs)
  ! local_product%Ni = lhs%Ni
  ! local_product%h = lhs%h
  ! local_product%nu = lhs%nu
  local_product%state = lhs%state * rhs%state
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction burgers_multiply_burgers

  pure function burgers_multiply_real(lhs, rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a Burgers field by a real scalar.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN)    :: lhs           !< Left hand side.
  real(R_P),      intent(IN)    :: rhs           !< Right hand side.
  class(integrand), allocatable :: product       !< Product.
  type(burgers),    allocatable :: local_product !< Temporary produtc.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(local_product, source=lhs)
  local_product%state = lhs%state * rhs
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction burgers_multiply_real

  pure function real_multiply_burgers(lhs, rhs) result(product)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a real scalar by a Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P),      intent(IN)    :: lhs           !< Left hand side.
  class(burgers), intent(IN)    :: rhs           !< Right hand side.
  class(integrand), allocatable :: product       !< Product.
  type(burgers),    allocatable :: local_product !< Temporary produtc.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(local_product, source=rhs)
  local_product%state = rhs%state * lhs
  call move_alloc(local_product, product)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_multiply_burgers

  pure function add_burgers(lhs, rhs) result(sum)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Burgers fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),   intent(IN)  :: lhs       !< Left hand side.
  class(integrand), intent(IN)  :: rhs       !< Right hand side.
  class(integrand), allocatable :: sum       !< Sum.
  type(burgers),    allocatable :: local_sum !< Temporary sum.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate (local_sum, source=lhs)
  select type(rhs)
    class is (burgers)
      local_sum%state = lhs%state + rhs%state
  endselect
  call move_alloc(local_sum, sum)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction add_burgers

  pure function sub_burgers(lhs, rhs) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Subtract two Burgers fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),   intent(IN)  :: lhs       !< Left hand side.
  class(integrand), intent(IN)  :: rhs       !< Right hand side.
  class(integrand), allocatable :: sub       !< Subtraction.
  type(burgers),    allocatable :: local_sub !< Temporary subtraction.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate (local_sub, source=lhs)
  select type(rhs)
    class is (burgers)
      local_sub%state = lhs%state - rhs%state
  endselect
  call move_alloc(local_sub, sub)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sub_burgers

  pure subroutine burgers_assign_burgers(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one Burgers field to another.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),   intent(INOUT) :: lhs !< Left hand side.
  class(integrand), intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(rhs)
    class is (burgers)
      lhs%Ni = rhs%Ni
      lhs%h = rhs%h
      lhs%nu = rhs%nu
      if (allocated(rhs%state)) lhs%state = rhs%state
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine burgers_assign_burgers

  pure subroutine burgers_assign_real(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one real to a Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(INOUT) :: lhs !< Left hand side.
  real(R_P),      intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(lhs%state)) lhs%state = rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine burgers_assign_real

  ! private methods
  pure function dBurgers_dx(self) result(derivative)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the first order spatial derivative of Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN) :: self       !< Burgers field.
  type(burgers), allocatable :: derivative !< Burgers field derivative.
  integer(I_P)               :: i          !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(derivative, source=self)
  do i=2, self%Ni - 1
    derivative%state(i) = (self%state(i+1) - self%state(i-1))/(2._R_P * self%h)
  enddo
  derivative%state(1) = (self%state(2) - self%state(self%Ni))/(2._R_P * self%h)
  derivative%state(self%Ni) = (self%state(1) - self%state(self%Ni-1))/(2._R_P * self%h)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction

  pure function d2Burgers_dx2(self) result(derivative)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the second order spatial derivative of Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN) :: self       !< Burgers field.
  type(burgers), allocatable :: derivative !< Burgers field derivative.
  integer(I_P)               :: i          !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(derivative, source=self)
  do i=2, self%Ni - 1
    derivative%state(i) = (self%state(i+1) - 2._R_P * self%state(i) + self%state(i-1))/(self%h**2)
  enddo
  derivative%state(1) = (self%state(2) - 2._R_P * self%state(1) + self%state(self%Ni))/(self%h**2)
  derivative%state(self%Ni) = (self%state(1) - 2._R_P * self%state(self%Ni) + self%state(self%Ni-1))/(self%h**2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction
endmodule type_burgers
