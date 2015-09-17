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
  !< It is a FOODiE integrand class concrete extension.
  !<
  !<### Burgers PDE equation
  !<The Burgers PDE equation is a non linear PDE widely used as numerical benchmark that can be applied to describe
  !<a wide range of different problems, from fluid dynamics to traffic flows, see [1].
  !<
  !<$$
  !<\begin{matrix}
  !<U_t = R(U)  \Leftrightarrow U_t = F(U)_x \\
  !<U = [u]\;\;\;
  !<F(U) = [-\frac{1}{2}u^2+\nu u_x]
  !<\end{matrix}
  !<$$
  !<
  !<This is the viscous Burgers' equation, \(\nu\) being the viscosity. This equation is of paramount relevance because it retains
  !<complex aspects such as the hyperbolic nature of the convection term (allowing discontinuous solutions) and the diffusive nature
  !<of the viscous term.
  !<
  !<#### Bibliography
  !<
  !<[1] *The partial differential equation ut + uux = nuxx*, Hopf, Eberhard, Communications on Pure and Applied Mathematics,
  !< vol 3, issue 3, doi 10.1002/cpa.3160030302, pp. 201--230, 1950.
  !<
  !<#### State variables organization
  !< State variable is organized as an array (rank 1) for whole physical domain.
  private
  integer(I_P)                           :: Ni=0      !< Number of grid nodes.
  integer(I_P)                           :: steps=0   !< Number of time steps stored.
  real(R_P)                              :: h=0._R_P  !< Space step discretization.
  real(R_P)                              :: nu=0._R_P !< Viscosity.
  real(R_P), dimension(:),   allocatable :: U         !< Integrand (state) variables, whole physical domain, [1:Ni].
  real(R_P), dimension(:,:), allocatable :: previous  !< Previous time steps states, [1:Ni,1:steps].
  contains
    ! auxiliary methods
    procedure, pass(self), public :: init             !< Init field.
    procedure, pass(self), public :: output           !< Extract Burgers field.
    procedure, pass(self), public :: dt => compute_dt !< Compute the current time step, by means of CFL condition.
    ! type_integrand deferred methods
    procedure, pass(self), public :: t => dBurgers_dt                                         !< Time derivative, residuals func.
    procedure, pass(self), public :: update_previous_steps                                    !< Update previous time steps.
    procedure, pass(self), public :: previous_step                                            !< Get a previous time step.
    procedure, pass(lhs),  public :: integrand_multiply_integrand => burgers_multiply_burgers !< Burgers * burgers operator.
    procedure, pass(lhs),  public :: integrand_multiply_real => burgers_multiply_real         !< Burgers * real operator.
    procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_burgers         !< Real * Burgers operator.
    procedure, pass(lhs),  public :: add => add_burgers                                       !< Burgers + Burgers operator.
    procedure, pass(lhs),  public :: sub => sub_burgers                                       !< Burgers - Burgers operator.
    procedure, pass(lhs),  public :: assign_integrand => burgers_assign_burgers               !< Burgers = Burgers.
    procedure, pass(lhs),  public :: assign_real => burgers_assign_real                       !< Burgers = real.
    ! private methods
    procedure, pass(self), private :: x  => dBurgers_dx   !< 1st derivative.
    procedure, pass(self), private :: xx => d2Burgers_dx2 !< 2nd derivative.
endtype burgers
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! auxiliary methods
  subroutine init(self, initial_state, Ni, h, nu, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Construct an initialized Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),             intent(INOUT) :: self          !< Burgers field.
  integer(I_P),               intent(IN)    :: Ni            !< Number of grid nodes.
  real(R_P), dimension(1:Ni), intent(IN)    :: initial_state !< Initial state of Burgers field domain.
  real(R_P),                  intent(IN)    :: h             !< Space step discretization.
  real(R_P),                  intent(IN)    :: nu            !< Viscosity.
  integer(I_P), optional,     intent(IN)    :: steps         !< Time steps stored.
  integer(I_P)                              :: s             !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = 0 ; if (present(steps)) self%steps = steps
  if (allocated(self%U)) deallocate(self%U) ; allocate(self%U(1:Ni))
  if (self%steps>0) then
    if (allocated(self%previous)) deallocate(self%previous) ; allocate(self%previous(1:Ni, 1:self%steps))
  endif
  self%U = initial_state
  if (self%steps>0) then
    do s=1, self%steps
      self%previous(:, s) = initial_state
    enddo
  endif
  self%Ni = Ni
  self%h = h
  self%nu = nu
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Burgers field state.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN)           :: self  !< Burgers field.
  real(R_P), dimension(:), allocatable :: state !< Burgers state variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  state = self%U
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  pure function compute_dt(self, CFL) result(dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the current time step, by means of CFL condition.
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

  ! type_integrand deferred methods
  pure function dBurgers_dt(self, n) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Burgers field, residuals function.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),         intent(IN) :: self      !< Burgers field.
  integer(I_P), optional, intent(IN) :: n         !< Time level.
  class(integrand), allocatable      :: dState_dt !< Burgers field time derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(burgers :: dState_dt)
  select type(dState_dt)
  class is(burgers)
    dState_dt = self
    dState_dt = self%xx(n=n) * self%nu
    dState_dt = dState_dt - self * self%x(n=n)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dBurgers_dt

  pure subroutine update_previous_steps(self, filter, weights)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Update previous time steps.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),             intent(INOUT) :: self       !< Burgers field.
  class(integrand), optional, intent(IN)    :: filter     !< Filter field displacement.
  real(R_P),        optional, intent(IN)    :: weights(:) !< Weights for filtering the steps.
  integer(I_P)                              :: s          !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (self%steps>0) then
    do s=1, self%steps - 1
      self%previous(:, s) = self%previous(:, s + 1)
    enddo
    self%previous(:, self%steps) = self%U
  endif
  if (present(filter).and.present(weights)) then
    select type(filter)
    class is(burgers)
      do s=1, self%steps
        self%previous(:, s) = self%previous(:, s) + filter%U * weights(s)
      enddo
    endselect
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_previous_steps

  pure function previous_step(self, n) result(previous)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Extract previous time solution of Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN)    :: self     !< Burgers field.
  integer(I_P),   intent(IN)    :: n        !< Time level.
  class(integrand), allocatable :: previous !< Previous time solution of Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(burgers :: previous)
  select type(previous)
  class is(burgers)
    previous = self
    previous%U = self%previous(:, n)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction previous_step

  pure function burgers_multiply_burgers(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a Burgers field by another one.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),   intent(IN)  :: lhs !< Left hand side.
  class(integrand), intent(IN)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(burgers :: opr)
  select type(opr)
  class is(burgers)
    opr = lhs
    select type(rhs)
    class is (burgers)
      opr%U = lhs%U * rhs%U
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction burgers_multiply_burgers

  pure function burgers_multiply_real(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a Burgers field by a real scalar.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers), intent(IN)    :: lhs !< Left hand side.
  real(R_P),      intent(IN)    :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(burgers :: opr)
  select type(opr)
  class is(burgers)
    opr = lhs
    opr%U = lhs%U * rhs
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction burgers_multiply_real

  pure function real_multiply_burgers(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a real scalar by a Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P),      intent(IN)    :: lhs !< Left hand side.
  class(burgers), intent(IN)    :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(burgers :: opr)
  select type(opr)
  class is(burgers)
    opr = rhs
    opr%U = rhs%U * lhs
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_multiply_burgers

  pure function add_burgers(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Burgers fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),   intent(IN)  :: lhs !< Left hand side.
  class(integrand), intent(IN)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(burgers :: opr)
  select type(opr)
  class is(burgers)
    opr = lhs
    select type(rhs)
    class is (burgers)
      opr%U = lhs%U + rhs%U
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction add_burgers

  pure function sub_burgers(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Subtract two Burgers fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),   intent(IN)  :: lhs !< Left hand side.
  class(integrand), intent(IN)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(burgers :: opr)
  select type(opr)
  class is(burgers)
    opr = lhs
    select type(rhs)
    class is (burgers)
      opr%U = lhs%U - rhs%U
    endselect
  endselect
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
    lhs%steps = rhs%steps
    lhs%h = rhs%h
    lhs%nu = rhs%nu
    if (allocated(rhs%U)) lhs%U = rhs%U
    if (allocated(rhs%previous)) lhs%previous = rhs%previous
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
  if (allocated(lhs%U)) lhs%U = rhs
  if (allocated(lhs%previous)) lhs%previous = rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine burgers_assign_real

  ! private methods
  pure function dBurgers_dx(self, n) result(derivative)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the first order spatial derivative of Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),    intent(IN) :: self       !< Burgers field.
  integer, optional, intent(IN) :: n          !< Time level.
  type(burgers)                 :: derivative !< Burgers field derivative.
  integer(I_P)                  :: i          !< Counter.
  integer                       :: dn         !< Time level, dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  derivative = self
  if (self%steps>=2) then ! self%previous should be used
    dn = self%steps ; if (present(n)) dn = n
    do i=2, self%Ni - 1
      derivative%U(i) = (self%previous(i+1, dn) - self%previous(i-1, dn))/(2._R_P * self%h)
    enddo
    derivative%U(1) = (self%previous(2, dn) - self%previous(self%Ni, dn))/(2._R_P * self%h)
    derivative%U(self%Ni) = (self%previous(1, dn) - self%previous(self%Ni-1, dn))/(2._R_P * self%h)
  else ! self%previous should not be used, use directly self%U
    do i=2, self%Ni - 1
      derivative%U(i) = (self%U(i+1) - self%U(i-1))/(2._R_P * self%h)
    enddo
    derivative%U(1) = (self%U(2) - self%U(self%Ni))/(2._R_P * self%h)
    derivative%U(self%Ni) = (self%U(1) - self%U(self%Ni-1))/(2._R_P * self%h)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction

  pure function d2Burgers_dx2(self, n) result(derivative)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the second order spatial derivative of Burgers field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(burgers),    intent(IN) :: self       !< Burgers field.
  integer, optional, intent(IN) :: n          !< Time level.
  type(burgers)                 :: derivative !< Burgers field derivative.
  integer(I_P)                  :: i          !< Counter.
  integer                       :: dn         !< Time level, dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  derivative = self
  if (self%steps>=2) then ! self%previous should be used
    dn = self%steps ; if (present(n)) dn = n
    do i=2, self%Ni - 1
      derivative%U(i) = (self%previous(i+1, dn) - 2._R_P * self%previous(i, dn) + self%previous(i-1, dn))/(self%h**2)
    enddo
    derivative%U(1) = (self%previous(2, dn) - 2._R_P * self%previous(1, dn) + self%previous(self%Ni, dn))/(self%h**2)
    derivative%U(self%Ni) = (self%previous(1, dn) - 2._R_P * self%previous(self%Ni, dn) + self%previous(self%Ni-1, dn))/&
                                (self%h**2)
  else ! self%previous should not be used, use directly self%U
    do i=2, self%Ni - 1
      derivative%U(i) = (self%U(i+1) - 2._R_P * self%U(i) + self%U(i-1))/(self%h**2)
    enddo
    derivative%U(1) = (self%U(2) - 2._R_P * self%U(1) + self%U(self%Ni))/(self%h**2)
    derivative%U(self%Ni) = (self%U(1) - 2._R_P * self%U(self%Ni) + self%U(self%Ni-1))/(self%h**2)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction
endmodule type_burgers
