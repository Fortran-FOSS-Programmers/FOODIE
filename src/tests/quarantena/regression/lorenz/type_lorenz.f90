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
  !< It is a FOODIE integrand class concrete extension.
  !<
  !<### Lorenz ODEs system
  !<The Lorenz' equations system [1] is a non linear system of pure ODEs that retains a reasonable-complex behaviour: such a
  !<system, for a certain parameters-region exhibits a chaotic dynamics useful for testing FOODIE solvers.
  !<
  !<The Lorenz' ODEs system can be written as:
  !<
  !<$$\begin{matrix}
  !< U_t = R(U)  \\
  !< U = \begin{bmatrix}
  !< v_1 \\
  !< v_2 \\
  !< v_3
  !< \end{bmatrix}\;\;\;
  !< R(U) = \begin{bmatrix}
  !< \sigma (v_2-v_1) \\
  !< v_1(\rho - v_3) -v_2 \\
  !< v_1 v_2 - \beta v_3
  !< \end{bmatrix}
  !<\end{matrix}$$
  !<
  !<The parameters set is constant and it is here selected as:
  !<
  !<$$\begin{matrix}
  !< \sigma = 10 \\
  !< \rho = 28 \\
  !< \beta = \frac{8}{3}
  !<\end{matrix}$$
  !<
  !<These values are chaos-inducing thus they magnify the eventual numerical inaccuracies of FOODIE solvers, see [2].
  !<
  !<#### Bibliography
  !<
  !<[1] *Deterministic Nonperiodic Flow*, Lorenz E.N., Journal of the Atmospheric Sciences, 1963, vol. 20, pp. 130--141,
  !<doi: http://dx.doi.org/10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2
  !<
  !<[2] *Scientific software design: the object-oriented way*, Rouson, Damian, Jim Xia, and Xiaofeng Xu,
  !<Cambridge University Press, 2011
  !<
  !<#### State variables organization
  !< State variables are organized as an array (rank 1) of reals of *dims* elements, in this case 3 elements.
  private
  integer(I_P)                           :: dims=0       !< Space dimensions.
  integer(I_P)                           :: steps=0      !< Number of time steps stored.
  real(R_P), dimension(:),   allocatable :: U            !< Integrand (state) variables, [1:dims].
  real(R_P), dimension(:,:), allocatable :: previous     !< Previous time steps states, [1:dims,1:steps].
  real(R_P)                              :: sigma=0._R_P !< Lorenz \(\sigma\).
  real(R_P)                              :: rho=0._R_P   !< Lorenz \(\rho\).
  real(R_P)                              :: beta=0._R_P  !< Lorenz \(\beta\).
  contains
    ! auxiliary methods
    procedure, pass(self), public :: init   !< Init field.
    procedure, pass(self), public :: output !< Extract Lorenz field.
    ! ADT integrand deferred methods
    procedure, pass(self), public :: t => dLorenz_dt                                        !< Time derivative, residuals function.
    procedure, pass(lhs),  public :: local_error => lorenz_local_error                      !< Local error.
    procedure, pass(self), public :: update_previous_steps                                  !< Update previous time steps.
    procedure, pass(self), public :: previous_step                                          !< Get a previous time step.
    procedure, pass(lhs),  public :: integrand_multiply_integrand => lorenz_multiply_lorenz !< Lorenz * Lorenz operator.
    procedure, pass(lhs),  public :: integrand_multiply_real => lorenz_multiply_real        !< Lorenz * real operator.
    procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_lorenz        !< Real * Lorenz operator.
    procedure, pass(lhs),  public :: add => add_lorenz                                      !< Lorenz + Lorenz operator.
    procedure, pass(lhs),  public :: sub => sub_lorenz                                      !< Lorenz - Lorenz.
    procedure, pass(lhs),  public :: assign_integrand => lorenz_assign_lorenz               !< Lorenz = Lorenz.
    procedure, pass(lhs),  public :: assign_real => lorenz_assign_real                      !< Lorenz = real.
endtype lorenz
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! auxiliary methods
  subroutine init(self, initial_state, sigma, rho, beta, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Construct an initialized Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),           intent(INOUT) :: self          !< Lorenz field.
  real(R_P), dimension(:), intent(IN)    :: initial_state !< Initial state of Lorenz field vector.
  real(R_P),               intent(IN)    :: sigma         !< Lorenz  \(\sigma\).
  real(R_P),               intent(IN)    :: rho           !< Lorenz  \(\rho\).
  real(R_P),               intent(IN)    :: beta          !< Lorenz  \(\beta\).
  integer(I_P), optional,  intent(IN)    :: steps         !< Time steps stored.
  integer(I_P)                           :: s             !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%dims = size(initial_state)
  self%steps = 0 ; if (present(steps)) self%steps = steps
  if (allocated(self%U)) deallocate(self%U) ; allocate(self%U(1:self%dims))
  if (self%steps>0) then
    if (allocated(self%previous)) deallocate(self%previous) ; allocate(self%previous(1:self%dims, 1:self%steps))
  endif
  self%U = initial_state
  if (self%steps>0) then
    do s=1, self%steps
      self%previous(:, s) = initial_state
    enddo
  endif
  self%sigma = sigma
  self%rho = rho
  self%beta = beta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  pure function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Lorenz field state.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(IN)            :: self  !< Lorenz field.
  real(R_P), dimension(:), allocatable :: state !< Lorenz state vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  state = self%U
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  ! ADT integrand deferred methods
  function dLorenz_dt(self, t) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),          intent(IN) :: self      !< Lorenz field.
  real(R_P),    optional, intent(IN) :: t         !< Time.
  class(integrand), allocatable      :: dState_dt !< Lorenz field time derivative.
  integer(I_P)                       :: dn        !< Time level, dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(lorenz :: dState_dt)
  select type(dState_dt)
  class is(lorenz)
    dState_dt = self
    dState_dt%U(1) = self%sigma * (self%U(2) - self%U(1))
    dState_dt%U(2) = self%U(1) * (self%rho - self%U(3)) - self%U(2)
    dState_dt%U(3) = self%U(1) * self%U(2) - self%beta * self%U(3)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dLorenz_dt

  subroutine update_previous_steps(self, filter, weights)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Update previous time steps.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),              intent(INOUT) :: self       !< Lorenz field.
  class(integrand), optional, intent(IN)    :: filter     !< Filter field displacement.
  real(R_P),        optional, intent(IN)    :: weights(:) !< Weights for filtering the steps.
  integer                                   :: s          !< Time steps counter.
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
    class is(lorenz)
      do s=1, self%steps
        self%previous(:, s) = self%previous(:, s) + filter%U * weights(s)
      enddo
    endselect
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_previous_steps

  function previous_step(self, n) result(previous)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Extract previous time solution of Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(IN)     :: self     !< Lorenz field.
  integer(I_P),  intent(IN)     :: n        !< Time level.
  class(integrand), allocatable :: previous !< Previous time solution of Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(lorenz :: previous)
  select type(previous)
  class is(lorenz)
    previous = self
    previous%U = self%previous(:, n)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction previous_step

  function lorenz_local_error(lhs, rhs) result(error)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Estimate local truncation error between 2 lorenz approximations.
  !<
  !< The estimation is done by norm L2 of U:
  !<
  !< $$ error = \sqrt{ \sum_i{ \frac{(lhs\%U_i - rhs\%U_i)^2}{lhs\%U_i^2} }} $$
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(IN) :: lhs   !< Left hand side.
  class(integrand),   intent(IN) :: rhs   !< Right hand side.
  real(R_P)                      :: error !< Error estimation.
  integer(I_P)                   :: i     !< Space counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(rhs)
  class is (lorenz)
    error = 0._R_P
    do i=1, lhs%dims
      error = error + (lhs%U(i) - rhs%U(i))**2/lhs%U(i)**2
    enddo
    error = sqrt(error)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction lorenz_local_error


  function lorenz_multiply_lorenz(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a lorenz field by another one.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),    intent(IN)  :: lhs !< Left hand side.
  class(integrand), intent(IN)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(lorenz :: opr)
  select type(opr)
  class is(lorenz)
    opr = lhs
    select type(rhs)
    class is (lorenz)
      opr%U = lhs%U * rhs%U
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction lorenz_multiply_lorenz

  function lorenz_multiply_real(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a Lorenz field by a real scalar.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(IN)     :: lhs !< Left hand side.
  real(R_P),     intent(IN)     :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(lorenz :: opr)
  select type(opr)
  class is(lorenz)
    opr = lhs
    opr%U = lhs%U * rhs
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction lorenz_multiply_real

  function real_multiply_lorenz(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a real scalar by a Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P),     intent(IN)     :: lhs !< Left hand side.
  class(lorenz), intent(IN)     :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(lorenz :: opr)
  select type(opr)
  class is(lorenz)
    opr = rhs
    opr%U = rhs%U * lhs
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_multiply_lorenz

  function add_lorenz(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Lorenz fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),    intent(IN)  :: lhs !< Left hand side.
  class(integrand), intent(IN)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(lorenz :: opr)
  select type(opr)
  class is(lorenz)
    opr = lhs
    select type(rhs)
    class is (lorenz)
      opr%U = lhs%U + rhs%U
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction add_Lorenz

  function sub_lorenz(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Subtract two Lorenz fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz),    intent(IN)  :: lhs !< Left hand side.
  class(integrand), intent(IN)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(lorenz :: opr)
  select type(opr)
  class is(lorenz)
    opr = lhs
    select type(rhs)
    class is (lorenz)
      opr%U = lhs%U - rhs%U
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sub_lorenz

  subroutine lorenz_assign_lorenz(lhs, rhs)
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
    if (allocated(rhs%U)) lhs%U = rhs%U
    if (allocated(rhs%previous)) lhs%previous = rhs%previous
    lhs%sigma = rhs%sigma
    lhs%rho = rhs%rho
    lhs%beta = rhs%beta
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine lorenz_assign_lorenz

  subroutine lorenz_assign_real(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one real to a Lorenz field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(lorenz), intent(INOUT) :: lhs !< Left hand side.
  real(R_P),     intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(lhs%U)) lhs%U = rhs
  if (allocated(lhs%previous)) lhs%previous = rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine lorenz_assign_real
endmodule type_lorenz
