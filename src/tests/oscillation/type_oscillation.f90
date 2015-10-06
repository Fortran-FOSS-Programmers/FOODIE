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
  !< It is a FOODiE integrand class concrete extension.
  !<
  !<### Oscillation ODEs system
  !<The (inertial) Oscillation equations system is a non linear system of pure ODEs and it can be written as:
  !<
  !<$$\begin{matrix}
  !< U_t = R(U)  \\
  !< U = \begin{bmatrix}
  !< v_1 \\
  !< v_2
  !< \end{bmatrix}\;\;\;
  !< R(U) = \begin{bmatrix}
  !< -f v_2 \\
  !< f v_1
  !< \end{bmatrix}
  !<\end{matrix}$$
  !<
  !<The frequency *f* is constant and it is here selected as *f=10^-4*.
  !<
  !<In the space *v1-v2* the path of the oscillation must be a circle, but for the frequency selected some ODE solvers are not
  !<stable leading to a wrong path, see [1].
  !<
  !<#### Bibliography
  !<
  !<[1] *Numerical Methods for Fluid Dynamics With Applications to Geophysics*, Dale R. Durran, Springer, 2010.
  !<
  !<#### State variables organization
  !< State variables are organized as an array (rank 1) of reals of *dims* elements, in this case 2 elements.
  private
  integer(I_P)                           :: dims=0   !< Space dimensions.
  real(R_P)                              :: f=0._R_P !< Oscillation frequency (Hz).
  real(R_P), dimension(:),   allocatable :: U        !< Integrand (state) variables, [1:dims].
  integer(I_P)                           :: steps=0  !< Number of time steps stored.
  real(R_P), dimension(:,:), allocatable :: previous !< Previous time steps states, [1:dims,1:steps].
  contains
    ! auxiliary methods
    procedure, pass(self), public :: init   !< Init field.
    procedure, pass(self), public :: output !< Extract Oscillation field.
    ! type_integrand deferred methods
    procedure, pass(self), public :: t => dOscillation_dt                                             !< Time derivative, residuals.
    procedure, pass(self), public :: update_previous_steps                                            !< Update previous time steps.
    procedure, pass(self), public :: previous_step                                                    !< Get a previous time step.
    procedure, pass(lhs),  public :: integrand_multiply_integrand => oscillation_multiply_oscillation !< Oscillation * oscillation.
    procedure, pass(lhs),  public :: integrand_multiply_real => oscillation_multiply_real             !< Oscillation * real.
    procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_oscillation             !< Real * Oscillation.
    procedure, pass(lhs),  public :: add => add_oscillation                                           !< Oscillation + Oscillation.
    procedure, pass(lhs),  public :: sub => sub_oscillation                                           !< Oscillation - Oscillation.
    procedure, pass(lhs),  public :: assign_integrand => oscillation_assign_oscillation               !< Oscillation = Oscillation.
    procedure, pass(lhs),  public :: assign_real => oscillation_assign_real                           !< Oscillation = real.
endtype oscillation
!-----------------------------------------------------------------------------------------------------------------------------------
contains
    ! auxiliary methods
  subroutine init(self, initial_state, f, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Construct an initialized Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation),      intent(INOUT) :: self          !< Oscillation field.
  real(R_P), dimension(:), intent(IN)    :: initial_state !< Initial state of the Oscillation field vector.
  real(R_P),               intent(IN)    :: f             !< Frequency.
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
  self%f = f
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  pure function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Oscillation field state.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN)       :: self  !< Oscillation field.
  real(R_P), dimension(:), allocatable :: state !< Oscillation state vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  state = self%U
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

    ! type_integrand deferred methods
  function dOscillation_dt(self, n, t) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation),     intent(IN) :: self      !< Oscillation field.
  integer(I_P), optional, intent(IN) :: n         !< Time level.
  real(R_P),    optional, intent(IN) :: t         !< Time.
  class(integrand),  allocatable     :: dState_dt !< Oscillation field time derivative.
  integer(I_P)                       :: dn        !< Time level, dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(oscillation :: dState_dt)
  select type(dState_dt)
  class is(oscillation)
    dState_dt = self
    if (self%steps>=2) then ! self%previous should be used
      dn = self%steps ; if (present(n)) dn = n
      dState_dt%U(1) = -self%f * self%previous(2, dn)
      dState_dt%U(2) =  self%f * self%previous(1, dn)
    else ! self%previous should not be used, use directly self%U
      dState_dt%U(1) = -self%f * self%U(2)
      dState_dt%U(2) =  self%f * self%U(1)
    endif
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dOscillation_dt

  subroutine update_previous_steps(self, filter, weights)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Update previous time steps.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation),         intent(INOUT) :: self       !< Oscillation field.
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
    class is(oscillation)
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
  !< Extract previous time solution of oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN) :: self     !< Oscillation field.
  integer(I_P),       intent(IN) :: n        !< Time level.
  class(integrand), allocatable  :: previous !< Previous time solution of oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(oscillation :: previous)
  select type(previous)
  class is(oscillation)
    previous = self
    previous%U = self%previous(:, n)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction previous_step

  function oscillation_multiply_oscillation(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a oscillation field by another one.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN) :: lhs !< Left hand side.
  class(integrand),   intent(IN) :: rhs !< Right hand side.
  class(integrand), allocatable  :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(oscillation :: opr)
  select type(opr)
  class is(oscillation)
    opr = lhs
    select type(rhs)
    class is (oscillation)
      opr%U = lhs%U * rhs%U
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction oscillation_multiply_oscillation

  function oscillation_multiply_real(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a Oscillation field by a real scalar.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN) :: lhs !< Left hand side.
  real(R_P),          intent(IN) :: rhs !< Right hand side.
  class(integrand), allocatable  :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(oscillation :: opr)
  select type(opr)
  class is(oscillation)
    opr = lhs
    opr%U = lhs%U * rhs
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction oscillation_multiply_real

  function real_multiply_oscillation(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a real scalar by a Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P),          intent(IN) :: lhs !< Left hand side.
  class(oscillation), intent(IN) :: rhs !< Right hand side.
  class(integrand), allocatable  :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(oscillation :: opr)
  select type(opr)
  class is(oscillation)
    opr = rhs
    opr%U = rhs%U * lhs
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_multiply_oscillation

  function add_oscillation(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Oscillation fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN) :: lhs !< Left hand side.
  class(integrand),   intent(IN) :: rhs !< Right hand side.
  class(integrand), allocatable  :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(oscillation :: opr)
  select type(opr)
  class is(oscillation)
    opr = lhs
    select type(rhs)
    class is (oscillation)
      opr%U = lhs%U + rhs%U
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction add_Oscillation

  function sub_oscillation(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Subtract two Oscillation fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN) :: lhs !< Left hand side.
  class(integrand),   intent(IN) :: rhs !< Right hand side.
  class(integrand), allocatable  :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(oscillation :: opr)
  select type(opr)
  class is(oscillation)
    opr = lhs
    select type(rhs)
    class is (oscillation)
      opr%U = lhs%U - rhs%U
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sub_Oscillation

  subroutine oscillation_assign_oscillation(lhs, rhs)
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
    if (allocated(rhs%U)) lhs%U = rhs%U
    if (allocated(rhs%previous)) lhs%previous = rhs%previous
    lhs%f = rhs%f
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine oscillation_assign_oscillation

  subroutine oscillation_assign_real(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one real to a Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(INOUT) :: lhs !< Left hand side.
  real(R_P),          intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(lhs%U)) lhs%U = rhs
  if (allocated(lhs%previous)) lhs%previous = rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine oscillation_assign_real
endmodule type_oscillation
