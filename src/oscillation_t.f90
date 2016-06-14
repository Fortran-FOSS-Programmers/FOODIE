!< Define Oscillation field that is a concrete extension of the abstract integrand type.
module oscillation_t
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define Oscillation field that is a concrete extension of the abstract integrand type.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie, only : integrand
use penf, only : R_P, I_P
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
  !< It is a FOODIE integrand class concrete extension.
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
  integer(I_P)                         :: dims=0   !< Space dimensions.
  real(R_P)                            :: f=0._R_P !< Oscillation frequency (Hz).
  real(R_P), dimension(:), allocatable :: U        !< Integrand (state) variables, [1:dims].
  contains
    ! auxiliary methods
    procedure, pass(self), public :: init   !< Init field.
    procedure, pass(self), public :: output !< Extract Oscillation field.
    ! ADT integrand deferred methods
    procedure, pass(self), public :: t => dOscillation_dt                                             !< Time derivative, residuals.
    procedure, pass(lhs),  public :: local_error => oscillation_local_error                           !<||Oscillation-oscillation||.
    procedure, pass(lhs),  public :: integrand_multiply_integrand => oscillation_multiply_oscillation !< Oscillation * oscillation.
    procedure, pass(lhs),  public :: integrand_multiply_real => oscillation_multiply_real             !< Oscillation * real.
    procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_oscillation             !< Real * Oscillation.
    procedure, pass(lhs),  public :: add => add_oscillation                                           !< Oscillation + Oscillation.
    procedure, pass(lhs),  public :: sub => sub_oscillation                                           !< Oscillation - Oscillation.
    procedure, pass(lhs),  public :: assign_integrand => oscillation_assign_oscillation               !< Oscillation = Oscillation.
endtype oscillation
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! auxiliary methods
  pure subroutine init(self, initial_state, frequency)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Construct an initialized Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation),      intent(INOUT) :: self          !< Oscillation field.
  real(R_P), dimension(:), intent(IN)    :: initial_state !< Initial state of the Oscillation field vector.
  real(R_P),               intent(IN)    :: frequency     !< Frequency of oscillation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%dims = size(initial_state)
  self%f = frequency
  if (allocated(self%U)) deallocate(self%U) ; allocate(self%U(1:self%dims)) ; self%U = initial_state
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

  ! ADT integrand deferred methods
#ifdef PURE
  pure function dOscillation_dt(self, t) result(dState_dt)
#else
  function dOscillation_dt(self, t) result(dState_dt)
#endif
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Oscillation field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation),     intent(IN) :: self      !< Oscillation field.
  real(R_P),    optional, intent(IN) :: t         !< Time.
  class(integrand),  allocatable     :: dState_dt !< Oscillation field time derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(oscillation :: dState_dt)
  select type(dState_dt)
  class is(oscillation)
    dState_dt = self
    dState_dt%U(1) = -self%f * self%U(2)
    dState_dt%U(2) =  self%f * self%U(1)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dOscillation_dt

#ifdef PURE
  pure function oscillation_local_error(lhs, rhs) result(error)
#else
  function oscillation_local_error(lhs, rhs) result(error)
#endif
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Estimate local truncation error between 2 oscillation approximations.
  !<
  !< The estimation is done by norm L2 of U:
  !<
  !< $$ error = \sqrt{ \sum_i{ \frac{(lhs\%U_i - rhs\%U_i)^2}{lhs\%U_i^2} }} $$
  !---------------------------------------------------------------------------------------------------------------------------------
  class(oscillation), intent(IN) :: lhs   !< Left hand side.
  class(integrand),   intent(IN) :: rhs   !< Right hand side.
  real(R_P)                      :: error !< Error estimation.
  integer(I_P)                   :: i     !< Space counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(rhs)
  class is (oscillation)
    error = 0._R_P
    do i=1, lhs%dims
      error = error + (lhs%U(i) - rhs%U(i))**2/lhs%U(i)**2
    enddo
    error = sqrt(error)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction oscillation_local_error

#ifdef PURE
  pure function oscillation_multiply_oscillation(lhs, rhs) result(opr)
#else
  function oscillation_multiply_oscillation(lhs, rhs) result(opr)
#endif
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

#ifdef PURE
  pure function oscillation_multiply_real(lhs, rhs) result(opr)
#else
  function oscillation_multiply_real(lhs, rhs) result(opr)
#endif
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

#ifdef PURE
  pure function real_multiply_oscillation(lhs, rhs) result(opr)
#else
  function real_multiply_oscillation(lhs, rhs) result(opr)
#endif
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

#ifdef PURE
  pure function add_oscillation(lhs, rhs) result(opr)
#else
  function add_oscillation(lhs, rhs) result(opr)
#endif
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

#ifdef PURE
  pure function sub_oscillation(lhs, rhs) result(opr)
#else
  function sub_oscillation(lhs, rhs) result(opr)
#endif
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

#ifdef PURE
  pure subroutine oscillation_assign_oscillation(lhs, rhs)
#else
  subroutine oscillation_assign_oscillation(lhs, rhs)
#endif
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
    lhs%f = rhs%f
    if (allocated(rhs%U)) lhs%U = rhs%U
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine oscillation_assign_oscillation
endmodule oscillation_t
