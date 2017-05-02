!< Define [[oscillator]], the Oscillation test field that is a concrete extension of the abstract integrand type.

module oscillation_oscillator
!< Define [[oscillator]], the Oscillation test field that is a concrete extension of the abstract integrand type.

use foodie, only : integrand_object
use penf, only : R_P, I_P

implicit none
private
public :: oscillator

type, extends(integrand_object) :: oscillator
  !< Oscillator, the Oscillation equations field.
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
  !< State variables are organized as an array (rank 1) of reals of *n=2* elements.
  private
  real(R_P) :: f=0._R_P                !< Oscillation frequency (Hz).
  real(R_P) :: U(1:2)=[0._R_P, 0._R_P] !< Integrand (state) variables.
  contains
    ! auxiliary methods
    procedure, pass(self), public :: init   !< Init field.
    procedure, pass(self), public :: output !< Extract Oscillation field.
    ! public deferred methods
    procedure, pass(self), public :: t => doscillator_dt !< Time derivative, residuals.
    ! operators
    procedure, pass(lhs), public :: local_error !<`||oscillator - oscillator||` operator.
    ! +
    procedure, pass(lhs), public :: integrand_add_integrand !< `+` operator.
    procedure, pass(lhs), public :: integrand_add_real      !< `+ real` operator.
    procedure, pass(rhs), public :: real_add_integrand      !< `real +` operator.
    ! *
    procedure, pass(lhs), public :: integrand_multiply_integrand   !< `*` operator.
    procedure, pass(lhs), public :: integrand_multiply_real        !< `* real` operator.
    procedure, pass(rhs), public :: real_multiply_integrand        !< `real *` operator.
    procedure, pass(lhs), public :: integrand_multiply_real_scalar !< `* real_scalar` operator.
    procedure, pass(rhs), public :: real_scalar_multiply_integrand !< `real_scalar *` operator.
    ! -
    procedure, pass(lhs), public :: integrand_sub_integrand !< `-` operator.
    procedure, pass(lhs), public :: integrand_sub_real      !< `- real` operator.
    procedure, pass(rhs), public :: real_sub_integrand      !< `real -` operator.
    ! =
    procedure, pass(lhs), public :: assign_integrand !< `=` operator.
    procedure, pass(lhs), public :: assign_real      !< `= real` operator.
endtype oscillator

contains
  ! auxiliary methods
  pure subroutine init(self, initial_state, frequency)
  !< Construct an initialized oscillator field.
  class(oscillator), intent(inout) :: self               !< Oscillation field.
  real(R_P),         intent(in)    :: initial_state(1:2) !< Initial state of the Oscillation field vector.
  real(R_P),         intent(in)    :: frequency          !< Frequency of oscillation.

  self%U = initial_state
  self%f = frequency
  endsubroutine init

  pure function output(self) result(state)
  !< Output the oscillator field state.
  class(oscillator), intent(in) :: self       !< Oscillation field.
  real(R_P)                     :: state(1:2) !< Oscillation state vector.

  state = self%U
  endfunction output

  ! deferred methods
  pure function doscillator_dt(self, t) result(dState_dt)
  !< Time derivative of oscillator field.
  class(oscillator), intent(in)           :: self         !< Oscillation field.
  real(R_P),         intent(in), optional :: t            !< Time.
  real(R_P), allocatable                  :: dState_dt(:) !< Oscillation field time derivative.

  dState_dt = [-self%f * self%U(2), &
                self%f * self%U(1)]
  endfunction doscillator_dt

  pure function local_error(lhs, rhs) result(error)
  !< Estimate local truncation error between 2 oscillation approximations.
  !<
  !< The estimation is done by norm L2 of U:
  !<
  !< $$ error = \sqrt{ \sum_i{ \frac{(lhs\%U_i - rhs\%U_i)^2}{lhs\%U_i^2} }} $$
  class(oscillator),       intent(in) :: lhs   !< Left hand side.
  class(integrand_object), intent(in) :: rhs   !< Right hand side.
  real(R_P)                           :: error !< Error estimation.
  integer(I_P)                        :: i     !< Space counter.

  select type(rhs)
  class is(oscillator)
    error = 0._R_P
    do i=1, size(lhs%U, dim=1)
      error = error + (lhs%U(i) - rhs%U(i)) ** 2 / lhs%U(i) ** 2
    enddo
    error = sqrt(error)
  endselect
  endfunction local_error

  ! +
  pure function integrand_add_integrand(lhs, rhs) result(opr)
  !< `+` operator.
  class(oscillator),       intent(in) :: lhs    !< Left hand side.
  class(integrand_object), intent(in) :: rhs    !< Right hand side.
  real(R_P), allocatable              :: opr(:) !< Operator result.

  select type(rhs)
  class is(oscillator)
    opr = lhs%U + rhs%U
  endselect
  endfunction integrand_add_integrand

  pure function integrand_add_real(lhs, rhs) result(opr)
  !< `+ real` operator.
  class(oscillator), intent(in) :: lhs     !< Left hand side.
  real(R_P),         intent(in) :: rhs(1:) !< Right hand side.
  real(R_P), allocatable        :: opr(:)  !< Operator result.

  opr = lhs%U + rhs
  endfunction integrand_add_real

  pure function real_add_integrand(lhs, rhs) result(opr)
  !< `real +` operator.
  real(R_P),         intent(in) :: lhs(1:) !< Left hand side.
  class(oscillator), intent(in) :: rhs     !< Left hand side.
  real(R_P), allocatable        :: opr(:)  !< Operator result.

  opr = lhs + rhs%U
  endfunction real_add_integrand

  ! *
  pure function integrand_multiply_integrand(lhs, rhs) result(opr)
  !< `*` operator.
  class(oscillator),       intent(in) :: lhs    !< Left hand side.
  class(integrand_object), intent(in) :: rhs    !< Right hand side.
  real(R_P), allocatable              :: opr(:) !< Operator result.

  select type(rhs)
  class is(oscillator)
    opr = lhs%U * rhs%U
  endselect
  endfunction integrand_multiply_integrand

  pure function integrand_multiply_real(lhs, rhs) result(opr)
  !< `* real_scalar` operator.
  class(oscillator), intent(in) :: lhs     !< Left hand side.
  real(R_P),         intent(in) :: rhs(1:) !< Right hand side.
  real(R_P), allocatable        :: opr(:)  !< Operator result.

  opr = lhs%U * rhs
  endfunction integrand_multiply_real

  pure function real_multiply_integrand(lhs, rhs) result(opr)
  !< `real_scalar *` operator.
  class(oscillator), intent(in) :: rhs     !< Right hand side.
  real(R_P),         intent(in) :: lhs(1:) !< Left hand side.
  real(R_P), allocatable        :: opr(:)  !< Operator result.

  opr = lhs * rhs%U
  endfunction real_multiply_integrand

  pure function integrand_multiply_real_scalar(lhs, rhs) result(opr)
  !< `* real_scalar` operator.
  class(oscillator), intent(in) :: lhs    !< Left hand side.
  real(R_P),         intent(in) :: rhs    !< Right hand side.
  real(R_P), allocatable        :: opr(:) !< Operator result.

  opr = lhs%U * rhs
  endfunction integrand_multiply_real_scalar

  pure function real_scalar_multiply_integrand(lhs, rhs) result(opr)
  !< `real_scalar *` operator.
  real(R_P),         intent(in) :: lhs    !< Left hand side.
  class(oscillator), intent(in) :: rhs    !< Right hand side.
  real(R_P), allocatable        :: opr(:) !< Operator result.

  opr = lhs * rhs%U
  endfunction real_scalar_multiply_integrand

   ! -
  pure function integrand_sub_integrand(lhs, rhs) result(opr)
  !< `-` operator.
  class(oscillator),       intent(in) :: lhs    !< Left hand side.
  class(integrand_object), intent(in) :: rhs    !< Right hand side.
  real(R_P), allocatable              :: opr(:) !< Operator result.

  select type(rhs)
  class is(oscillator)
    opr = lhs%U - rhs%U
  endselect
  endfunction integrand_sub_integrand

  pure function integrand_sub_real(lhs, rhs) result(opr)
  !< `- real` operator.
  class(oscillator), intent(in) :: lhs     !< Left hand side.
  real(R_P),         intent(in) :: rhs(1:) !< Right hand side.
  real(R_P), allocatable        :: opr(:)  !< Operator result.

  opr = lhs%U - rhs
  endfunction integrand_sub_real

  pure function real_sub_integrand(lhs, rhs) result(opr)
  !< `real -` operator.
  class(oscillator), intent(in) :: rhs     !< Left hand side.
  real(R_P),         intent(in) :: lhs(1:) !< Left hand side.
  real(R_P), allocatable        :: opr(:)  !< Operator result.

  opr = lhs - rhs%U
  endfunction real_sub_integrand

  ! =
  pure subroutine assign_integrand(lhs, rhs)
  !< `=` operator.
  class(oscillator),       intent(inout) :: lhs !< Left hand side.
  class(integrand_object), intent(in)    :: rhs !< Right hand side.

  select type(rhs)
  class is(oscillator)
    lhs%U = rhs%U
    lhs%f = rhs%f
  endselect
  endsubroutine assign_integrand

  pure subroutine assign_real(lhs, rhs)
  !< `= real` operator.
  class(oscillator), intent(inout) :: lhs     !< Left hand side.
  real(R_P),         intent(in)    :: rhs(1:) !< Right hand side.

  lhs%U = rhs
  endsubroutine assign_real
endmodule oscillation_oscillator
