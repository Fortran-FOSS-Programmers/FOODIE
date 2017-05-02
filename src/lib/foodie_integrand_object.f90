!< Define the abstract type *integrand* for building FOODIE ODE integrators.

module foodie_integrand_object
!< Define the abstract type *integrand* for building FOODIE ODE integrators.

use penf, only : I_P, R_P

implicit none
private
public :: integrand_object

type, abstract :: integrand_object
  !< Abstract type for building FOODIE ODE integrators.
#ifdef CAF
  class(*), allocatable :: dummy_to_allow_extensions[:] !< Dummy member to allow concrete extensions with coarray members.
#endif
  contains
    ! public deferred procedures that concrete integrand-field must implement
    procedure(time_derivative), pass(self), deferred, public :: t !< Time derivative, residuals.
    ! operators
    procedure(local_error_operator), pass(lhs), deferred, public :: local_error !< `||integrand - integrand||` operator.
    generic, public :: operator(.lterror.) => local_error !< Estimate local truncation error.
    ! +
    procedure(symmetric_operator), pass(lhs), deferred, public :: integrand_add_integrand !< `+` operator.
    procedure(integrand_op_real),  pass(lhs), deferred, public :: integrand_add_real      !< `+ real` operator.
    procedure(real_op_integrand),  pass(rhs), deferred, public :: real_add_integrand      !< `real +` operator.
    generic, public :: operator(+) => integrand_add_integrand, &
                                      integrand_add_real,      &
                                      real_add_integrand !< Overloading `+` operator.
    ! *
    procedure(symmetric_operator),       pass(lhs), deferred, public :: integrand_multiply_integrand   !< `*` operator.
    procedure(integrand_op_real),        pass(lhs), deferred, public :: integrand_multiply_real        !< `* real` operator.
    procedure(real_op_integrand),        pass(rhs), deferred, public :: real_multiply_integrand        !< `real *` operator.
    procedure(integrand_op_real_scalar), pass(lhs), deferred, public :: integrand_multiply_real_scalar !< `* real_scalar` operator.
    procedure(real_scalar_op_integrand), pass(rhs), deferred, public :: real_scalar_multiply_integrand !< `real_scalar *` operator.
    generic, public :: operator(*) => integrand_multiply_integrand,   &
                                      integrand_multiply_real,        &
                                      real_multiply_integrand,        &
                                      integrand_multiply_real_scalar, &
                                      real_scalar_multiply_integrand !< Overloading `*` operator.
    ! -
    procedure(symmetric_operator), pass(lhs), deferred, public :: integrand_sub_integrand !< `-` operator.
    procedure(integrand_op_real),  pass(lhs), deferred, public :: integrand_sub_real      !< `- real` operator.
    procedure(real_op_integrand),  pass(rhs), deferred, public :: real_sub_integrand      !< `real -` operator.
    generic, public :: operator(-) => integrand_sub_integrand, &
                                      integrand_sub_real,      &
                                      real_sub_integrand !< Overloading `-` operator.
    ! =
    procedure(assignment_integrand), pass(lhs), deferred, public :: assign_integrand !< `=` operator.
    procedure(assignment_real),      pass(lhs), deferred, public :: assign_real      !< `= real` operator.
    generic, public :: assignment(=) => assign_integrand, assign_real !< Overloading `=` assignament.
endtype integrand_object

abstract interface
  !< Abstract type bound procedures necessary for implementing a concrete extension of [[integrand_object]].

  function time_derivative(self, t) result(dState_dt)
  !< Time derivative function of integrand class, i.e. the residuals function.
  import :: integrand_object, R_P
  class(integrand_object), intent(in) :: self         !< Integrand field.
  real(R_P),    optional,  intent(in) :: t            !< Time.
  real(R_P), allocatable              :: dState_dt(:) !< Result of the time derivative function of integrand field.
  endfunction time_derivative

  ! operators
  function local_error_operator(lhs, rhs) result(error)
  !< Estimate local truncation error between 2 solution approximations.
  import :: integrand_object, R_P
  class(integrand_object), intent(in) :: lhs   !< Left hand side.
  class(integrand_object), intent(in) :: rhs   !< Right hand side.
  real(R_P)                           :: error !< Error estimation.
  endfunction local_error_operator

  pure function integrand_op_real(lhs, rhs) result(operator_result)
  !< Asymmetric type operator `integrand.op.real`.
  import :: integrand_object, R_P
  class(integrand_object), intent(in) :: lhs                !< Left hand side.
  real(R_P),               intent(in) :: rhs(1:)            !< Right hand side.
  real(R_P), allocatable              :: operator_result(:) !< Operator result.
  endfunction integrand_op_real

  pure function real_op_integrand(lhs, rhs) result(operator_result)
  !< Asymmetric type operator `real.op.integrand`.
  import :: integrand_object, R_P
  class(integrand_object), intent(in) :: rhs                !< Right hand side.
  real(R_P),               intent(in) :: lhs(1:)            !< Left hand side.
  real(R_P), allocatable              :: operator_result(:) !< Operator result.
  endfunction real_op_integrand

  pure function integrand_op_real_scalar(lhs, rhs) result(operator_result)
  !< Asymmetric type operator `integrand.op.real`.
  import :: integrand_object, R_P
  class(integrand_object), intent(in) :: lhs                !< Left hand side.
  real(R_P),               intent(in) :: rhs                !< Right hand side.
  real(R_P), allocatable              :: operator_result(:) !< Operator result.
  endfunction integrand_op_real_scalar

  pure function real_scalar_op_integrand(lhs, rhs) result(operator_result)
  !< Asymmetric type operator `real.op.integrand`.
  import :: integrand_object, R_P
  real(R_P),               intent(in) :: lhs                !< Left hand side.
  class(integrand_object), intent(in) :: rhs                !< Right hand side.
  real(R_P), allocatable              :: operator_result(:) !< Operator result.
  endfunction real_scalar_op_integrand

  pure function symmetric_operator(lhs, rhs) result(operator_result)
  !< Symmetric type operator integrand.op.integrand.
  import :: integrand_object, R_P
  class(integrand_object), intent(in) :: lhs                !< Left hand side.
  class(integrand_object), intent(in) :: rhs                !< Right hand side.
  real(R_P), allocatable              :: operator_result(:) !< Operator result.
  endfunction symmetric_operator

  pure subroutine assignment_integrand(lhs, rhs)
  !< Symmetric assignment integrand = integrand.
  import :: integrand_object
  class(integrand_object), intent(inout) :: lhs !< Left hand side.
  class(integrand_object), intent(in)    :: rhs !< Right hand side.
  endsubroutine assignment_integrand

  pure subroutine assignment_real(lhs, rhs)
  !< Symmetric assignment integrand = integrand.
  import :: integrand_object, R_P
  class(integrand_object), intent(inout) :: lhs     !< Left hand side.
  real(R_P),               intent(in)    :: rhs(1:) !< Right hand side.
  endsubroutine assignment_real
endinterface
endmodule foodie_integrand_object
