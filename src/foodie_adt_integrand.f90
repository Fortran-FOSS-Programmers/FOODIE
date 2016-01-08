!< Define the abstract type *integrand* for building FOODIE ODE integrators.
module foodie_adt_integrand
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define the abstract type *integrand* for building FOODIE ODE integrators.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie_kinds, only : R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: integrand
  !< Abstract type for building FOODIE ODE integrators.
#ifdef CAF
  class(*), allocatable :: dummy_to_allow_extensions[:] !< Dummy member to allow concrete extensions with coarray members.
#endif
  contains
    ! public deferred procedures that concrete integrand-field must implement
    procedure(time_derivative),      pass(self), deferred, public :: t !< Time derivative, residuals.
    ! operators
    procedure(local_error_operator), pass(lhs),  deferred, public :: local_error                  !< ||integrand - integrand||.
    procedure(symmetric_operator),   pass(lhs),  deferred, public :: integrand_multiply_integrand !< Integrand * integrand operator.
    procedure(integrand_op_real),    pass(lhs),  deferred, public :: integrand_multiply_real      !< Integrand * real operator.
    procedure(real_op_integrand),    pass(rhs),  deferred, public :: real_multiply_integrand      !< Real * integrand operator.
    procedure(symmetric_operator),   pass(lhs),  deferred, public :: add                          !< Integrand + integrand operator.
    procedure(symmetric_operator),   pass(lhs),  deferred, public :: sub                          !< Integrand - integrand operator.
    procedure(assignment_integrand), pass(lhs),  deferred, public :: assign_integrand             !< Integrand = integrand.
    ! operators overloading
    generic, public :: operator(.lterror.) => local_error              !< Estimate local truncation error.
    generic, public :: operator(+) => add                              !< Overloading + operator.
    generic, public :: operator(-) => sub                              !< Overloading - operator.
    generic, public :: operator(*) => integrand_multiply_integrand, &
                                      real_multiply_integrand, &
                                      integrand_multiply_real          !< Overloading * operator.
    generic, public :: assignment(=) => assign_integrand               !< Overloading = assignament.
endtype integrand

abstract interface
  !< Abstract type bound procedures necessary for implementing a concrete extension of the class(integrand).
  function time_derivative(self, t) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative function of integrand class, i.e. the residuals function.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand, R_P
  class(integrand),       intent(IN) :: self      !< Integrand field.
  real(R_P),    optional, intent(IN) :: t         !< Time.
  class(integrand), allocatable      :: dState_dt !< Result of the time derivative function of integrand field.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction time_derivative

  function local_error_operator(lhs, rhs) result(error)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Estimate local truncation error between 2 solution approximations.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand, R_P
  class(integrand), intent(IN)  :: lhs   !< Left hand side.
  class(integrand), intent(IN)  :: rhs   !< Right hand side.
  real(R_P)                     :: error !< Error estimation.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction local_error_operator

  function integrand_op_real(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Asymmetric type operator integrand.op.real.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand, R_P
  class(integrand), intent(IN)  :: lhs              !< Left hand side.
  real(R_P),        intent(IN)  :: rhs              !< Right hand side.
  class(integrand), allocatable :: operator_result  !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction integrand_op_real

  function real_op_integrand(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Asymmetric type operator real.op.integrand.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand, R_P
  real(R_P),        intent(IN)  :: lhs             !< Left hand side.
  class(integrand), intent(IN)  :: rhs             !< Right hand side.
  class(integrand), allocatable :: operator_result !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_op_integrand

  function symmetric_operator(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Symmetric type operator integrand.op.integrand.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand
  class(integrand), intent(IN)  :: lhs             !< Left hand side.
  class(integrand), intent(IN)  :: rhs             !< Right hand side.
  class(integrand), allocatable :: operator_result !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction symmetric_operator

#ifdef PURE
  ! declare this pure should not be necessary, but this is a temporary workaround for a possible GNU gfortran bug.
  pure subroutine assignment_integrand(lhs, rhs)
#else
  subroutine assignment_integrand(lhs, rhs)
#endif
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Symmetric assignment integrand = integrand.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand
  class(integrand), intent(INOUT) :: lhs !< Left hand side.
  class(integrand), intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assignment_integrand
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule foodie_adt_integrand
