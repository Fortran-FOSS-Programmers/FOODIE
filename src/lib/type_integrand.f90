!< Define the abstract type *integrand* for building FOODiE ODE integrators.
module type_integrand
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define the abstract type *integrand* for building FOODiE ODE integrators.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: integrand
  !< Abstract type for building FOODiE ODE integrators.
  contains
    ! public deferred procedures that concrete interpolators must implement
    procedure(time_derivative),      pass(self), deferred, public :: t                            !< Time derivative, residuals.
    procedure(symmetric_operator),   pass(lhs),  deferred, public :: integrand_multiply_integrand !< Integrand * integrand operator.
    procedure(integrand_op_real),    pass(lhs),  deferred, public :: integrand_multiply_real      !< Integrand * real operator.
    procedure(real_op_integrand),    pass(rhs),  deferred, public :: real_multiply_integrand      !< Real * integrand operator.
    procedure(symmetric_operator),   pass(lhs),  deferred, public :: add                          !< Integrand + integrand oprator.
    procedure(assignment_integrand), pass(lhs),  deferred, public :: assign_integrand             !< Integrand = integrand.
    procedure(assignment_real),      pass(lhs),  deferred, public :: assign_real                  !< Integrand = real.
    ! operators overloading
    generic, public :: operator(+) => add                              !< Overloading + operator.
    generic, public :: operator(*) => integrand_multiply_integrand, &
                                      real_multiply_integrand, &
                                      integrand_multiply_real          !< Overloading * operator.
    generic, public :: assignment(=) => assign_integrand, assign_real  !< Overloading = assignament.
endtype integrand

abstract interface
  !< Abstract type bound procedures necessary for implementing the class(integrand).
  pure function time_derivative(self, n) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative function of integrand class, i.e. the residuals function.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand
  class(integrand),  intent(IN) :: self      !< Integrand field.
  integer, optional, intent(IN) :: n         !< Time level.
  class(integrand), allocatable :: dState_dt !< Result of the time derivative function of integrand field.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction time_derivative

  pure function integrand_op_real(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Asymmetric type operator integrand.op.real.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand, R_P
  class(integrand), intent(IN)  :: lhs              !< Left hand side.
  real(R_P),        intent(IN)  :: rhs              !< Right hand side.
  class(integrand), allocatable :: operator_result  !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction integrand_op_real

  pure function real_op_integrand(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Asymmetric type operator real.op.integrand.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand, R_P
  real(R_P),        intent(IN)  :: lhs             !< Left hand side.
  class(integrand), intent(IN)  :: rhs             !< Right hand side.
  class(integrand), allocatable :: operator_result !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_op_integrand

  pure function symmetric_operator(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Symmetric type operator integrand.op.integrand.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand
  class(integrand), intent(IN)  :: lhs             !< Left hand side.
  class(integrand), intent(IN)  :: rhs             !< Right hand side.
  class(integrand), allocatable :: operator_result !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction symmetric_operator

  pure subroutine assignment_integrand(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Symmetric assignment integrand = integrand.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand
  class(integrand), intent(INOUT) :: lhs !< Left hand side.
  class(integrand), intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assignment_integrand

  pure subroutine assignment_real(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Asymmetric assignment integrand = real.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand, R_P
  class(integrand), intent(INOUT) :: lhs !< Left hand side.
  real(R_P),        intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assignment_real
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_integrand
