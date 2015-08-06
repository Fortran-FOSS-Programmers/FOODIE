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
    procedure(time_derivative     ), pass(self), deferred :: t
    procedure(asymmetric_operator ), pass(lhs),  deferred :: multiply
    procedure(symmetric_operator  ), pass(lhs),  deferred :: add
    procedure(symmetric_assignment), pass(lhs),  deferred :: assign
    generic :: operator(+) => add
    generic :: operator(*) => multiply
    generic :: assignment(=) => assign
endtype integrand

abstract interface
  function time_derivative(self) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative function of integrand class, i.e. the residuals function.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand
  class(integrand), intent(in)  :: self
  class(integrand), allocatable :: dState_dt
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction time_derivative

  function asymmetric_operator(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Asymmetric type operator integrand.op.real.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand, R_P
  class(integrand), intent(in)  :: lhs
  real(R_P),        intent(in)  :: rhs
  class(integrand), allocatable :: operator_result
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction asymmetric_operator

  function symmetric_operator(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Symmetric type operator integrand.op.integrand.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand
  class(integrand), intent(in)  :: lhs, rhs
  class(integrand), allocatable :: operator_result
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction symmetric_operator

  subroutine symmetric_assignment(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Symmetric assignamente integrand = integrand.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: integrand
  class(integrand) ,intent(inout) :: lhs
  class(integrand) ,intent(in)    :: rhs
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine symmetric_assignment
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_integrand
