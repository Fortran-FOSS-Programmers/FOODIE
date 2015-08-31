!< FOODiE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 3rd order accutate.
module foodie_integrator_adams_bashforth
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODiE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 3rd order accutate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the Adams-Bashforth class scheme implemented is:
!<
!< $$ U^{n+s} = U^{n+s-1} +\Delta t \left[ \sum_{i=1}^{n+s}{ b_i \cdot R(t^{n+i-1}, U^{n+i-1}) } \right] $$
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit. The coefficients \(\b\) define the actual scheme, that is selected accordingly to the number of
!< **steps** used.
!<
!< Currently the following schemes are available:
!<##### 1 step, Explicit Forward Euler, 1st order
!< This scheme is TVD and reverts to Explicit Forward Euler, it being 1st order.
!< $$b = \left[1\right]$$
!<
!<##### 2 steps
!< This scheme is 2nd order.
!< $$b = \left[ {\begin{array}{*{20}{c}} b_1 & b_2 \end{array}} \right] =
!<       \left[ {\begin{array}{*{20}{c}} -\frac{1}{2} & \frac{3}{2} \end{array}} \right]$$
!<
!<##### 3 steps
!< This scheme is 3rd order.
!< $$b = \left[ {\begin{array}{*{20}{c}} b_1 & b_2 & b_3 \end{array}} \right] =
!<       \left[ {\begin{array}{*{20}{c}} \frac{5}{12} & -\frac{4}{3} & \frac{23}{12} \end{array}} \right]$$
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P
use type_integrand, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: adams_bashforth_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: adams_bashforth_integrator
  !< FOODiE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 3rd order accutate.
  !<
  !< @note The integrator must be created or initialized (initialize the *b* coeficients) before used.
  real(R_P), allocatable:: b(:)   !< \(b\) coefficients.
  contains
    procedure, pass(self), public :: destroy   !< Destroy the integrator.
    procedure, pass(self), public :: init      !< Initialize (create) the integrator.
    procedure, pass(self), public :: integrate !< Integrate integrand field.
    final                         :: finalize  !< Finalize object.
endtype adams_bashforth_integrator
interface adams_bashforth_integrator
  !< Overload adams_bashforth_integrator name with the constructor function *create*.
  module procedure create
endinterface adams_bashforth_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function create(steps) result(ab)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual Adams-Bashforth integrator: initialize the *b* coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)         :: steps !< Number of time steps used.
  type(adams_bashforth_integrator) :: ab    !< Actual AB integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call ab%init(steps=steps)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction create

  ! public methods
  elemental subroutine init(self, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual Adams-Bashforth integrator: initialize the *b* coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_integrator), intent(INOUT) :: self  !< AB integrator.
  integer(I_P), intent(IN)                         :: steps !< Number of time steps used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%b)) deallocate(self%b) ; allocate(self%b(1:steps)) ; self%b = 0._R_P
  select case(steps)
  case(1)
    ! AB(1) Forward-Euler
    self%b(1) = 1._R_P
  case(2)
    ! AB(2)
    self%b(1) = -0.5_R_P
    self%b(2) = 1.5_R_P
  case(3)
    ! AB(3)
    self%b(1) =  5._R_P/12._R_P
    self%b(2) = -4._R_P/3._R_P
    self%b(3) = 23._R_P/12._R_P
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy the integrator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_integrator), intent(INOUT) :: self !< AB integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%b)) deallocate(self%b)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine integrate(self, field, steps, dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with Adams-Bashforth class scheme.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_integrator), intent(IN)    :: self  !< Actual AB integrator.
  class(integrand),                  intent(INOUT) :: field !< Field to be integrated.
  integer(I_P),                      intent(IN)    :: steps !< Number of time steps used.
  real(R_P),                         intent(in)    :: dt    !< Time step.
  integer(I_P)                                     :: s     !< Steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do s=1, steps
    field = field +  field%t(n=s) * (dt * self%b(s))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate

  ! private methods
  elemental subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Finalize object.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(adams_bashforth_integrator), intent(INOUT) :: self !< AB integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule foodie_integrator_adams_bashforth
