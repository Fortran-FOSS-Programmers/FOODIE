!< FOODIE integrator: provide a predictor-corrector class of Adams-Bashforth-Moutlon multi-step schemes, from 1st to 4rd order
!< accurate.

module foodie_integrator_adams_bashforth_moulton
!< FOODIE integrator: provide a predictor-corrector class of Adams-Bashforth-Moutlon multi-step schemes, from 1st to 4rd order
!< accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the Adams-Bashforth-Moulton class scheme implemented is:
!<
!<##### predictor
!< $$ U^{n+N_s^p} = U^{n+N_s^p-1} +\Delta t \left[ \sum_{s=1}^{N_s^p}{ b_s^p \cdot R(t^{n+s-1}, U^{n+s-1}) } \right] $$
!<##### corrector
!< $$ U^{n+N_s^c} = U^{n+N_s^c-1} +\Delta t \left[ \sum_{s=0}^{N_s^c}{ b_s^c \cdot R(t^{n+s}, U^{n+s}) } \right] $$
!<
!<where \(N_s^p\) is the number of previous steps considered for the predictor and \(N_s^c\) is the number of previous steps
!<considered for the corrector.
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The coefficients \(b_s^{p,c}\) define the actual scheme, that is selected accordingly to the number of
!< **steps** used. The predictor and corrector schemes should have the same formal order of accuracy, thus \(N_s^p=N_s^c+1\)
!< should hold.
!<
!< Currently, the following schemes are available:
!<##### P=AB(1)-C=AM(0) step, Explicit/Implicit Farward/Backward Euler, 1st order
!< This scheme is TVD and reverts to Explicit/Implicit Farward/Backward Euler, it being 1st order.
!< The *b* coefficient is:
!< $$b^p = \left[b_1\right] = \left[1\right]$$
!< $$b^c = \left[b_0\right] = \left[1\right]$$
!< The scheme is:
!< $$ U^{n+1,p} = U^{n} + \Delta t R(t^{n},U^{n}) $$
!< $$ U^{n+1,c} = U^{n} + \Delta t R(t^{n+1},U^{n+1,p}) $$
!<
!<##### P=AB(2)-C=AM(1) steps
!< This scheme is 2nd order.
!< The *b* coefficients are:
!< $$b^p = \left[ {\begin{array}{*{20}{c}} b_1 & b_2 \end{array}} \right] =
!<         \left[ {\begin{array}{*{20}{c}} -\frac{1}{2} & \frac{3}{2} \end{array}} \right]$$
!< $$b^c = \left[ {\begin{array}{*{20}{c}} b_0 & b_1 \end{array}} \right] =
!<         \left[ {\begin{array}{*{20}{c}} \frac{1}{2} & \frac{1}{2} \end{array}} \right]$$
!< The scheme is:
!< $$ U^{n+2,p} = U^{n+1} +\Delta t \left[ \frac{3}{2} R(t^{n+1}, U^{n+1})-\frac{1}{2} R(t^{n}, U^{n})  \right] $$
!< $$ U^{n+2,c} = U^{n+1} +\Delta t \left[ \frac{1}{2} R(t^{n+2,p}, U^{n+2})+\frac{1}{2} R(t^{n+1}, U^{n+1}) \right] $$
!<
!<##### P=AB(3)-C=AM(2) steps
!< This scheme is 3rd order.
!< The *b* coefficients are:
!< $$b^p = \left[ {\begin{array}{*{20}{c}} b_1 & b_2 & b_3 \end{array}} \right] =
!<         \left[ {\begin{array}{*{20}{c}} \frac{5}{12} & -\frac{4}{3} & \frac{23}{12} \end{array}} \right]$$
!< $$b^c = \left[ {\begin{array}{*{20}{c}} b_0 & b_1 & b_2 \end{array}} \right] =
!<         \left[ {\begin{array}{*{20}{c}} -\frac{1}{12} & \frac{2}{3} & \frac{5}{12} \end{array}} \right]$$
!< The scheme is:
!< $$ U^{n+3,p} = U^{n+2} +\Delta t \left[ \frac{23}{12}R(t^{n+2}, U^{n+2}) - \frac{4}{3}R(t^{n+1}, U^{n+1})
!< +\frac{5}{12} R(t^{n}, U^{n})  \right] $$
!< $$ U^{n+3,c} = U^{n+2} +\Delta t \left[ \frac{5}{12}R(t^{n+3}, U^{n+3,p}) + \frac{2}{3}R(t^{n+2}, U^{n+2})
!< -\frac{1}{12} R(t^{n+1}, U^{n+1})  \right] $$
!<
!<##### P=AB(4)-C=AM(3) steps
!< This scheme is 4th order.
!< The *b* coefficients are:
!< $$b^p = \left[ {\begin{array}{*{20}{c}} b_1 & b_2 & b_3 & b_4 \end{array}} \right] =
!<         \left[ {\begin{array}{*{20}{c}} -\frac{9}{24} & \frac{37}{24} & -\frac{59}{24} & \frac{55}{24} \end{array}} \right]$$
!< $$b^c = \left[ {\begin{array}{*{20}{c}} b_0 & b_1 & b_2 & b_3 \end{array}} \right] =
!<         \left[ {\begin{array}{*{20}{c}} \frac{1}{24} & -\frac{5}{24} & \frac{19}{24} & \frac{9}{24} \end{array}} \right]$$
!< The scheme is:
!< $$ U^{n+4,p} = U^{n+3} +\Delta t \left[ \frac{55}{24}R(t^{n+3}, U^{n+3}) - \frac{59}{24}R(t^{n+2}, U^{n+2})
!< +\frac{37}{24} R(t^{n+1}, U^{n+1}) - \frac{9}{24} R(t^{n}, U^{n}) \right] $$
!< $$ U^{n+4,c} = U^{n+3} +\Delta t \left[ \frac{9}{24}R(t^{n+4}, U^{n+4,p}) + \frac{19}{24}R(t^{n+3}, U^{n+3})
!< -\frac{5}{24} R(t^{n+2}, U^{n+2}) + \frac{1}{24} R(t^{n+1}, U^{n+1}) \right] $$
!<
!<#### Bibliography
!<

use foodie_adt_integrand, only : integrand
use foodie_error_codes, only : ERROR_BAD_STEPS_NUMBER, ERROR_MISSING_STEPS_NUMBER
use foodie_kinds, only : I_P, R_P
use foodie_integrator_adams_bashforth, only : integrator_adams_bashforth
use foodie_integrator_adams_moulton, only : integrator_adams_moulton
use foodie_integrator_object, only : integrator_object
use foodie_utils, only : is_admissible

implicit none
private
public :: integrator_adams_bashforth_moulton

character(len=99), parameter :: supported_steps='1-16' !< List of supported steps number. Valid format is `1-2,4,9-23...`.
integer(I_P),      parameter :: min_ss=1               !< Minimum number of steps supported.
integer(I_P),      parameter :: max_ss=16              !< Maximum number of steps supported.

type, extends(integrator_object) :: integrator_adams_bashforth_moulton
  !< FOODIE integrator: provide an explicit class of Adams-Bashforth-Moulton multi-step schemes, from 1st to 4rd order accurate.
  !<
  !< @note The integrator must be created or initialized (predictor and corrector schemes selection) before used.
  private
  integer(I_P)                     :: steps=-1  !< Number of time steps.
  type(integrator_adams_bashforth) :: predictor !< Predictor solver.
  type(integrator_adams_moulton)   :: corrector !< Corrector solver.
  contains
    ! deferred methods
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    ! public methods
    procedure, pass(self) :: destroy      !< Destroy the integrator.
    procedure, pass(self) :: init         !< Initialize (create) the integrator.
    procedure, pass(self) :: integrate    !< Integrate integrand field.
    procedure, nopass     :: is_supported !< Check if the queried number of steps is supported or not.
    procedure, nopass     :: min_steps    !< Return the minimum number of steps supported.
    procedure, nopass     :: max_steps    !< Return the maximum number of steps supported.
endtype integrator_adams_bashforth_moulton

contains
  ! deferred methods
  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_adams_bashforth_moulton), intent(in)           :: self             !< Integrator.
  character(*),                              intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                                   :: desc             !< Description.
  character(len=:), allocatable                                   :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                                     :: NL=new_line('a') !< New line character.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = ''
  desc = desc//prefix_//'Adams-Bashforth-Moulton multi-step schemes class'//NL
  desc = desc//prefix_//'  Supported steps numbers: ['//trim(adjustl(supported_steps))//']'
  endfunction description

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_adams_bashforth_moulton), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),                  intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  select type(rhs)
  class is (integrator_adams_bashforth_moulton)
    lhs%steps     = rhs%steps
    lhs%predictor = rhs%predictor
    lhs%corrector = rhs%corrector
  endselect
  endsubroutine integr_assign_integr

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_adams_bashforth_moulton), intent(inout) :: self !< Integrator.

  call self%destroy_abstract
  self%steps = -1
  call self%predictor%destroy
  call self%corrector%destroy
  endsubroutine destroy

  subroutine init(self, steps)
  !< Create the actual Adams-Bashforth-Moulton integrator: initialize the *b* coefficients.
  class(integrator_adams_bashforth_moulton), intent(inout) :: self  !< Integrator.
  integer(I_P),                              intent(in)    :: steps !< Number of time steps used.

  if (self%is_supported(steps)) then
    call self%destroy
    self%steps = steps
    call self%predictor%init(steps=steps)
    call self%corrector%init(steps=steps-1)
  else
    call self%trigger_error(error=ERROR_BAD_STEPS_NUMBER,                           &
                            error_message='bad (unsupported) number of time steps', &
                            is_severe=.true.)
  endif
  endsubroutine init

  subroutine integrate(self, U, previous, Dt, t, iterations)
  !< Integrate field with Adams-Bashforth-Moulton class scheme.
  class(integrator_adams_bashforth_moulton), intent(in)           :: self         !< Integrator.
  class(integrand),                          intent(inout)        :: U            !< Field to be integrated.
  class(integrand),                          intent(inout)        :: previous(1:) !< Previous time steps solutions of integrand.
  real(R_P),                                 intent(in)           :: Dt           !< Time steps.
  real(R_P),                                 intent(in)           :: t(:)         !< Times.
  integer(I_P),                              intent(in), optional :: iterations   !< Fixed point iterations of AM scheme.

  call self%predictor%integrate(U=U, previous=previous, Dt=Dt, t=t, autoupdate=.false.)
  call self%corrector%integrate(U=U, previous=previous(2:), Dt=Dt, t=t, iterations=iterations, autoupdate=.false.)
  call self%predictor%update_previous(U=U, previous=previous)
  endsubroutine integrate

  elemental function is_supported(steps)
  !< Check if the queried number of steps is supported or not.
  integer(I_P), intent(in) :: steps        !< Number of time steps used.
  logical                  :: is_supported !< Is true is the steps number is in *supported_steps*.

  is_supported = is_admissible(n=steps, adm_range=trim(supported_steps))
  endfunction is_supported

  pure function min_steps()
  !< Return the minimum number of steps supported.
  integer(I_P) :: min_steps !< Minimum number of steps supported.

  min_steps = min_ss
  endfunction min_steps

  pure function max_steps()
  !< Return the maximum number of steps supported.
  integer(I_P) :: max_steps !< Maximum number of steps supported.

  max_steps = max_ss
  endfunction max_steps
endmodule foodie_integrator_adams_bashforth_moulton
