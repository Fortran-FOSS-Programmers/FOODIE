!< FOODIE integrator: provide a predictor-corrector class of Adams-Bashforth-Moutlon multi-step schemes, from 1st to 4rd order
!< accurate.
module foodie_integrator_adams_bashforth_moulton
!-----------------------------------------------------------------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie_adt_integrand, only : integrand
use foodie_kinds, only : I_P, R_P
use foodie_integrator_adams_bashforth, only : adams_bashforth_integrator
use foodie_integrator_adams_moulton, only : adams_moulton_integrator
use foodie_utils, only : is_admissible
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: adams_bashforth_moulton_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
character(len=99), parameter :: supported_steps='1-16' !< List of supported steps number. Valid format is `1-2,4,9-23...`.
integer(I_P),      parameter :: min_ss=1               !< Minimum number of steps supported.
integer(I_P),      parameter :: max_ss=16              !< Maximum number of steps supported.

type :: adams_bashforth_moulton_integrator
  !< FOODIE integrator: provide an explicit class of Adams-Bashforth-Moulton multi-step schemes, from 1st to 4rd order accurate.
  !<
  !< @note The integrator must be created or initialized (predictor and corrector schemes selection) before used.
  !<
  !< ### List of errors status
  !<+ error=0 => no error;
  !<+ error=1 => bad (unsupported) number of required time steps;
  private
  integer(I_P)                     :: steps=-1  !< Number of time steps.
  type(adams_bashforth_integrator) :: predictor !< Predictor solver.
  type(adams_moulton_integrator)   :: corrector !< Corrector solver.
  integer(I_P)                     :: error=0   !< Error status flag: trap occurrences of errors.
  contains
    private
    procedure, pass(self), public :: init         !< Initialize (create) the integrator.
    procedure, pass(self), public :: destroy      !< Destroy the integrator.
    procedure, pass(self), public :: integrate    !< Integrate integrand field.
    procedure, nopass,     public :: min_steps    !< Return the minimum number of steps supported.
    procedure, nopass,     public :: max_steps    !< Return the maximum number of steps supported.
    procedure, nopass,     public :: is_supported !< Check if the queried number of steps is supported or not.
endtype adams_bashforth_moulton_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  elemental subroutine init(self, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual Adams-Bashforth-Moulton integrator: initialize the *b* coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_moulton_integrator), intent(INOUT) :: self  !< ABM integrator.
  integer(I_P), intent(IN)                                 :: steps !< Number of time steps used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (self%is_supported(steps)) then
    self%steps = steps
    call self%predictor%init(steps=steps)
    call self%corrector%init(steps=steps-1)
    self%error = 0
  else
    ! bad (unsupported) number of required time steps
    self%error = 1
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy the integrator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_moulton_integrator), intent(INOUT) :: self !< ABM integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = -1
  call self%predictor%destroy
  call self%corrector%destroy
  self%error = 0
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine integrate(self, U, previous, Dt, t, iterations)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with Adams-Bashforth-Moulton class scheme.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_moulton_integrator), intent(IN)    :: self         !< Actual ABM integrator.
  class(integrand),                          intent(INOUT) :: U            !< Field to be integrated.
  class(integrand),                          intent(INOUT) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                                 intent(IN)    :: Dt           !< Time steps.
  real(R_P),                                 intent(IN)    :: t(:)         !< Times.
  integer(I_P), optional,                    intent(IN)    :: iterations   !< Fixed point iterations of AM scheme.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%predictor%integrate(U=U, previous=previous, Dt=Dt, t=t, autoupdate=.false.)
  call self%corrector%integrate(U=U, previous=previous(2:), Dt=Dt, t=t, iterations=iterations, autoupdate=.false.)
  call self%predictor%update_previous(U=U, previous=previous)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate

  pure function min_steps()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the minimum number of steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P) :: min_steps !< Minimum number of steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  min_steps = min_ss
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction min_steps

  pure function max_steps()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the maximum number of steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P) :: max_steps !< Maximum number of steps supported.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  max_steps = max_ss
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction max_steps

  elemental function is_supported(steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Check if the queried number of steps is supported or not.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN) :: steps        !< Number of time steps used.
  logical                  :: is_supported !< Is true is the steps number is in *supported_steps*.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_supported = is_admissible(n=steps, adm_range=trim(supported_steps))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_supported
endmodule foodie_integrator_adams_bashforth_moulton
