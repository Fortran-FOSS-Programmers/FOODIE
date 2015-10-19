!< FOODIE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 4rd order accutate.
module foodie_integrator_adams_bashforth
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 4rd order accutate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the Adams-Bashforth class scheme implemented is:
!<
!< $$ U^{n+N_s} = U^{n+N_s-1} +\Delta t \left[ \sum_{s=1}^{N_s}{ b_s \cdot R(t^{n+s-1}, U^{n+s-1}) } \right] $$
!<
!<where \(N_s\) is the number of previous steps considered.
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit. The coefficients \(b_s\) define the actual scheme, that is selected accordingly to the number of
!< **steps** used.
!<
!< Currently, the following schemes are available:
!<##### 1 step, Explicit Forward Euler, 1st order
!< This scheme is TVD and reverts to Explicit Forward Euler, it being 1st order.
!< The *b* coefficient is:
!< $$b = \left[b_1\right] = \left[1\right]$$
!< The scheme is:
!< $$ U^{n+1} = U^n + \Delta t R(t^n,U^n) $$
!<
!<##### 2 steps
!< This scheme is 2nd order.
!< The *b* coefficients are:
!< $$b = \left[ {\begin{array}{*{20}{c}} b_1 & b_2 \end{array}} \right] =
!<       \left[ {\begin{array}{*{20}{c}} -\frac{1}{2} & \frac{3}{2} \end{array}} \right]$$
!< The scheme is:
!< $$ U^{n+2} = U^{n+1} +\Delta t \left[ \frac{3}{2} R(t^{n+1}, U^{n+1})-\frac{1}{2} R(t^{n}, U^{n})  \right] $$
!<
!<##### 3 steps
!< This scheme is 3rd order.
!< The *b* coefficients are:
!< $$b = \left[ {\begin{array}{*{20}{c}} b_1 & b_2 & b_3 \end{array}} \right] =
!<       \left[ {\begin{array}{*{20}{c}} \frac{5}{12} & -\frac{4}{3} & \frac{23}{12} \end{array}} \right]$$
!< The scheme is:
!< $$ U^{n+3} = U^{n+2} +\Delta t \left[ \frac{23}{12}R(t^{n+2}, U^{n+2}) - \frac{4}{3}R(t^{n+1}, U^{n+1})
!< +\frac{5}{12} R(t^{n}, U^{n})  \right] $$
!<
!<##### 4 steps
!< This scheme is 4th order.
!< The *b* coefficients are:
!< $$b = \left[ {\begin{array}{*{20}{c}} b_1 & b_2 & b_3 & b_4 \end{array}} \right] =
!<       \left[ {\begin{array}{*{20}{c}} -\frac{9}{24} & \frac{37}{24} & -\frac{59}{24} & \frac{55}{24} \end{array}} \right]$$
!< The scheme is:
!< $$ U^{n+4} = U^{n+3} +\Delta t \left[ \frac{55}{24}R(t^{n+3}, U^{n+3}) - \frac{59}{24}R(t^{n+2}, U^{n+2})
!< +\frac{37}{24} R(t^{n+1}, U^{n+1}) - \frac{9}{24} R(t^{n}, U^{n}) \right] $$
!<
!<#### Bibliography
!<
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie_kinds, only : R_P, I_P
use foodie_adt_integrand, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: adams_bashforth_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: adams_bashforth_integrator
  !< FOODIE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 4rd order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *b* coefficients) before used.
  private
  integer(I_P)           :: steps=0 !< Number of time steps.
  real(R_P), allocatable :: b(:)    !< \(b\) coefficients.
  contains
    procedure, pass(self), public :: destroy         !< Destroy the integrator.
    procedure, pass(self), public :: init            !< Initialize (create) the integrator.
    procedure, pass(self), public :: integrate       !< Integrate integrand field.
    procedure, pass(self), public :: update_previous !< Cyclic update previous time steps.
    final                         :: finalize        !< Finalize object.
endtype adams_bashforth_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  elemental subroutine init(self, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual Adams-Bashforth integrator: initialize the *b* coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_integrator), intent(INOUT) :: self  !< AB integrator.
  integer(I_P), intent(IN)                         :: steps !< Number of time steps used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = steps
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
  case(4)
    ! AB(4)
    self%b(1) =  -3._R_P/8._R_P
    self%b(2) =  37._R_P/24._R_P
    self%b(3) = -59._R_P/24._R_P
    self%b(4) =  55._R_P/24._R_P
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
  self%steps = 0
  if (allocated(self%b)) deallocate(self%b)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine integrate(self, U, previous, Dt, t, autoupdate)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with Adams-Bashforth class scheme.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_integrator), intent(IN)    :: self         !< Actual AB integrator.
  class(integrand),                  intent(INOUT) :: U            !< Field to be integrated.
  class(integrand),                  intent(INOUT) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                         intent(IN)    :: Dt           !< Time steps.
  real(R_P),                         intent(IN)    :: t(:)         !< Times.
  logical, optional,                 intent(IN)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  logical                                          :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  integer(I_P)                                     :: s            !< Steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  do s=1, self%steps
    U = U + previous(s)%t(t=t(s)) * (Dt * self%b(s))
  enddo
  if (autoupdate_) call self%update_previous(U=U, previous=previous)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate

  subroutine update_previous(self, U, previous)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cyclic update previous time steps.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_integrator), intent(IN)    :: self         !< Actual AB integrator.
  class(integrand),                  intent(IN)    :: U            !< Field to be integrated.
  class(integrand),                  intent(INOUT) :: previous(1:) !< Previous time steps solutions of integrand field.
  integer(I_P)                                     :: s            !< Steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do s=1, self%steps - 1
    previous(s) = previous(s + 1)
  enddo
  previous(self%steps) = U
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_previous

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
