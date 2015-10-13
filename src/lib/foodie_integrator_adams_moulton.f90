!< FOODiE integrator: provide an implicit class of Adams-Moutlon multi-step schemes, from 1st to 4rd order accutate.
module foodie_integrator_adams_moulton
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODiE integrator: provide an implicit class of Adams-Moutlon multi-step schemes, from 1st to 4rd order accutate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the Adams-Moulton class scheme implemented is:
!<
!< $$ U^{n+N_s} = U^{n+N_s-1} +\Delta t \left[ \sum_{s=0}^{N_s-1}{ b_s \cdot R(t^{n+s}, U^{n+s}) } +
!< b_{N_S}\cdot R(t^{n+N_s}, U^{n+N_s}) \right] $$
!<
!<where \(N_s\) is the number of previous steps considered.
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are implicit. The coefficients \(b_s\) define the actual scheme, that is selected accordingly to the number of
!< **steps** used.
!<
!< Currently, the following schemes are available:
!<##### 0 step, Implicit Backward Euler, 1st order
!< This scheme is TVD and reverts to Implicit Backward Euler, it being 1st order.
!< The *b* coefficient is:
!< $$b = \left[b_0\right] = \left[1\right]$$
!< The scheme is:
!< $$ U^{n} = U^{n-1} + \Delta t R(t^{n},U^{n}) $$
!<
!<##### 1 step
!< This scheme is 2nd order.
!< The *b* coefficients are:
!< $$b = \left[ {\begin{array}{*{20}{c}} b_0 & b_1 \end{array}} \right] =
!<       \left[ {\begin{array}{*{20}{c}} \frac{1}{2} & \frac{1}{2} \end{array}} \right]$$
!< The scheme is:
!< $$ U^{n+1} = U^{n} +\Delta t \left[ \frac{1}{2} R(t^{n+1}, U^{n+1})+\frac{1}{2} R(t^{n}, U^{n}) \right] $$
!<
!<##### 2 steps
!< This scheme is 3rd order.
!< The *b* coefficients are:
!< $$b = \left[ {\begin{array}{*{20}{c}} b_0 & b_1 & b_2 \end{array}} \right] =
!<       \left[ {\begin{array}{*{20}{c}} -\frac{1}{12} & \frac{2}{3} & \frac{5}{12} \end{array}} \right]$$
!< The scheme is:
!< $$ U^{n+2} = U^{n+1} +\Delta t \left[ \frac{5}{12}R(t^{n+2}, U^{n+2}) + \frac{2}{3}R(t^{n+1}, U^{n+1})
!< -\frac{1}{12} R(t^{n}, U^{n})  \right] $$
!<
!<##### 3 steps
!< This scheme is 4th order.
!< The *b* coefficients are:
!< $$b = \left[ {\begin{array}{*{20}{c}} b_0 & b_1 & b_2 & b_3 \end{array}} \right] =
!<       \left[ {\begin{array}{*{20}{c}} \frac{1}{24} & -\frac{5}{24} & \frac{19}{24} & \frac{9}{24} \end{array}} \right]$$
!< The scheme is:
!< $$ U^{n+3} = U^{n+2} +\Delta t \left[ \frac{9}{24}R(t^{n+3}, U^{n+3}) + \frac{19}{24}R(t^{n+2}, U^{n+2})
!< -\frac{5}{24} R(t^{n+1}, U^{n+1}) + \frac{1}{24} R(t^{n}, U^{n}) \right] $$
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
public :: adams_moulton_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: adams_moulton_integrator
  !< FOODiE integrator: provide an explicit class of Adams-Moulton multi-step schemes, from 1st to 3rd order accutate.
  !<
  !< @note The integrator must be created or initialized (initialize the *b* coeficients) before used.
  private
  integer(I_P)           :: steps=-1 !< Number of time steps.
  real(R_P), allocatable :: b(:)     !< \(b\) coefficients.
  contains
    procedure, pass(self), public :: destroy         !< Destroy the integrator.
    procedure, pass(self), public :: init            !< Initialize (create) the integrator.
    procedure, pass(self), public :: integrate       !< Integrate integrand field.
    procedure, pass(self), public :: update_previous !< Cyclic update previous time steps.
    final                         :: finalize        !< Finalize object.
endtype adams_moulton_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  elemental subroutine init(self, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual Adams-Moulton integrator: initialize the *b* coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_moulton_integrator), intent(INOUT) :: self  !< AB integrator.
  integer(I_P), intent(IN)                       :: steps !< Number of time steps used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = steps
  if (allocated(self%b)) deallocate(self%b) ; allocate(self%b(0:steps)) ; self%b = 0._R_P
  select case(steps)
  case(0)
    ! AM(0) Bacward-Euler
    self%b(0) = 1._R_P
  case(1)
    ! AM(1)
    self%b(0) = 0.5_R_P
    self%b(1) = 0.5_R_P
  case(2)
    ! AM(2)
    self%b(0) = -1._R_P/12._R_P
    self%b(1) =  2._R_P/3._R_P
    self%b(2) =  5._R_P/12._R_P
  case(3)
    ! AM(3)
    self%b(0) =  1._R_P/24._R_P
    self%b(1) = -5._R_P/24._R_P
    self%b(2) = 19._R_P/24._R_P
    self%b(3) =  3._R_P/8._R_P
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy the integrator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_moulton_integrator), intent(INOUT) :: self !< AB integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = -1
  if (allocated(self%b)) deallocate(self%b)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine integrate(self, U, previous, Dt, t, autoupdate)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with Adams-Moulton class scheme.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_moulton_integrator), intent(IN)    :: self         !< Actual AB integrator.
  class(integrand),                intent(INOUT) :: U            !< Field to be integrated.
  class(integrand),                intent(INOUT) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                       intent(IN)    :: Dt           !< Time steps.
  real(R_P),                       intent(IN)    :: t(:)         !< Times.
  logical, optional,               intent(IN)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  logical                                        :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  integer(I_P)                                   :: s            !< Steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  if (self%steps>0) then
    U = previous(self%steps) + U%t(t=t(self%steps) + Dt) * (Dt * self%b(self%steps))
    do s=0, self%steps - 1
      U = U + previous(s+1)%t(t=t(s+1)) * (Dt * self%b(s))
    enddo
    if (autoupdate_) call self%update_previous(U=U, previous=previous)
  else
    U = U + U%t(t=t(1)) * (Dt * self%b(0))
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate

  subroutine update_previous(self, U, previous)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cyclic update previous time steps.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_moulton_integrator), intent(IN)    :: self         !< Actual AB integrator.
  class(integrand),                intent(IN)    :: U            !< Field to be integrated.
  class(integrand),                intent(INOUT) :: previous(1:) !< Previous time steps solutions of integrand field.
  integer(I_P)                                   :: s            !< Steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (self%steps>0) then
    do s=0, self%steps - 2
      previous(s + 1) = previous(s + 2)
    enddo
    previous(self%steps) = U
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_previous

  ! private methods
  elemental subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Finalize object.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(adams_moulton_integrator), intent(INOUT) :: self !< AB integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule foodie_integrator_adams_moulton
