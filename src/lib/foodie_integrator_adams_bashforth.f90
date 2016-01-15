!< FOODIE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 16th order accurate.
module foodie_integrator_adams_bashforth
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 16th order accurate.
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
!< The schemes are explicit. The coefficients *b* define the actual scheme, that is selected accordingly to the number of
!< **steps** used. Currently, the schemes provided have steps number in *[1, 16]*. Note that the scheme using only 1 step reverts
!< to Explicit Forward Euler. The formal order of accuracy varies consistently in *[1st, 16th]* order.
!<
!<#### Bibliography
!< [1] *Cowell Type Numerical Integration As Applied to Satellite Orbit Computation*, J. L. Maury Jr.,
!< G. P. Segal, X-553-69-46, April 1969, [NASA-TM-X-63542](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19690017325.pdf).
!<
!< [2] *Linear multistep method*, [wikipedia article](https://en.wikipedia.org/wiki/Linear_multistep_method).
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie_adt_integrand, only : integrand
use foodie_kinds, only : I_P, R_P
use foodie_utils, only : is_admissible
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: adams_bashforth_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
character(len=99), parameter :: supported_steps='1-16' !< List of supported steps number. Valid format is `1-2,4,9-23...`.
integer(I_P),      parameter :: min_ss=1               !< Minimum number of steps supported.
integer(I_P),      parameter :: max_ss=16              !< Maximum number of steps supported.

type :: adams_bashforth_integrator
  !< FOODIE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 16th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *b* coefficients) before used.
  !<
  !< ### List of errors status
  !<+ error=0 => no error;
  !<+ error=1 => bad (unsupported) number of required time steps;
  private
  integer(I_P)                :: steps=0  !< Number of time steps.
  real(R_P), allocatable      :: b(:)     !< *b* coefficients.
  integer(I_P)                :: error=0  !< Error status flag: trap occurrences of errors.
  contains
    private
    procedure, pass(self), public :: init            !< Initialize (create) the integrator.
    procedure, pass(self), public :: destroy         !< Destroy the integrator.
    procedure, pass(self), public :: integrate       !< Integrate integrand field.
    procedure, pass(self), public :: update_previous !< Cyclic update previous time steps.
    procedure, nopass,     public :: min_steps       !< Return the minimum number of steps supported.
    procedure, nopass,     public :: max_steps       !< Return the maximum number of steps supported.
    procedure, nopass,     public :: is_supported    !< Check if the queried number of steps is supported or not.
    final                         :: finalize        !< Finalize object.
endtype adams_bashforth_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  elemental subroutine init(self, steps)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual Adams-Bashforth integrator: initialize the *b* coefficients.
  !<
  !< @note If the integrator is initialized with a bad (unsupported) number of required time steps the initialization fails and
  !< the integrator error status is updated consistently for external-provided errors handling.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_bashforth_integrator), intent(INOUT) :: self  !< AB integrator.
  integer(I_P),                      intent(IN)    :: steps !< Number of time steps used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (self%is_supported(steps)) then
    self%steps = steps
    if (allocated(self%b)) deallocate(self%b) ; allocate(self%b(1:steps)) ; self%b = 0.0_R_P
    select case(steps)
    case(1)
      ! AB(1) Forward-Euler
      self%b(1) = 1.0_R_P
    case(2)
      self%b(1) = -1.0_R_P/2.0_R_P
      self%b(2) = 3.0_R_P/2.0_R_P
    case(3)
      self%b(1) = 5.0_R_P/12.0_R_P
      self%b(2) = -16.0_R_P/12.0_R_P
      self%b(3) = 23.0_R_P/12.0_R_P
    case(4)
      self%b(1) = -9.0_R_P/24.0_R_P
      self%b(2) = 37.0_R_P/24.0_R_P
      self%b(3) = -59.0_R_P/24.0_R_P
      self%b(4) = 55.0_R_P/24.0_R_P
    case(5)
      self%b(1) = 251.0_R_P/720.0_R_P
      self%b(2) = -1274.0_R_P/720.0_R_P
      self%b(3) = 2616.0_R_P/720.0_R_P
      self%b(4) = -2774.0_R_P/720.0_R_P
      self%b(5) = 1901.0_R_P/720.0_R_P
    case(6)
      self%b(1) = -475.0_R_P/1440.0_R_P
      self%b(2) = 2877.0_R_P/1440.0_R_P
      self%b(3) = -7298.0_R_P/1440.0_R_P
      self%b(4) = 9982.0_R_P/1440.0_R_P
      self%b(5) = -7923.0_R_P/1440.0_R_P
      self%b(6) = 4277.0_R_P/1440.0_R_P
    case(7)
      self%b(1) = 19087.0_R_P/60480.0_R_P
      self%b(2) = -134472.0_R_P/60480.0_R_P
      self%b(3) = 407139.0_R_P/60480.0_R_P
      self%b(4) = -688256.0_R_P/60480.0_R_P
      self%b(5) = 705549.0_R_P/60480.0_R_P
      self%b(6) = -447288.0_R_P/60480.0_R_P
      self%b(7) = 198721.0_R_P/60480.0_R_P
    case(8)
      self%b(1) = -36799.0_R_P/120960.0_R_P
      self%b(2) = 295767.0_R_P/120960.0_R_P
      self%b(3) = -1041723.0_R_P/120960.0_R_P
      self%b(4) = 2102243.0_R_P/120960.0_R_P
      self%b(5) = -2664477.0_R_P/120960.0_R_P
      self%b(6) = 2183877.0_R_P/120960.0_R_P
      self%b(7) = -1152169.0_R_P/120960.0_R_P
      self%b(8) = 434241.0_R_P/120960.0_R_P
    case(9)
      self%b(1) = 1070017.0_R_P/3628800.0_R_P
      self%b(2) = -9664106.0_R_P/3628800.0_R_P
      self%b(3) = 38833486.0_R_P/3628800.0_R_P
      self%b(4) = -91172642.0_R_P/3628800.0_R_P
      self%b(5) = 137968480.0_R_P/3628800.0_R_P
      self%b(6) = -139855262.0_R_P/3628800.0_R_P
      self%b(7) = 95476786.0_R_P/3628800.0_R_P
      self%b(8) = -43125206.0_R_P/3628800.0_R_P
      self%b(9) = 14097247.0_R_P/3628800.0_R_P
    case(10)
      self%b(1) = -2082753.0_R_P/7257600.0_R_P
      self%b(2) = 20884811.0_R_P/7257600.0_R_P
      self%b(3) = -94307320.0_R_P/7257600.0_R_P
      self%b(4) = 252618224.0_R_P/7257600.0_R_P
      self%b(5) = -444772162.0_R_P/7257600.0_R_P
      self%b(6) = 538363838.0_R_P/7257600.0_R_P
      self%b(7) = -454661776.0_R_P/7257600.0_R_P
      self%b(8) = 265932680.0_R_P/7257600.0_R_P
      self%b(9) = -104995189.0_R_P/7257600.0_R_P
      self%b(10) = 30277247.0_R_P/7257600.0_R_P
    case(11)
      self%b(1) = 134211265.0_R_P/479001600.0_R_P
      self%b(2) = -1479574348.0_R_P/479001600.0_R_P
      self%b(3) = 7417904451.0_R_P/479001600.0_R_P
      self%b(4) = -22329634920.0_R_P/479001600.0_R_P
      self%b(5) = 44857168434.0_R_P/479001600.0_R_P
      self%b(6) = -63176201472.0_R_P/479001600.0_R_P
      self%b(7) = 63716378958.0_R_P/479001600.0_R_P
      self%b(8) = -46113029016.0_R_P/479001600.0_R_P
      self%b(9) = 23591063805.0_R_P/479001600.0_R_P
      self%b(10) = -8271795124.0_R_P/479001600.0_R_P
      self%b(11) = 2132509567.0_R_P/479001600.0_R_P
    case(12)
      self%b(1) = -262747265.0_R_P/958003200.0_R_P
      self%b(2) = 3158642445.0_R_P/958003200.0_R_P
      self%b(3) = -17410248271.0_R_P/958003200.0_R_P
      self%b(4) = 58189107627.0_R_P/958003200.0_R_P
      self%b(5) = -131365867290.0_R_P/958003200.0_R_P
      self%b(6) = 211103573298.0_R_P/958003200.0_R_P
      self%b(7) = -247741639374.0_R_P/958003200.0_R_P
      self%b(8) = 214139355366.0_R_P/958003200.0_R_P
      self%b(9) = -135579356757.0_R_P/958003200.0_R_P
      self%b(10) = 61633227185.0_R_P/958003200.0_R_P
      self%b(11) = -19433810163.0_R_P/958003200.0_R_P
      self%b(12) = 4527766399.0_R_P/958003200.0_R_P
    case(13)
      self%b(1) = 703604254357.0_R_P/2615348736000.0_R_P
      self%b(2) = -9160551085734.0_R_P/2615348736000.0_R_P
      self%b(3) = 55060974662412.0_R_P/2615348736000.0_R_P
      self%b(4) = -202322913738370.0_R_P/2615348736000.0_R_P
      self%b(5) = 507140369728425.0_R_P/2615348736000.0_R_P
      self%b(6) = -915883387152444.0_R_P/2615348736000.0_R_P
      self%b(7) = 1226443086129408.0_R_P/2615348736000.0_R_P
      self%b(8) = -1233589244941764.0_R_P/2615348736000.0_R_P
      self%b(9) = 932884546055895.0_R_P/2615348736000.0_R_P
      self%b(10) = -524924579905150.0_R_P/2615348736000.0_R_P
      self%b(11) = 214696591002612.0_R_P/2615348736000.0_R_P
      self%b(12) = -61497552797274.0_R_P/2615348736000.0_R_P
      self%b(13) = 13064406523627.0_R_P/2615348736000.0_R_P
    case(14)
      self%b(1) = -1382741929621.0_R_P/5230697472000.0_R_P
      self%b(2) = 19382853593787.0_R_P/5230697472000.0_R_P
      self%b(3) = -126174972681906.0_R_P/5230697472000.0_R_P
      self%b(4) = 505586141196430.0_R_P/5230697472000.0_R_P
      self%b(5) = -1393306307155755.0_R_P/5230697472000.0_R_P
      self%b(6) = 2793869602879077.0_R_P/5230697472000.0_R_P
      self%b(7) = -4204551925534524.0_R_P/5230697472000.0_R_P
      self%b(8) = 4825671323488452.0_R_P/5230697472000.0_R_P
      self%b(9) = -4246767353305755.0_R_P/5230697472000.0_R_P
      self%b(10) = 2854429571790805.0_R_P/5230697472000.0_R_P
      self%b(11) = -1445313351681906.0_R_P/5230697472000.0_R_P
      self%b(12) = 537247052515662.0_R_P/5230697472000.0_R_P
      self%b(13) = -140970750679621.0_R_P/5230697472000.0_R_P
      self%b(14) = 27511554976875.0_R_P/5230697472000.0_R_P
    case(15)
      self%b(1) = 8164168737599.0_R_P/31384184832000.0_R_P
      self%b(2) = -122594813904112.0_R_P/31384184832000.0_R_P
      self%b(3) = 859236476684231.0_R_P/31384184832000.0_R_P
      self%b(4) = -3728807256577472.0_R_P/31384184832000.0_R_P
      self%b(5) = 11205849753515179.0_R_P/31384184832000.0_R_P
      self%b(6) = -24704503655607728.0_R_P/31384184832000.0_R_P
      self%b(7) = 41280216336284259.0_R_P/31384184832000.0_R_P
      self%b(8) = -53246738660646912.0_R_P/31384184832000.0_R_P
      self%b(9) = 53471026659940509.0_R_P/31384184832000.0_R_P
      self%b(10) = -41825269932507728.0_R_P/31384184832000.0_R_P
      self%b(11) = 25298910337081429.0_R_P/31384184832000.0_R_P
      self%b(12) = -11643637530577472.0_R_P/31384184832000.0_R_P
      self%b(13) = 3966421670215481.0_R_P/31384184832000.0_R_P
      self%b(14) = -960122866404112.0_R_P/31384184832000.0_R_P
      self%b(15) = 173233498598849.0_R_P/31384184832000.0_R_P
    case(16)
      self%b(1) = -16088129229375.0_R_P/62768369664000.0_R_P
      self%b(2) = 257650275915823.0_R_P/62768369664000.0_R_P
      self%b(3) = -1934443196892599.0_R_P/62768369664000.0_R_P
      self%b(4) = 9038571752734087.0_R_P/62768369664000.0_R_P
      self%b(5) = -29417910911251819.0_R_P/62768369664000.0_R_P
      self%b(6) = 70724351582843483.0_R_P/62768369664000.0_R_P
      self%b(7) = -129930094104237331.0_R_P/62768369664000.0_R_P
      self%b(8) = 186087544263596643.0_R_P/62768369664000.0_R_P
      self%b(9) = -210020588912321949.0_R_P/62768369664000.0_R_P
      self%b(10) = 187463140112902893.0_R_P/62768369664000.0_R_P
      self%b(11) = -131963191940828581.0_R_P/62768369664000.0_R_P
      self%b(12) = 72558117072259733.0_R_P/62768369664000.0_R_P
      self%b(13) = -30607373860520569.0_R_P/62768369664000.0_R_P
      self%b(14) = 9622096909515337.0_R_P/62768369664000.0_R_P
      self%b(15) = -2161567671248849.0_R_P/62768369664000.0_R_P
      self%b(16) = 362555126427073.0_R_P/62768369664000.0_R_P
    endselect
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
  class(adams_bashforth_integrator), intent(INOUT) :: self !< AB integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = 0
  if (allocated(self%b)) deallocate(self%b)
  self%error = 0
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

  elemental function min_steps()
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

  elemental function max_steps()
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
