!< FOODIE integrator: provide an implicit class of Adams-Moutlon multi-step schemes, from 1st to 16th order accurate.
module foodie_integrator_adams_moulton
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE integrator: provide an implicit class of Adams-Moutlon multi-step schemes, from 1st to 16th order accurate.
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
!< The schemes are implicit. The coefficients *b* define the actual scheme, that is selected accordingly to the number of
!< **steps** used. Currently, the schemes provided have steps number in *[0, 15]*. Note that the scheme using only 1 step reverts
!< to Implciti Backwarad Euler. The formal order of accuracy varies consistently in *[1st, 16th]* order.
!<
!<#### Bibliography
!< [1] *Cowell Type Numerical Integration As Applied to Satellite Orbit Computation*, J. L. Maury Jr.,
!< G. P. Segal, X-553-69-46, April 1969, [NASA-TM-X-63542](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19690017325.pdf).
!<
!< [2] *Linear multistep method*, [wikipedia article](https://en.wikipedia.org/wiki/Linear_multistep_method).
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie_adt_integrand, only : integrand
use foodie_kinds, only : R_P, I_P
use foodie_utils, only : is_admissible
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: adams_moulton_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
character(len=99), parameter :: supported_steps='0-15' !< List of supported steps number. Valid format is `1-2,4,9-23...`.

type :: adams_moulton_integrator
  !< FOODIE integrator: provide an explicit class of Adams-Moulton multi-step schemes, from 1st to 16th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *b* coefficients) before used.
  !<
  !< ### List of errors status
  !<+ error=0 => no error;
  !<+ error=1 => bad (unsupported) number of required time steps;
  private
  integer(I_P)           :: steps=-1 !< Number of time steps.
  real(R_P), allocatable :: b(:)     !< \(b\) coefficients.
  integer(I_P)           :: error=0  !< Error status flag: trap occurrences of errors.
  contains
    procedure, pass(self), public :: init            !< Initialize (create) the integrator.
    procedure, pass(self), public :: destroy         !< Destroy the integrator.
    procedure, pass(self), public :: integrate       !< Integrate integrand field.
    procedure, pass(self), public :: update_previous !< Cyclic update previous time steps.
    procedure, nopass,     public :: is_supported    !< Check if the queried number of steps is supported or not.
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
  integer(I_P),                    intent(IN)    :: steps !< Number of time steps used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (self%is_supported(steps)) then
    self%steps = steps
    if (allocated(self%b)) deallocate(self%b) ; allocate(self%b(0:steps)) ; self%b = 0.0_R_P
    select case(steps)
    case(0)
      ! AM(0) Bacward-Euler
      self%b(0) = 1.0_R_P
    case(1)
      self%b(0) = 1.0_R_P/2.0_R_P
      self%b(1) = 1.0_R_P/2.0_R_P
    case(2)
      self%b(0) = -1.0_R_P/12.0_R_P
      self%b(1) = 8.0_R_P/12.0_R_P
      self%b(2) = 5.0_R_P/12.0_R_P
    case(3)
      self%b(0) = 1.0_R_P/24.0_R_P
      self%b(1) = -5.0_R_P/24.0_R_P
      self%b(2) = 19.0_R_P/24.0_R_P
      self%b(3) = 9.0_R_P/24.0_R_P
    case(4)
      self%b(0) = -19.0_R_P/720.0_R_P
      self%b(1) = 106.0_R_P/720.0_R_P
      self%b(2) = -264.0_R_P/720.0_R_P
      self%b(3) = 646.0_R_P/720.0_R_P
      self%b(4) = 251.0_R_P/720.0_R_P
    case(5)
      self%b(0) = 27.0_R_P/1440.0_R_P
      self%b(1) = -173.0_R_P/1440.0_R_P
      self%b(2) = 482.0_R_P/1440.0_R_P
      self%b(3) = -798.0_R_P/1440.0_R_P
      self%b(4) = 1427.0_R_P/1440.0_R_P
      self%b(5) = 475.0_R_P/1440.0_R_P
    case(6)
      self%b(0) = -863.0_R_P/60480.0_R_P
      self%b(1) = 6312.0_R_P/60480.0_R_P
      self%b(2) = -20211.0_R_P/60480.0_R_P
      self%b(3) = 37504.0_R_P/60480.0_R_P
      self%b(4) = -46461.0_R_P/60480.0_R_P
      self%b(5) = 65112.0_R_P/60480.0_R_P
      self%b(6) = 19087.0_R_P/60480.0_R_P
    case(7)
      self%b(0) = 1375.0_R_P/120960.0_R_P
      self%b(1) = -11351.0_R_P/120960.0_R_P
      self%b(2) = 41499.0_R_P/120960.0_R_P
      self%b(3) = -88547.0_R_P/120960.0_R_P
      self%b(4) = 123133.0_R_P/120960.0_R_P
      self%b(5) = -121797.0_R_P/120960.0_R_P
      self%b(6) = 139849.0_R_P/120960.0_R_P
      self%b(7) = 36799.0_R_P/120960.0_R_P
    case(8)
      self%b(0) = -33953.0_R_P/3628800.0_R_P
      self%b(1) = 312874.0_R_P/3628800.0_R_P
      self%b(2) = -1291214.0_R_P/3628800.0_R_P
      self%b(3) = 3146338.0_R_P/3628800.0_R_P
      self%b(4) = -5033120.0_R_P/3628800.0_R_P
      self%b(5) = 5595358.0_R_P/3628800.0_R_P
      self%b(6) = -4604594.0_R_P/3628800.0_R_P
      self%b(7) = 4467094.0_R_P/3628800.0_R_P
      self%b(8) = 1070017.0_R_P/3628800.0_R_P
    case(9)
      self%b(0) = 57281.0_R_P/7257600.0_R_P
      self%b(1) = -583435.0_R_P/7257600.0_R_P
      self%b(2) = 2687864.0_R_P/7257600.0_R_P
      self%b(3) = -7394032.0_R_P/7257600.0_R_P
      self%b(4) = 13510082.0_R_P/7257600.0_R_P
      self%b(5) = -17283646.0_R_P/7257600.0_R_P
      self%b(6) = 16002320.0_R_P/7257600.0_R_P
      self%b(7) = -11271304.0_R_P/7257600.0_R_P
      self%b(8) = 9449717.0_R_P/7257600.0_R_P
      self%b(9) = 2082753.0_R_P/7257600.0_R_P
    case(10)
      self%b(0) = -3250433.0_R_P/479001600.0_R_P
      self%b(1) = 36284876.0_R_P/479001600.0_R_P
      self%b(2) = -184776195.0_R_P/479001600.0_R_P
      self%b(3) = 567450984.0_R_P/479001600.0_R_P
      self%b(4) = -1170597042.0_R_P/479001600.0_R_P
      self%b(5) = 1710774528.0_R_P/479001600.0_R_P
      self%b(6) = -1823311566.0_R_P/479001600.0_R_P
      self%b(7) = 1446205080.0_R_P/479001600.0_R_P
      self%b(8) = -890175549.0_R_P/479001600.0_R_P
      self%b(9) = 656185652.0_R_P/479001600.0_R_P
      self%b(10) = 134211265.0_R_P/479001600.0_R_P
    case(11)
      self%b(0) = 5675265.0_R_P/958003200.0_R_P
      self%b(1) = -68928781.0_R_P/958003200.0_R_P
      self%b(2) = 384709327.0_R_P/958003200.0_R_P
      self%b(3) = -1305971115.0_R_P/958003200.0_R_P
      self%b(4) = 3007739418.0_R_P/958003200.0_R_P
      self%b(5) = -4963166514.0_R_P/958003200.0_R_P
      self%b(6) = 6043521486.0_R_P/958003200.0_R_P
      self%b(7) = -5519460582.0_R_P/958003200.0_R_P
      self%b(8) = 3828828885.0_R_P/958003200.0_R_P
      self%b(9) = -2092490673.0_R_P/958003200.0_R_P
      self%b(10) = 1374799219.0_R_P/958003200.0_R_P
      self%b(11) = 262747265.0_R_P/958003200.0_R_P
    case(12)
      self%b(0) = -13695779093.0_R_P/2615348736000.0_R_P
      self%b(1) = 179842822566.0_R_P/2615348736000.0_R_P
      self%b(2) = -1092096992268.0_R_P/2615348736000.0_R_P
      self%b(3) = 4063327863170.0_R_P/2615348736000.0_R_P
      self%b(4) = -10344711794985.0_R_P/2615348736000.0_R_P
      self%b(5) = 19058185652796.0_R_P/2615348736000.0_R_P
      self%b(6) = -26204344465152.0_R_P/2615348736000.0_R_P
      self%b(7) = 27345870698436.0_R_P/2615348736000.0_R_P
      self%b(8) = -21847538039895.0_R_P/2615348736000.0_R_P
      self%b(9) = 13465774256510.0_R_P/2615348736000.0_R_P
      self%b(10) = -6616420957428.0_R_P/2615348736000.0_R_P
      self%b(11) = 3917551216986.0_R_P/2615348736000.0_R_P
      self%b(12) = 703604254357.0_R_P/2615348736000.0_R_P
    case(13)
      self%b(0) = 24466579093.0_R_P/5230697472000.0_R_P
      self%b(1) = -345457086395.0_R_P/5230697472000.0_R_P
      self%b(2) = 2268078814386.0_R_P/5230697472000.0_R_P
      self%b(3) = -9181635605134.0_R_P/5230697472000.0_R_P
      self%b(4) = 25620259777835.0_R_P/5230697472000.0_R_P
      self%b(5) = -52177910882661.0_R_P/5230697472000.0_R_P
      self%b(6) = 80101021029180.0_R_P/5230697472000.0_R_P
      self%b(7) = -94393338653892.0_R_P/5230697472000.0_R_P
      self%b(8) = 86180228689563.0_R_P/5230697472000.0_R_P
      self%b(9) = -61188680131285.0_R_P/5230697472000.0_R_P
      self%b(10) = 33928990133618.0_R_P/5230697472000.0_R_P
      self%b(11) = -15141235084110.0_R_P/5230697472000.0_R_P
      self%b(12) = 8153167962181.0_R_P/5230697472000.0_R_P
      self%b(13) = 1382741929621.0_R_P/5230697472000.0_R_P
    case(14)
      self%b(0) = -132282840127.0_R_P/31384184832000.0_R_P
      self%b(1) = 1998759236336.0_R_P/31384184832000.0_R_P
      self%b(2) = -14110480969927.0_R_P/31384184832000.0_R_P
      self%b(3) = 61759426692544.0_R_P/31384184832000.0_R_P
      self%b(4) = -187504936597931.0_R_P/31384184832000.0_R_P
      self%b(5) = 418551804601264.0_R_P/31384184832000.0_R_P
      self%b(6) = -710312834197347.0_R_P/31384184832000.0_R_P
      self%b(7) = 934600833490944.0_R_P/31384184832000.0_R_P
      self%b(8) = -963605400824733.0_R_P/31384184832000.0_R_P
      self%b(9) = 781911618071632.0_R_P/31384184832000.0_R_P
      self%b(10) = -499547203754837.0_R_P/31384184832000.0_R_P
      self%b(11) = 251724894607936.0_R_P/31384184832000.0_R_P
      self%b(12) = -102885148956217.0_R_P/31384184832000.0_R_P
      self%b(13) = 50770967534864.0_R_P/31384184832000.0_R_P
      self%b(14) = 8164168737599.0_R_P/31384184832000.0_R_P
    case(15)
      self%b(0) = 240208245823.0_R_P/62768369664000.0_R_P
      self%b(1) = -3867689367599.0_R_P/62768369664000.0_R_P
      self%b(2) = 29219384284087.0_R_P/62768369664000.0_R_P
      self%b(3) = -137515713789319.0_R_P/62768369664000.0_R_P
      self%b(4) = 451403108933483.0_R_P/62768369664000.0_R_P
      self%b(5) = -1096355235402331.0_R_P/62768369664000.0_R_P
      self%b(6) = 2039345879546643.0_R_P/62768369664000.0_R_P
      self%b(7) = -2966365730265699.0_R_P/62768369664000.0_R_P
      self%b(8) = 3414941728852893.0_R_P/62768369664000.0_R_P
      self%b(9) = -3129453071993581.0_R_P/62768369664000.0_R_P
      self%b(10) = 2285168598349733.0_R_P/62768369664000.0_R_P
      self%b(11) = -1326978663058069.0_R_P/62768369664000.0_R_P
      self%b(12) = 612744541065337.0_R_P/62768369664000.0_R_P
      self%b(13) = -230992163723849.0_R_P/62768369664000.0_R_P
      self%b(14) = 105145058757073.0_R_P/62768369664000.0_R_P
      self%b(15) = 16088129229375.0_R_P/62768369664000.0_R_P
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
  class(adams_moulton_integrator), intent(INOUT) :: self !< AB integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = -1
  if (allocated(self%b)) deallocate(self%b)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine integrate(self, U, previous, Dt, t, iterations, autoupdate)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with Adams-Moulton class scheme.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(adams_moulton_integrator), intent(IN)    :: self         !< Actual AB integrator.
  class(integrand),                intent(INOUT) :: U            !< Field to be integrated.
  class(integrand),                intent(INOUT) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                       intent(IN)    :: Dt           !< Time steps.
  real(R_P),                       intent(IN)    :: t(:)         !< Times.
  integer(I_P), optional,          intent(IN)    :: iterations   !< Fixed point iterations.
  logical,      optional,          intent(IN)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  logical                                        :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  class(integrand), allocatable                  :: delta        !< Delta RHS for fixed point iterations.
  integer(I_P)                                   :: s            !< Steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  if (self%steps>0) then
    if (present(iterations)) then ! perform fixed point iterations
      allocate(delta, source=previous(self%steps))
      do s=0, self%steps - 1
        delta = delta + previous(s+1)%t(t=t(s+1)) * (Dt * self%b(s))
      enddo
      do s=1, iterations
        U = delta + U%t(t=t(self%steps) + Dt) * (Dt * self%b(self%steps))
      enddo
    else
      U = previous(self%steps) + U%t(t=t(self%steps) + Dt) * (Dt * self%b(self%steps))
      do s=0, self%steps - 1
        U = U + previous(s+1)%t(t=t(s+1)) * (Dt * self%b(s))
      enddo
    endif
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
  type(adams_moulton_integrator), intent(INOUT) :: self !< AB integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule foodie_integrator_adams_moulton
