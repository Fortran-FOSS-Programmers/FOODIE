#include "preprocessor_macros.h"
!< FOODIE integrator: provide an implicit class of Adams-Moutlon multi-step schemes, from 1st to 16th order accurate.

module foodie_integrator_adams_moulton
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

use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrand_object, only : integrand_object
use foodie_integrator_multistep_object, only : integrator_multistep_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_adams_moulton

character(len=99), parameter :: class_name_='adams_moulton'                         !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:16)=[trim(class_name_)//'_0 ', &
                                                          trim(class_name_)//'_1 ', &
                                                          trim(class_name_)//'_2 ', &
                                                          trim(class_name_)//'_3 ', &
                                                          trim(class_name_)//'_4 ', &
                                                          trim(class_name_)//'_5 ', &
                                                          trim(class_name_)//'_6 ', &
                                                          trim(class_name_)//'_7 ', &
                                                          trim(class_name_)//'_8 ', &
                                                          trim(class_name_)//'_9 ', &
                                                          trim(class_name_)//'_10', &
                                                          trim(class_name_)//'_11', &
                                                          trim(class_name_)//'_12', &
                                                          trim(class_name_)//'_13', &
                                                          trim(class_name_)//'_14', &
                                                          trim(class_name_)//'_15'] !< List of supported schemes.

logical, parameter :: has_fast_mode_=.true.  !< Flag to check if integrator provides *fast mode* integrate.

type, extends(integrator_multistep_object) :: integrator_adams_moulton
  !< FOODIE integrator: provide an explicit class of Adams-Moulton multi-step schemes, from 1st to 16th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *b* coefficients) before used.
   !<
   !< @note The time steps `Dt(1:steps)` passed to the integrate methods must be identical: this integrator supports only
   !< fixed time steps.
  private
  real(R_P), allocatable :: b(:) !< \(b\) coefficients.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: has_fast_mode        !< Return .true. if the integrator class has *fast mode* integrate.
    procedure, pass(self) :: integrate            !< Integrate integrand field.
    procedure, pass(self) :: integrate_fast       !< Integrate integrand field, fast mode.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy    !< Destroy the integrator.
    procedure, pass(self) :: initialize !< Initialize (create) the integrator.
endtype integrator_adams_moulton

contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_adams_moulton), intent(in) :: self       !< Integrator.
  character(len=99)                           :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  elemental function has_fast_mode(self)
  !< Return .true. if the integrator class has *fast mode* integrate.
  class(integrator_adams_moulton), intent(in) :: self          !< Integrator.
  logical                                     :: has_fast_mode !< Inquire result.

  has_fast_mode = has_fast_mode_
  endfunction has_fast_mode

  _PURE_ subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_adams_moulton), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),        intent(in)    :: rhs !< Right hand side.

  call lhs%assign_multistep(rhs=rhs)
  select type(rhs)
  class is (integrator_adams_moulton)
                          lhs%steps = rhs%steps
    if (allocated(rhs%b)) lhs%b     = rhs%b
  endselect
  endsubroutine integr_assign_integr

  subroutine integrate(self, U, Dt, t)
  !< Integrate field with Adams-Moulton class scheme.
  class(integrator_adams_moulton), intent(inout) :: self  !< Integrator.
  class(integrand_object),         intent(inout) :: U     !< Field to be integrated.
  real(R_P),                       intent(in)    :: Dt    !< Time step.
  real(R_P),                       intent(in)    :: t     !< Time.
  class(integrand_object), allocatable           :: delta !< Delta RHS for fixed point iterations.
  integer(I_P)                                   :: s     !< Steps counter.

  if ( self%steps>0 ) then
    if (self%iterations > 0) then ! perform fixed point iterations
      allocate(delta, mold=U)
      delta = self%previous(self%steps)
      do s=0, self%steps - 1
        delta = delta + (self%previous(s+1)%t(t=self%t(s+1)) * (Dt * self%b(s)))
      enddo
      do s=1, self%iterations
        U = delta + (U%t(t=self%t(self%steps) + Dt) * (Dt * self%b(self%steps)))
      enddo
    else
      U = self%previous(self%steps) + (U%t(t=self%t(self%steps) + Dt) * (Dt * self%b(self%steps)))
      do s=0, self%steps - 1
        U = U + (self%previous(s+1)%t(t=self%t(s+1)) * (Dt * self%b(s)))
      enddo
    endif
    if (self%autoupdate) call self%update_previous(U=U, previous=self%previous, Dt=Dt, t=t, previous_t=self%t)
  else
    U = U + (U%t(t=t) * (Dt * self%b(0)))
  endif
  endsubroutine integrate

  subroutine integrate_fast(self, U, Dt, t)
  !< Integrate field with Adams-Moulton class scheme, fast mode.
  class(integrator_adams_moulton), intent(inout) :: self  !< Integrator.
  class(integrand_object),         intent(inout) :: U     !< Field to be integrated.
  real(R_P),                       intent(in)    :: Dt    !< Time step.
  real(R_P),                       intent(in)    :: t     !< Time.
  class(integrand_object), allocatable           :: delta !< Delta RHS for fixed point iterations.
  integer(I_P)                                   :: s     !< Steps counter.

  if (self%steps>0) then
    if (self%iterations > 0) then ! perform fixed point iterations
      allocate(delta, mold=U)
      delta = self%previous(self%steps)
      do s=0, self%steps - 1
        self%buffer = self%previous(s+1)
        call self%buffer%t_fast(t=self%t(s+1))
        call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%b(s))
        call delta%add_fast(lhs=delta, rhs=self%buffer)
      enddo
      do s=1, self%iterations
        self%buffer = U
        call self%buffer%t_fast(t=self%t(self%steps) + Dt)
        call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%b(self%steps))
        call U%add_fast(lhs=delta, rhs=self%buffer)
      enddo
    else
      self%buffer = U
      call self%buffer%t_fast(t=self%t(self%steps) + Dt)
      call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%b(self%steps))
      call U%add_fast(lhs=self%previous(self%steps), rhs=self%buffer)
      do s=0, self%steps - 1
        self%buffer = self%previous(s+1)
        call self%buffer%t_fast(t=self%t(s+1))
        call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%b(s))
        call U%add_fast(lhs=U, rhs=self%buffer)
      enddo
    endif
    if (self%autoupdate) call self%update_previous(U=U, previous=self%previous, Dt=Dt, t=t, previous_t=self%t)
  else
    self%buffer = U
    call self%buffer%t_fast(t=t)
    call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%b(0))
    call U%add_fast(lhs=U, rhs=self%buffer)
  endif
  endsubroutine integrate_fast

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_adams_moulton), intent(in) :: self         !< Integrator.
  character(*),                    intent(in) :: scheme       !< Selected scheme.
  logical                                     :: is_supported !< Inquire result.
  integer(I_P)                                :: s            !< Counter.

  is_supported = .false.
  do s=lbound(supported_schemes_, dim=1), ubound(supported_schemes_, dim=1)
    if (trim(adjustl(scheme)) == trim(adjustl(supported_schemes_(s)))) then
      is_supported = .true.
      return
    endif
  enddo
  endfunction is_supported

  pure function supported_schemes(self) result(schemes)
  !< Return the list of supported schemes.
  class(integrator_adams_moulton), intent(in) :: self       !< Integrator.
  character(len=99), allocatable              :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_adams_moulton), intent(INOUT) :: self !< Integrator.

  call self%destroy_multistep
  if (allocated(self%b)) deallocate(self%b)
  endsubroutine destroy

  subroutine initialize(self, scheme, iterations, autoupdate, U, stop_on_fail)
  !< Create the actual Adams-Moulton integrator: initialize the *b* coefficients.
  class(integrator_adams_moulton), intent(inout)        :: self         !< Integrator.
  character(*),                    intent(in)           :: scheme       !< Selected scheme.
  integer(I_P),                    intent(in), optional :: iterations   !< Implicit iterations.
  logical,                         intent(in), optional :: autoupdate   !< Enable cyclic autoupdate of previous time steps.
  class(integrand_object),         intent(in), optional :: U            !< Integrand molding prototype.
  logical,                         intent(in), optional :: stop_on_fail !< Stop execution if initialization fail.

  if (self%is_supported(scheme=scheme)) then
     call self%destroy
     self%description_ = trim(adjustl(scheme))
    select case(trim(adjustl(scheme)))
    case('adams_moulton_0')
      self%steps = 0 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
      self%b(0) = 1.0_R_P
    case('adams_moulton_1')
      self%steps = 1 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
      self%b(0) = 1.0_R_P/2.0_R_P
      self%b(1) = 1.0_R_P/2.0_R_P
    case('adams_moulton_2')
      self%steps = 2 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
      self%b(0) = -1.0_R_P/12.0_R_P
      self%b(1) = 8.0_R_P/12.0_R_P
      self%b(2) = 5.0_R_P/12.0_R_P
    case('adams_moulton_3')
      self%steps = 3 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
      self%b(0) = 1.0_R_P/24.0_R_P
      self%b(1) = -5.0_R_P/24.0_R_P
      self%b(2) = 19.0_R_P/24.0_R_P
      self%b(3) = 9.0_R_P/24.0_R_P
    case('adams_moulton_4')
      self%steps = 4 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
      self%b(0) = -19.0_R_P/720.0_R_P
      self%b(1) = 106.0_R_P/720.0_R_P
      self%b(2) = -264.0_R_P/720.0_R_P
      self%b(3) = 646.0_R_P/720.0_R_P
      self%b(4) = 251.0_R_P/720.0_R_P
    case('adams_moulton_5')
      self%steps = 5 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
      self%b(0) = 27.0_R_P/1440.0_R_P
      self%b(1) = -173.0_R_P/1440.0_R_P
      self%b(2) = 482.0_R_P/1440.0_R_P
      self%b(3) = -798.0_R_P/1440.0_R_P
      self%b(4) = 1427.0_R_P/1440.0_R_P
      self%b(5) = 475.0_R_P/1440.0_R_P
    case('adams_moulton_6')
      self%steps = 6 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
      self%b(0) = -863.0_R_P/60480.0_R_P
      self%b(1) = 6312.0_R_P/60480.0_R_P
      self%b(2) = -20211.0_R_P/60480.0_R_P
      self%b(3) = 37504.0_R_P/60480.0_R_P
      self%b(4) = -46461.0_R_P/60480.0_R_P
      self%b(5) = 65112.0_R_P/60480.0_R_P
      self%b(6) = 19087.0_R_P/60480.0_R_P
    case('adams_moulton_7')
      self%steps = 7 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
      self%b(0) = 1375.0_R_P/120960.0_R_P
      self%b(1) = -11351.0_R_P/120960.0_R_P
      self%b(2) = 41499.0_R_P/120960.0_R_P
      self%b(3) = -88547.0_R_P/120960.0_R_P
      self%b(4) = 123133.0_R_P/120960.0_R_P
      self%b(5) = -121797.0_R_P/120960.0_R_P
      self%b(6) = 139849.0_R_P/120960.0_R_P
      self%b(7) = 36799.0_R_P/120960.0_R_P
    case('adams_moulton_8')
      self%steps = 8 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
      self%b(0) = -33953.0_R_P/3628800.0_R_P
      self%b(1) = 312874.0_R_P/3628800.0_R_P
      self%b(2) = -1291214.0_R_P/3628800.0_R_P
      self%b(3) = 3146338.0_R_P/3628800.0_R_P
      self%b(4) = -5033120.0_R_P/3628800.0_R_P
      self%b(5) = 5595358.0_R_P/3628800.0_R_P
      self%b(6) = -4604594.0_R_P/3628800.0_R_P
      self%b(7) = 4467094.0_R_P/3628800.0_R_P
      self%b(8) = 1070017.0_R_P/3628800.0_R_P
    case('adams_moulton_9')
      self%steps = 9 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_moulton_10')
      self%steps = 10 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_moulton_11')
      self%steps = 11 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_moulton_12')
      self%steps = 12 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_moulton_13')
      self%steps = 13 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_moulton_14')
      self%steps = 14 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_moulton_15')
      self%steps = 15 ; allocate(self%b(0:self%steps)) ; self%b = 0.0_R_P
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
     self%autoupdate = .true. ; if (present(autoupdate)) self%autoupdate = autoupdate
     self%iterations = 1 ; if (present(iterations)) self%iterations = iterations
     self%registers = self%steps
     if (present(U)) call self%allocate_integrand_members(U=U)
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=stop_on_fail)
  endif
  endsubroutine initialize
endmodule foodie_integrator_adams_moulton
