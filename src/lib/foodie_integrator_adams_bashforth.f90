!< FOODIE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 16th order accurate.

module foodie_integrator_adams_bashforth
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

use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrand_object, only : integrand_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_adams_bashforth

character(len=99), parameter :: class_name_='adams_bashforth'                       !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:16)=[trim(class_name_)//'_1 ', &
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
                                                          trim(class_name_)//'_15', &
                                                          trim(class_name_)//'_16'] !< List of supported schemes.

logical, parameter :: has_fast_mode_=.true.  !< Flag to check if integrator provides *fast mode* integrate.
logical, parameter :: is_multistage_=.false. !< Flag to check if integrator is multistage.
logical, parameter :: is_multistep_=.true.   !< Flag to check if integrator is multistep.

type, extends(integrator_object) :: integrator_adams_bashforth
  !< FOODIE integrator: provide an explicit class of Adams-Bashforth multi-step schemes, from 1st to 16th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *b* coefficients) before used.
  private
  integer(I_P)           :: steps=0 !< Number of time steps.
  real(R_P), allocatable :: b(:)    !< *b* coefficients.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(self) :: has_fast_mode        !< Return .true. if the integrator class has *fast mode* integrate.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: is_multistage        !< Return .true. for multistage integrator.
    procedure, pass(self) :: is_multistep         !< Return .true. for multistep integrator.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: stages_number        !< Return number of stages used.
    procedure, pass(self) :: steps_number         !< Return number of steps used.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy         !< Destroy the integrator.
    procedure, pass(self) :: initialize      !< Initialize (create) the integrator.
    procedure, pass(self) :: integrate       !< Integrate integrand field.
    procedure, pass(self) :: integrate_fast  !< Integrate integrand field, fast mode.
    procedure, pass(self) :: update_previous !< Cyclic update previous time steps.
endtype integrator_adams_bashforth

contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_adams_bashforth), intent(in) :: self       !< Integrator.
  character(len=99)                             :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_adams_bashforth), intent(in)           :: self             !< Integrator.
  character(*),                      intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                           :: desc             !< Description.
  character(len=:), allocatable                           :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                             :: NL=new_line('a') !< New line character.
  integer(I_P)                                            :: s                !< Counter.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = ''
  desc = desc//prefix_//'Adams-Bashforth multi-step schemes class'//NL
  desc = desc//prefix_//'  Supported schemes:'//NL
  do s=lbound(supported_schemes_, dim=1), ubound(supported_schemes_, dim=1) - 1
    desc = desc//prefix_//'    + '//supported_schemes_(s)//NL
  enddo
  desc = desc//prefix_//'    + '//supported_schemes_(ubound(supported_schemes_, dim=1))
  endfunction description

  elemental function has_fast_mode(self)
  !< Return .true. if the integrator class has *fast mode* integrate.
  class(integrator_adams_bashforth), intent(in) :: self          !< Integrator.
  logical                                       :: has_fast_mode !< Inquire result.

  has_fast_mode = has_fast_mode_
  endfunction has_fast_mode

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_adams_bashforth), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),          intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  select type(rhs)
  class is (integrator_adams_bashforth)
                          lhs%steps = rhs%steps
    if (allocated(rhs%b)) lhs%b     = rhs%b
  endselect
  endsubroutine integr_assign_integr

  elemental function is_multistage(self)
  !< Return .true. for multistage integrator.
  class(integrator_adams_bashforth), intent(in) :: self          !< Integrator.
  logical                                       :: is_multistage !< Inquire result.

  is_multistage = is_multistage_
  endfunction is_multistage

  elemental function is_multistep(self)
  !< Return .true. for multistage integrator.
  class(integrator_adams_bashforth), intent(in) :: self         !< Integrator.
  logical                                       :: is_multistep !< Inquire result.

  is_multistep = is_multistep_
  endfunction is_multistep

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_adams_bashforth), intent(in) :: self         !< Integrator.
  character(*),                      intent(in) :: scheme       !< Selected scheme.
  logical                                       :: is_supported !< Inquire result.
  integer(I_P)                                  :: s            !< Counter.

  is_supported = .false.
  do s=lbound(supported_schemes_, dim=1), ubound(supported_schemes_, dim=1)
    if (trim(adjustl(scheme)) == trim(adjustl(supported_schemes_(s)))) then
      is_supported = .true.
      return
    endif
  enddo
  endfunction is_supported

  elemental function stages_number(self)
  !< Return number of stages used.
  class(integrator_adams_bashforth), intent(in) :: self          !< Integrator.
  integer(I_P)                                  :: stages_number !< Number of stages used.

  stages_number = 0
  endfunction stages_number

  elemental function steps_number(self)
  !< Return number of steps used.
  class(integrator_adams_bashforth), intent(in) :: self         !< Integrator.
  integer(I_P)                                  :: steps_number !< Number of steps used.

  steps_number = self%steps
  endfunction steps_number

  pure function supported_schemes(self) result(schemes)
  !< Return the list of supported schemes.
  class(integrator_adams_bashforth), intent(in) :: self       !< Integrator.
  character(len=99), allocatable                :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_adams_bashforth), intent(inout) :: self !< Integrator.

  call self%destroy_abstract
  self%steps = 0
  if (allocated(self%b)) deallocate(self%b)
  endsubroutine destroy

  subroutine initialize(self, scheme)
  !< Create the actual Adams-Bashforth integrator: initialize the *b* coefficients.
  !<
  !< @note If the integrator is initialized with a bad (unsupported) number of required time steps the initialization fails and
  !< the integrator error status is updated consistently for external-provided errors handling.
  class(integrator_adams_bashforth), intent(inout) :: self   !< Integrator.
  character(*),                      intent(in)    :: scheme !< Selected scheme.

  if (self%is_supported(scheme=scheme)) then
    call self%destroy
    select case(trim(adjustl(scheme)))
    case('adams_bashforth_1')
      self%steps = 1 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%b(1) = 1.0_R_P
    case('adams_bashforth_2')
      self%steps = 2 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%b(1) = -1.0_R_P/2.0_R_P
      self%b(2) = 3.0_R_P/2.0_R_P
    case('adams_bashforth_3')
      self%steps = 3 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%b(1) = 5.0_R_P/12.0_R_P
      self%b(2) = -16.0_R_P/12.0_R_P
      self%b(3) = 23.0_R_P/12.0_R_P
    case('adams_bashforth_4')
      self%steps = 4 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%b(1) = -9.0_R_P/24.0_R_P
      self%b(2) = 37.0_R_P/24.0_R_P
      self%b(3) = -59.0_R_P/24.0_R_P
      self%b(4) = 55.0_R_P/24.0_R_P
    case('adams_bashforth_5')
      self%steps = 5 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%b(1) = 251.0_R_P/720.0_R_P
      self%b(2) = -1274.0_R_P/720.0_R_P
      self%b(3) = 2616.0_R_P/720.0_R_P
      self%b(4) = -2774.0_R_P/720.0_R_P
      self%b(5) = 1901.0_R_P/720.0_R_P
    case('adams_bashforth_6')
      self%steps = 6 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%b(1) = -475.0_R_P/1440.0_R_P
      self%b(2) = 2877.0_R_P/1440.0_R_P
      self%b(3) = -7298.0_R_P/1440.0_R_P
      self%b(4) = 9982.0_R_P/1440.0_R_P
      self%b(5) = -7923.0_R_P/1440.0_R_P
      self%b(6) = 4277.0_R_P/1440.0_R_P
    case('adams_bashforth_7')
      self%steps = 7 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%b(1) = 19087.0_R_P/60480.0_R_P
      self%b(2) = -134472.0_R_P/60480.0_R_P
      self%b(3) = 407139.0_R_P/60480.0_R_P
      self%b(4) = -688256.0_R_P/60480.0_R_P
      self%b(5) = 705549.0_R_P/60480.0_R_P
      self%b(6) = -447288.0_R_P/60480.0_R_P
      self%b(7) = 198721.0_R_P/60480.0_R_P
    case('adams_bashforth_8')
      self%steps = 8 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%b(1) = -36799.0_R_P/120960.0_R_P
      self%b(2) = 295767.0_R_P/120960.0_R_P
      self%b(3) = -1041723.0_R_P/120960.0_R_P
      self%b(4) = 2102243.0_R_P/120960.0_R_P
      self%b(5) = -2664477.0_R_P/120960.0_R_P
      self%b(6) = 2183877.0_R_P/120960.0_R_P
      self%b(7) = -1152169.0_R_P/120960.0_R_P
      self%b(8) = 434241.0_R_P/120960.0_R_P
    case('adams_bashforth_9')
      self%steps = 9 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%b(1) = 1070017.0_R_P/3628800.0_R_P
      self%b(2) = -9664106.0_R_P/3628800.0_R_P
      self%b(3) = 38833486.0_R_P/3628800.0_R_P
      self%b(4) = -91172642.0_R_P/3628800.0_R_P
      self%b(5) = 137968480.0_R_P/3628800.0_R_P
      self%b(6) = -139855262.0_R_P/3628800.0_R_P
      self%b(7) = 95476786.0_R_P/3628800.0_R_P
      self%b(8) = -43125206.0_R_P/3628800.0_R_P
      self%b(9) = 14097247.0_R_P/3628800.0_R_P
    case('adams_bashforth_10')
      self%steps = 10 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_bashforth_11')
      self%steps = 11 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_bashforth_12')
      self%steps = 12 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_bashforth_13')
      self%steps = 13 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_bashforth_14')
      self%steps = 14 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_bashforth_15')
      self%steps = 15 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
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
    case('adams_bashforth_16')
      self%steps = 16 ; allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
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
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=.true.)
  endif
  endsubroutine initialize

  subroutine integrate(self, U, previous, Dt, t, autoupdate)
  !< Integrate field with Adams-Bashforth class scheme.
  class(integrator_adams_bashforth), intent(in)    :: self         !< Integrator.
  class(integrand_object),           intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),           intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                         intent(in)    :: Dt           !< Time steps.
  real(R_P),                         intent(in)    :: t(:)         !< Times.
  logical, optional,                 intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  logical                                          :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  integer(I_P)                                     :: s            !< Steps counter.

  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  do s=1, self%steps
    U = U + (previous(s)%t(t=t(s)) * (Dt * self%b(s)))
  enddo
  if (autoupdate_) call self%update_previous(U=U, previous=previous)
  endsubroutine integrate

  subroutine integrate_fast(self, U, previous, buffer, Dt, t, autoupdate)
  !< Integrate field with Adams-Bashforth class scheme.
  class(integrator_adams_bashforth), intent(in)    :: self         !< Integrator.
  class(integrand_object),           intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),           intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  class(integrand_object),           intent(inout) :: buffer       !< Temporary buffer for doing fast operation.
  real(R_P),                         intent(in)    :: Dt           !< Time steps.
  real(R_P),                         intent(in)    :: t(:)         !< Times.
  logical, optional,                 intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  logical                                          :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  integer(I_P)                                     :: s            !< Steps counter.

  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  do s=1, self%steps
    buffer = previous(s)
    call buffer%t_fast(t=t(s))
    call buffer%multiply_fast(lhs=buffer, rhs=Dt * self%b(s))
    call U%add_fast(lhs=U, rhs=buffer)
  enddo
  if (autoupdate_) call self%update_previous(U=U, previous=previous)
  endsubroutine integrate_fast

  subroutine update_previous(self, U, previous)
  !< Cyclic update previous time steps.
  class(integrator_adams_bashforth), intent(in)    :: self         !< Integrator.
  class(integrand_object),           intent(in)    :: U            !< Field to be integrated.
  class(integrand_object),           intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  integer(I_P)                                     :: s            !< Steps counter.

  do s=1, self%steps - 1
    previous(s) = previous(s + 1)
  enddo
  previous(self%steps) = U
  endsubroutine update_previous
endmodule foodie_integrator_adams_bashforth
