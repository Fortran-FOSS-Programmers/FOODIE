!< FOODIE integrator: provide an explicit class of Linear SSP Runge-Kutta schemes, from 1st to s-th order accurate.

module foodie_integrator_runge_kutta_lssp
!< FOODIE integrator: provide an explicit class of Linear SSP Runge-Kutta schemes, from 1st to s-th order accurate.
!<
!< The integrators provided have the Total Variation Diminishing (TVD) property or the Strong Stability Preserving (SSP)
!< one. The schemes are explicit.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the class of schemes implemented are written in the following 2 forms:
!<
!<#### Ns stages, Ns-1 order
!<
!< $$ U^0 = U^n $$
!< $$ U^k = U^{k-1} + \frac{1}{2} \Delta t R(U^{k-1}) \quad k=1,..,Ns $$
!< $$ U^{Ns} = U^{Ns} + \frac{1}{2} \Delta t R(U^{Ns}) $$
!< $$ U^{n+1} = \sum_{k=1}^{Ns}\alpha_k U^k $$
!<
!< where Ns is the number of stages used and \(U^k\) is the \(k^{th}\) stage. The minimum number of stages is 2. The
!< coefficients \(\alpha\) are computed by means of the following recursive algorith:
!<
!<##### Computation of coefficients
!<
!<+ allocate coefficients array, `allocate(alpha(1:Ns))`;
!<+ initialize the first 2 elements with the coefficients of the two stages first order methods, `alpha(1)=0, alpha(2)=1`;
!<+ for `i` in [3,Ns] do:
!<  + `alpha(i) = (2 / i) * alpha(i-1)`
!<  + for `j` in [i-1,2,-1] do:
!<    + `alpha(j) = (2 / (j-1)) * alpha(j-1)`
!<  + for `j` in [2,i] do:
!<    + `alpha(1) = 1 - alpha(j)`
!<
!<#### Ns stages, Ns order
!<
!< $$ U^0 = U^n $$
!< $$ U^k = U^{k-1} + \Delta t R(U^{k-1}) \quad k=1,..,Ns $$
!< $$ U^{Ns} = U^{Ns} +\Delta t R(U^{Ns}) $$
!< $$ U^{n+1} = \sum_{k=1}^{Ns}\alpha_k U^k $$
!<
!< where Ns is the number of stages used and \(U^k\) is the \(k^{th}\) stage. The minimum number of stages is 1. The
!< coefficients \(\alpha\) are computed by means of the following recursive algorith:
!<
!<##### Computation of coefficients
!<
!<+ allocate coefficients array, `allocate(alpha(1:Ns))`;
!<+ initialize the first element with the coefficient of the one stage first order method, `alpha(1)=1`;
!<+ for `i` in [2,Ns] do:
!<  + `alpha(i) = 1 / i!`
!<  + for `j` in [i-1,2,-1] do:
!<    + `alpha(j) = (1 / (j-1)) * alpha(j-1)`
!<  + for `j` in [2,i] do:
!<    + `alpha(1) = 1 - alpha(j)`
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!<#### Bibliography
!<
!< [2] *Strong Stability Preserving Runge-Kutta and Multistep Time Discretizations*, S. Gottlieb, D. Ketcheson, C.W. Shu,
!< 2011, 978-981-4289-26-9, doi:10.1142/7498, World Scientific Publishing Co. Pte. Ltd.

use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrand_object, only : integrand_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_runge_kutta_lssp

character(len=99), parameter :: class_name_='runge_kutta_lssp'                                     !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:2)=[trim(class_name_)//'_stages_s_order_s_1', &
                                                         trim(class_name_)//'_stages_s_order_s  '] !< List of supported schemes.

logical, parameter :: has_fast_mode_=.false. !< Flag to check if integrator provides *fast mode* integrate.

type, extends(integrator_object) :: integrator_runge_kutta_lssp
  !< FOODIE integrator: provide an explicit class of Linear SSP Runge-Kutta schemes, from 1st to s-th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the RK coefficients) before used.
  private
  integer(I_P), public                    :: stages=0                          !< Number of stages.
  real(R_P), allocatable                  :: alpha(:)                          !< \(\alpha\) coefficients.
  procedure(integrate_interface), pointer :: integrate_ => integrate_order_s_1 !< Integrate integrand field.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(self) :: has_fast_mode        !< Return .true. if the integrator class has *fast mode* integrate.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy    !< Destroy the integrator.
    procedure, pass(self) :: initialize !< Initialize (create) the integrator.
    procedure, pass(self) :: integrate  !< Integrate integrand field.
    ! private methods
    procedure, pass(self), private :: initialize_order_s_1 !< Integrate integrator for (s-1)-th order formula.
    procedure, pass(self), private :: initialize_order_s   !< Integrate integrator for s-th order formula.
    procedure, pass(self), private :: integrate_order_s_1  !< Integrate integrand field by (s-1)-th order formula.
    procedure, pass(self), private :: integrate_order_s    !< Integrate integrand field by s-th order formula.
endtype integrator_runge_kutta_lssp

abstract interface
  !< Abstract interfaces of [[integrator_runge_kutta_lssp]] methods.
  subroutine integrate_interface(self, U, stage, Dt, t)
  !< Integrate field with Linear SSP Runge-Kutta scheme.
  import :: integrand_object, integrator_runge_kutta_lssp, R_P
  class(integrator_runge_kutta_lssp), intent(in)    :: self      !< Integrator.
  class(integrand_object),            intent(inout) :: U         !< Field to be integrated.
  class(integrand_object),            intent(inout) :: stage(1:) !< Runge-Kutta stages [1:stages].
  real(R_P),                          intent(in)    :: Dt        !< Time steps.
  real(R_P),                          intent(in)    :: t         !< Times.
  endsubroutine integrate_interface
endinterface

contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_runge_kutta_lssp), intent(in) :: self       !< Integrator.
  character(len=99)                              :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_runge_kutta_lssp), intent(in)           :: self             !< Integrator.
  character(*),                       intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                            :: desc             !< Description.
  character(len=:), allocatable                            :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                              :: NL=new_line('a') !< New line character.
  integer(I_P)                                             :: s                !< Counter.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = ''
  desc = desc//prefix_//'Linear SSP Runge-Kutta multi-stage schemes class'//NL
  desc = desc//prefix_//'  Supported schemes:'//NL
  do s=lbound(supported_schemes_, dim=1), ubound(supported_schemes_, dim=1) - 1
    desc = desc//prefix_//'    + '//supported_schemes_(s)//NL
  enddo
  desc = desc//prefix_//'    + '//supported_schemes_(ubound(supported_schemes_, dim=1))
  endfunction description

  elemental function has_fast_mode(self)
  !< Return .true. if the integrator class has *fast mode* integrate.
  class(integrator_runge_kutta_lssp), intent(in) :: self          !< Integrator.
  logical                                        :: has_fast_mode !< Inquire result.

  has_fast_mode = has_fast_mode_
  endfunction has_fast_mode

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_runge_kutta_lssp), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),           intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  select type(rhs)
  class is (integrator_runge_kutta_lssp)
                              lhs%stages = rhs%stages
    if (allocated(rhs%alpha)) lhs%alpha  = rhs%alpha
  endselect
  endsubroutine integr_assign_integr

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_runge_kutta_lssp), intent(in) :: self         !< Integrator.
  character(*),                       intent(in) :: scheme       !< Selected scheme.
  logical                                        :: is_supported !< Inquire result.
  integer(I_P)                                   :: s            !< Counter.

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
  class(integrator_runge_kutta_lssp), intent(in) :: self       !< Integrator.
  character(len=99), allocatable                 :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_runge_kutta_lssp), intent(inout) :: self !< Integrator.

  call self%destroy_abstract
  self%stages = 0
  if (allocated(self%alpha)) deallocate(self%alpha)
  self%integrate_ => integrate_order_s_1
  endsubroutine destroy

  subroutine initialize(self, scheme, stages, stop_on_fail)
  !< Create the actual RK integrator: initialize the Butcher' table coefficients.
  class(integrator_runge_kutta_lssp), intent(inout)        :: self         !< Integrator.
  character(*),                       intent(in)           :: scheme       !< Selected scheme.
  integer(I_P),                       intent(in), optional :: stages       !< Stages number.
  logical,                            intent(in), optional :: stop_on_fail !< Stop execution if initialization fail.

  if (self%is_supported(scheme=scheme)) then
    call self%destroy
    select case(trim(adjustl(scheme)))
    case('runge_kutta_lssp_stages_s_order_s_1')
      self%integrate_ => integrate_order_s_1
      self%stages = 2 ; if (present(stages)) self%stages = stages
      if (self%stages < 2) then
        error stop 'error: the number of stages of "runge_kutta_lssp_stages_s_order_s_1" must be >=2'
      endif
      allocate(self%alpha(1:self%stages)) ; self%alpha = 0._R_P
      call self%initialize_order_s_1
    case('runge_kutta_lssp_stages_s_order_s')
      self%integrate_ => integrate_order_s
      self%stages = 1 ; if (present(stages)) self%stages = stages
      if (self%stages < 1) then
        error stop 'error: the number of stages of "runge_kutta_lssp_stages_s_order_s" must be >=1'
      endif
      allocate(self%alpha(1:self%stages)) ; self%alpha = 0._R_P
      call self%initialize_order_s
    endselect
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=stop_on_fail)
  endif
  endsubroutine initialize

  subroutine integrate(self, U, stage, Dt, t)
  !< Integrate integrand field by Linear SSP Runge-Kutta methods.
  !<
  !< @note This method can be used **after** the integrator is created (i.e. the RK coefficients are initialized).
  class(integrator_runge_kutta_lssp), intent(in)    :: self      !< Integrator.
  class(integrand_object),            intent(inout) :: U         !< Field to be integrated.
  class(integrand_object),            intent(inout) :: stage(1:) !< Runge-Kutta stages [1:stages].
  real(R_P),                          intent(in)    :: Dt        !< Time step.
  real(R_P),                          intent(in)    :: t         !< Time.

  call self%integrate_(U=U, stage=stage, Dt=Dt, t=t)
  endsubroutine integrate

  ! private methods
  elemental subroutine initialize_order_s_1(self)
  !< Initialize integrator for (s-1)-th order formula.
  class(integrator_runge_kutta_lssp), intent(inout) :: self !< Integrator.
  integer(I_P)                                      :: i    !< Counter.
  integer(I_P)                                      :: j    !< Counter.

  associate(alpha=>self%alpha)
    alpha(1) = 0._R_P
    alpha(2) = 1._R_P
    do i=3, self%stages
      alpha(i) = 2._R_P / i * alpha(i-1)
      do j=i-1, 2, -1
        alpha(j) = 2._R_P / (j-1) * alpha(j-1)
      enddo
      alpha(1) = 1._R_P
      do j=2, i
        alpha(1) = alpha(1) - alpha(j)
      enddo
    enddo
  endassociate
  endsubroutine initialize_order_s_1

  elemental subroutine initialize_order_s(self)
  !< Initialize integrator for s-th order formula.
  class(integrator_runge_kutta_lssp), intent(inout) :: self !< Integrator.
  integer(I_P)                                      :: i    !< Counter.
  integer(I_P)                                      :: j    !< Counter.

  associate(alpha=>self%alpha)
    alpha(1) = 1._R_P
    do i=2, self%stages
      alpha(i) = 1._R_P / product([(j, j=1,i)])
      do j=i-1, 2, -1
        alpha(j) = 1._R_P / (j-1) * alpha(j-1)
      enddo
      alpha(1) = 1._R_P
      do j=2, i
        alpha(1) = alpha(1) - alpha(j)
      enddo
    enddo
  endassociate

  endsubroutine initialize_order_s

  subroutine integrate_order_s_1(self, U, stage, Dt, t)
  !< Integrate integrand field by (s-1)-th order formula.
  !<
  !< @note This method can be used **after** the integrator is created (i.e. the RK coefficients are initialized).
  class(integrator_runge_kutta_lssp), intent(in)    :: self      !< Integrator.
  class(integrand_object),            intent(inout) :: U         !< Field to be integrated.
  class(integrand_object),            intent(inout) :: stage(1:) !< Runge-Kutta stages [1:stages].
  real(R_P),                          intent(in)    :: Dt        !< Time step.
  real(R_P),                          intent(in)    :: t         !< Time.
  integer(I_P)                                      :: s         !< First stages counter.

  ! computing stages
  stage(1) = U
  do s=2, self%stages
    stage(s) = stage(s-1) + (stage(s-1)%t(t=t) * (Dt * 0.5_R_P))
  enddo
  stage(self%stages) = stage(self%stages) + (stage(self%stages)%t(t=t) * (Dt * 0.5_R_P))
  ! computing new time step
  U = U * 0._R_P
  do s=1, self%stages
    U = U + (stage(s) * self%alpha(s))
  enddo
  endsubroutine integrate_order_s_1

  subroutine integrate_order_s(self, U, stage, Dt, t)
  !< Integrate integrand field by s-th order formula.
  !<
  !< @note This method can be used **after** the integrator is created (i.e. the RK coefficients are initialized).
  class(integrator_runge_kutta_lssp), intent(in)    :: self      !< Integrator.
  class(integrand_object),            intent(inout) :: U         !< Field to be integrated.
  class(integrand_object),            intent(inout) :: stage(1:) !< Runge-Kutta stages [1:stages].
  real(R_P),                          intent(in)    :: Dt        !< Time step.
  real(R_P),                          intent(in)    :: t         !< Time.
  integer(I_P)                                      :: s         !< First stages counter.

  ! computing stages
  stage(1) = U
  do s=2, self%stages
    stage(s) = stage(s-1) + (stage(s-1)%t(t=t) * Dt)
  enddo
  stage(self%stages) = stage(self%stages) + (stage(self%stages)%t(t=t) * Dt)
  ! computing new time step
  U = U * 0._R_P
  do s=1, self%stages
    U = U + (stage(s) * self%alpha(s))
  enddo
  endsubroutine integrate_order_s
endmodule foodie_integrator_runge_kutta_lssp
