!< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property, from
!< 2nd to 3rd order accurate.

module foodie_integrator_lmm_ssp
!< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property, from
!< 2nd to 3rd order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the LMM-SSP class scheme implemented is:
!<
!< $$ U^{n+N_s} = \sum_{s=1}^{N_s}{\left[a_s U^{n+s-1} + \Delta t b_s \cdot R(t^{n+s-1}, U^{n+s-1}) \right]} $$
!<
!<where \(N_s\) is the number of previous steps considered.
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit. The coefficients *a,b* define the actual scheme, that is selected accordingly to the number of
!< **steps** used. Currently, the schemes provided have steps number in *[3, 5]*. The formal order of accuracy varies
!< consistently in *[2nd, 3rd]* order.
!<
!<#### Bibliography
!< [1] *Strong Stability Preserving Runge-Kutta and Multistep Time Discretizations*, S. Gottlieb, D. Ketcheson, C.W. Shu,
!< 2011, 978-981-4289-26-9, doi:10.1142/7498, World Scientific Publishing Co. Pte. Ltd.

use foodie_adt_integrand, only : integrand
use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_lmm_ssp

character(len=99), parameter :: class_name_='lmm_ssp'                                           !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:3)=[trim(class_name_)//'_steps_3_order_2', &
                                                         trim(class_name_)//'_steps_4_order_3', &
                                                         trim(class_name_)//'_steps_5_order_3'] !< List of supported schemes.

type, extends(integrator_object) :: integrator_lmm_ssp
  !< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property,
  !< from 2nd to 3rd order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *a,b* coefficients) before used.
  private
  integer(I_P), public   :: steps=0 !< Number of time steps.
  real(R_P), allocatable :: a(:)    !< *a* coefficients.
  real(R_P), allocatable :: b(:)    !< *b* coefficients.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy         !< Destroy the integrator.
    procedure, pass(self) :: initialize      !< Initialize (create) the integrator.
    procedure, pass(self) :: integrate       !< Integrate integrand field.
    procedure, pass(self) :: update_previous !< Cyclic update previous time steps.
endtype integrator_lmm_ssp

contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_lmm_ssp), intent(in) :: self       !< Integrator.
  character(len=99)                     :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_lmm_ssp), intent(in)           :: self             !< Integrator.
  character(*),              intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                   :: desc             !< Description.
  character(len=:), allocatable                   :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                     :: NL=new_line('a') !< New line character.
  integer(I_P)                                    :: s                !< Counter.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = ''
  desc = desc//prefix_//'Strong Stability preserving Linear-Multistep-Methods class'//NL
  desc = desc//prefix_//'  Supported schemes:'//NL
  do s=lbound(supported_schemes_, dim=1), ubound(supported_schemes_, dim=1) - 1
    desc = desc//prefix_//'    + '//supported_schemes_(s)//NL
  enddo
  desc = desc//prefix_//'    + '//supported_schemes_(ubound(supported_schemes_, dim=1))
  endfunction description

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_lmm_ssp), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),  intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  select type(rhs)
  class is (integrator_lmm_ssp)
                          lhs%steps = rhs%steps
    if (allocated(rhs%a)) lhs%a     = rhs%a
    if (allocated(rhs%b)) lhs%b     = rhs%b
  endselect
  endsubroutine integr_assign_integr

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_lmm_ssp), intent(in) :: self         !< Integrator.
  character(*),              intent(in) :: scheme       !< Selected scheme.
  logical                               :: is_supported !< Inquire result.
  integer(I_P)                          :: s            !< Counter.

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
  class(integrator_lmm_ssp), intent(in) :: self       !< Integrator.
  character(len=99), allocatable        :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_lmm_ssp), intent(inout) :: self !< Integrator.

  call self%destroy_abstract
  self%steps = 0
  if (allocated(self%a)) deallocate(self%a)
  if (allocated(self%b)) deallocate(self%b)
  endsubroutine destroy

  subroutine initialize(self, scheme)
  !< Create the actual LMM-SSP integrator: initialize the *a,b* coefficients.
  !<
  !< @note If the integrator is initialized with a bad (unsupported) number of required time steps the initialization fails and
  !< the integrator error status is updated consistently for external-provided errors handling.
  class(integrator_lmm_ssp), intent(inout) :: self   !< Integrator.
  character(*),              intent(in)    :: scheme !< Selected scheme.

  if (self%is_supported(scheme=scheme)) then
    call self%destroy
    select case(trim(adjustl(scheme)))
    case('lmm_ssp_steps_3_order_2')
      self%steps = 3
      allocate(self%a(1:self%steps)) ; self%a = 0.0_R_P
      allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%a(1) = 1._R_P/4._R_P
      self%a(2) = 0._R_P
      self%a(3) = 3._R_P/4._R_P

      self%b(1) = 0._R_P
      self%b(2) = 0._R_P
      self%b(3) = 3._R_P/2._R_P
    case('lmm_ssp_steps_4_order_3')
      self%steps = 4
      allocate(self%a(1:self%steps)) ; self%a = 0.0_R_P
      allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%a(1) = 11._R_P/27._R_P
      self%a(2) = 0._R_P
      self%a(3) = 0._R_P
      self%a(4) = 16._R_P/27._R_P

      self%b(1) = 12._R_P/27._R_P
      self%b(2) = 0._R_P
      self%b(3) = 0._R_P
      self%b(4) = 16._R_P/9._R_P
    case('lmm_ssp_steps_5_order_3')
      self%steps = 5
      allocate(self%a(1:self%steps)) ; self%a = 0.0_R_P
      allocate(self%b(1:self%steps)) ; self%b = 0.0_R_P
      self%a(1) = 7._R_P/32._R_P
      self%a(2) = 0._R_P
      self%a(3) = 0._R_P
      self%a(4) = 0._R_P
      self%a(5) = 25._R_P/32._R_P

      self%b(1) = 5._R_P/16._R_P
      self%b(2) = 0._R_P
      self%b(3) = 0._R_P
      self%b(4) = 0._R_P
      self%b(5) = 25._R_P/16._R_P
    endselect
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=.true.)
  endif
  endsubroutine initialize

  subroutine integrate(self, U, previous, Dt, t, autoupdate)
  !< Integrate field with LMM-SSP class scheme.
  class(integrator_lmm_ssp), intent(in)    :: self         !< Integrator.
  class(integrand),          intent(inout) :: U            !< Field to be integrated.
  class(integrand),          intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                 intent(in)    :: Dt           !< Time steps.
  real(R_P),                 intent(in)    :: t(:)         !< Times.
  logical, optional,         intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  logical                                  :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  integer(I_P)                             :: s            !< Steps counter.

  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  U = U * 0._R_P
  do s=1, self%steps
    if (self%a(s) /= 0._R_P) U = U + previous(s) * self%a(s)
    if (self%b(s) /= 0._R_P) U = U + previous(s)%t(t=t(s)) * (Dt * self%b(s))
  enddo
  if (autoupdate_) call self%update_previous(U=U, previous=previous)
  endsubroutine integrate

  subroutine update_previous(self, U, previous)
  !< Cyclic update previous time steps.
  class(integrator_lmm_ssp), intent(in)    :: self         !< Integrator.
  class(integrand),          intent(in)    :: U            !< Field to be integrated.
  class(integrand),          intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  integer(I_P)                             :: s            !< Steps counter.

  do s=1, self%steps - 1
    previous(s) = previous(s + 1)
  enddo
  previous(self%steps) = U
  endsubroutine update_previous
endmodule foodie_integrator_lmm_ssp
