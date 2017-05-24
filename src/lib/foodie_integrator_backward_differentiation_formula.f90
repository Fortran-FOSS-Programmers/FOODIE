!< FOODIE integrator: provide an implicit class of Backward Differentiation Formula schemes, from 1st to 6th order accurate.

module foodie_integrator_backward_differentiation_formula
!< FOODIE integrator: provide an implicit class of Backward Differentiation Formula schemes, from 1st to 6th order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the Backward Differentiation Formula class scheme implemented is:
!<
!< $$ U^{n+N_s} + \sum_{s=1}^{N_s} \alpha_s U^{n+N_s-s} = \Delta t \left[ \beta R(t^{n+N_s}, U^{n+N_s}) \right] $$
!<
!< where \(N_s\) is the number of previous steps considered.
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are implicit. The coefficients \(\alpha_s\) and \(\beta\) define the actual scheme, that is selected accordingly
!< to the number of *steps* used.
!<
!< Currently, the following schemes are available:
!<
!< | ` Step ` | ` beta `   | ` alpha 1 `  | ` alpha 2 ` | ` alpha 3 `  | ` alpha 4 ` | ` alpha 5 ` | ` alpha 6 ` |
!< |----------|------------|--------------|-------------|--------------|-------------|-------------|-------------|
!< | ` 1 `    | `  1 `     | ` -1 `       |             |              |             |             |             |
!< | ` 2 `    | ` 2/3 `    | ` -4/3 `     | ` 1/3 `     |              |             |             |             |
!< | ` 3 `    | ` 6/11 `   | ` -18/11 `   | ` 9/11  `   | ` -2/11 `    |             |             |             |
!< | ` 4 `    | ` 12/25 `  | ` -48/25 `   | ` 36/25 `   | ` -16/25 `   | ` 3/25 `    |             |             |
!< | ` 5 `    | ` 60/137 ` | ` -300/137 ` | ` 300/137 ` | ` -200/137 ` | ` 75/137 `  | ` -12/137 ` |             |
!< | ` 6 `    | ` 60/147 ` | ` -360/147 ` | ` 450/147 ` | ` -400/147 ` | ` 225/147 ` | ` -72/147 ` | ` 10/147 `  |
!<
!<#### Bibliography
!<

use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrand_object, only : integrand_object
use foodie_integrator_multistep_implicit_object, only : integrator_multistep_implicit_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_back_df

character(len=99), parameter :: class_name_='back_df'                             !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:6)=[trim(class_name_)//'_1', &
                                                         trim(class_name_)//'_2', &
                                                         trim(class_name_)//'_3', &
                                                         trim(class_name_)//'_4', &
                                                         trim(class_name_)//'_5', &
                                                         trim(class_name_)//'_6'] !< List of supported schemes.

logical, parameter :: has_fast_mode_=.true.  !< Flag to check if integrator provides *fast mode* integrate.

type, extends(integrator_multistep_implicit_object) :: integrator_back_df
  !< FOODIE integrator: provide an implicit class of Backward-Differentiation-Formula multi-step schemes, from 1st to 6th order
  !< accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *alpha* and *beta* coefficients) before used.
  private
  real(R_P), allocatable :: a(:)      !< \(\alpha\) coefficients.
  real(R_P)              :: b=0.0_R_P !< \(\beta\) coefficient.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(self) :: has_fast_mode        !< Return .true. if the integrator class has *fast mode* integrate.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: integrate            !< Integrate integrand field.
    procedure, pass(self) :: integrate_fast       !< Integrate integrand field, fast mode.
    procedure, pass(self) :: integrate_ub         !< Integrate integrand field, unbuffered.
    procedure, pass(self) :: integrate_ub_fast    !< Integrate integrand field, unbuffered, fast mode.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy    !< Destroy the integrator.
    procedure, pass(self) :: initialize !< Initialize (create) the integrator.
endtype integrator_back_df

contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_back_df), intent(in) :: self       !< Integrator.
  character(len=99)                     :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_back_df), intent(in)           :: self             !< Integrator.
  character(*),              intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                   :: desc             !< Description.
  character(len=:), allocatable                   :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                     :: NL=new_line('a') !< New line character.
  integer(I_P)                                    :: s                !< Counter.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = ''
  desc = desc//prefix_//'Backward-Differentiation-Formula multi-step schemes class'//NL
  desc = desc//prefix_//'  Supported schemes:'//NL
  do s=lbound(supported_schemes_, dim=1), ubound(supported_schemes_, dim=1) - 1
    desc = desc//prefix_//'    + '//supported_schemes_(s)//NL
  enddo
  desc = desc//prefix_//'    + '//supported_schemes_(ubound(supported_schemes_, dim=1))
  endfunction description

  elemental function has_fast_mode(self)
  !< Return .true. if the integrator class has *fast mode* integrate.
  class(integrator_back_df), intent(in) :: self          !< Integrator.
  logical                               :: has_fast_mode !< Inquire result.

  has_fast_mode = has_fast_mode_
  endfunction has_fast_mode

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_back_df), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),  intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  select type(rhs)
  class is (integrator_back_df)
                          lhs%steps = rhs%steps
    if (allocated(rhs%a)) lhs%a     = rhs%a
                          lhs%b     = rhs%b
  endselect
  endsubroutine integr_assign_integr

  subroutine integrate(self, U, Dt, t, iterations, autoupdate)
  !< Integrate field with BDF class scheme.
  !<
  !< @note This method uses integrand previous-steps-buffer stored inside integrator.
  class(integrator_back_df),  intent(inout) :: self        !< Integrator.
  class(integrand_object),    intent(inout) :: U           !< Field to be integrated.
  real(R_P),                  intent(in)    :: Dt          !< Time steps.
  real(R_P),                  intent(in)    :: t(:)        !< Times.
  integer(I_P), optional,     intent(in)    :: iterations  !< Fixed point iterations.
  logical,      optional,     intent(in)    :: autoupdate  !< Perform cyclic autoupdate of previous time steps.

  call self%integrate_ub(U=U, previous=self%previous, Dt=Dt, t=t, iterations=iterations, autoupdate=autoupdate)
  endsubroutine integrate

  subroutine integrate_fast(self, U, Dt, t, iterations, autoupdate)
  !< Integrate field with BDF class scheme.
  !<
  !< @note This method uses integrand previous-steps-buffer stored inside integrator.
  class(integrator_back_df),  intent(inout) :: self         !< Integrator.
  class(integrand_object),    intent(inout) :: U            !< Field to be integrated.
  real(R_P),                  intent(in)    :: Dt           !< Time steps.
  real(R_P),                  intent(in)    :: t(:)         !< Times.
  integer(I_P), optional,     intent(in)    :: iterations   !< Fixed point iterations.
  logical,      optional,     intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.

  call self%integrate_ub_fast(U=U, previous=self%previous, Dt=Dt, t=t, iterations=iterations, autoupdate=autoupdate)
  endsubroutine integrate_fast

  subroutine integrate_ub(self, U, previous, Dt, t, iterations, autoupdate)
  !< Integrate field with BDF class scheme.
  class(integrator_back_df),  intent(in)    :: self         !< Integrator.
  class(integrand_object),    intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),    intent(inout) :: previous(1:) !< Integrand.
  real(R_P),                  intent(in)    :: Dt           !< Time steps.
  real(R_P),                  intent(in)    :: t(:)         !< Times.
  integer(I_P), optional,     intent(in)    :: iterations   !< Fixed point iterations.
  logical,      optional,     intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  integer(I_P)                              :: iterations_  !< Fixed point iterations.
  logical                                   :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  class(integrand_object), allocatable      :: delta        !< Delta RHS for fixed point iterations.
  integer(I_P)                              :: s            !< Steps counter.

  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  iterations_ = 1 ; if (present(iterations)) iterations_ = iterations
  allocate(delta, mold=U)
  delta = previous(self%steps) * (-self%a(self%steps))
  do s=1, self%steps - 1
    delta = delta + (previous(s) * (-self%a(s)))
  enddo
  do s=1, iterations_
    U = delta + (U%t(t=t(self%steps) + Dt) * (Dt * self%b))
  enddo
  if (autoupdate_) call self%update_previous(U=U, previous=previous, is_like_explicit=.true.)
  endsubroutine integrate_ub

  subroutine integrate_ub_fast(self, U, previous, Dt, t, iterations, autoupdate)
  !< Integrate field with BDF class scheme.
  class(integrator_back_df),  intent(inout) :: self         !< Integrator.
  class(integrand_object),    intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),    intent(inout) :: previous(1:) !< Integrand.
  real(R_P),                  intent(in)    :: Dt           !< Time steps.
  real(R_P),                  intent(in)    :: t(:)         !< Times.
  integer(I_P), optional,     intent(in)    :: iterations   !< Fixed point iterations.
  logical,      optional,     intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  integer(I_P)                              :: iterations_  !< Fixed point iterations.
  logical                                   :: autoupdate_  !< Perform cyclic autoupdate of previous time steps, dummy var.
  class(integrand_object), allocatable      :: delta        !< Delta RHS for fixed point iterations.
  integer(I_P)                              :: s            !< Steps counter.

  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  iterations_ = 1 ; if (present(iterations)) iterations_ = iterations
  allocate(delta, mold=U)
  call delta%multiply_fast(lhs=previous(self%steps), rhs=-self%a(self%steps))
  do s=1, self%steps - 1
    call self%buffer%multiply_fast(lhs=previous(s), rhs=-self%a(s))
    call delta%add_fast(lhs=delta, rhs=self%buffer)
  enddo
  do s=1, iterations
    self%buffer = U
    call self%buffer%t_fast(t=t(self%steps) + Dt)
    call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%b)
    call U%add_fast(lhs=delta, rhs=self%buffer)
  enddo
  if (autoupdate_) call self%update_previous(U=U, previous=previous, is_like_explicit=.true.)
  endsubroutine integrate_ub_fast

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_back_df), intent(in) :: self         !< Integrator.
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
  class(integrator_back_df), intent(in) :: self       !< Integrator.
  character(len=99), allocatable        :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_back_df), intent(inout) :: self !< Integrator.

  call self%destroy_multistep
  if (allocated(self%a)) deallocate(self%a)
  self%b = 0.0_R_P
  endsubroutine destroy

  subroutine initialize(self, scheme, U, stop_on_fail)
  !< Create the actual BDF integrator: initialize the *alpha* and *beta* coefficients.
  class(integrator_back_df), intent(inout)        :: self         !< Integrator.
  character(*),              intent(in)           :: scheme       !< Selected scheme.
  class(integrand_object),   intent(in), optional :: U            !< Integrand molding prototype.
  logical,                   intent(in), optional :: stop_on_fail !< Stop execution if initialization fail.

  if (self%is_supported(scheme=scheme)) then
    call self%destroy
    select case(trim(adjustl(scheme)))
    case('back_df_1')
      self%steps = 1 ; allocate(self%a(1:self%steps)) ; self%a = 0.0_R_P
      self%a(1) = -1.0_R_P
      self%b = 1.0_R_P
    case('back_df_2')
      self%steps = 2 ; allocate(self%a(1:self%steps)) ; self%a = 0.0_R_P
      self%a(1) = 1.0_R_P/3.0_R_P
      self%a(2) = -4.0_R_P/3.0_R_P
      self%b = 2.0_R_P/3.0_R_P
    case('back_df_3')
      self%steps = 3 ; allocate(self%a(1:self%steps)) ; self%a = 0.0_R_P
      self%a(1) = -2.0_R_P/11.0_R_P
      self%a(2) = 9.0_R_P/11.0_R_P
      self%a(3) = -18.0_R_P/11.0_R_P
      self%b = 6.0_R_P/11.0_R_P
    case('back_df_4')
      self%steps = 4 ; allocate(self%a(1:self%steps)) ; self%a = 0.0_R_P
      self%a(1) = 3.0_R_P/25.0_R_P
      self%a(2) = -16.0_R_P/25.0_R_P
      self%a(3) = 36.0_R_P/25.0_R_P
      self%a(4) = -48.0_R_P/25.0_R_P
      self%b = 12.0_R_P/25.0_R_P
    case('back_df_5')
      self%steps = 5 ; allocate(self%a(1:self%steps)) ; self%a = 0.0_R_P
      self%a(1) = -12.0_R_P/137.0_R_P
      self%a(2) = 75.0_R_P/137.0_R_P
      self%a(3) = -200.0_R_P/137.0_R_P
      self%a(4) = 300.0_R_P/137.0_R_P
      self%a(5) = -300.0_R_P/137.0_R_P
      self%b = 60.0_R_P/137.0_R_P
    case('back_df_6')
      self%steps = 6 ; allocate(self%a(1:self%steps)) ; self%a = 0.0_R_P
      self%a(1) = 10.0_R_P/147.0_R_P
      self%a(2) = -72.0_R_P/147.0_R_P
      self%a(3) = 225.0_R_P/147.0_R_P
      self%a(4) = -400.0_R_P/147.0_R_P
      self%a(5) = 450.0_R_P/147.0_R_P
      self%a(6) = -360.0_R_P/147.0_R_P
      self%b = 60.0_R_P/147.0_R_P
    case default
    endselect
    self%registers = self%steps
    if (present(U)) call self%allocate_integrand_members(U=U)
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=stop_on_fail)
  endif
  endsubroutine initialize
endmodule foodie_integrator_backward_differentiation_formula
