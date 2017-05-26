!< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property and
!< variable stepsize (VSS), from 2nd to 3rd order accurate.

module foodie_integrator_lmm_ssp_vss
!< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property and
!< variable stepsize (VSS), from 2nd to 3rd order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the LMM-SSP class scheme implemented is:
!<
!<#### Second order formula
!<
!< $$ U^{n+N_s} = \frac{1}{\Omega_{N_s-1}^2} U^n + \frac{\Omega_{N_s-1}^2 - 1}{\Omega_{N_s-1}^2} U^{n+N_s-1} +
!<    \frac{\Omega_{N_s-1} + 1}{\Omega_{N_s-1}} \Delta t^{n+N_s} R(U^{n+N_s-1}) $$
!<
!<#### Third order formula
!<
!< $$ U^{n+N_s} = \frac{3 \Omega_{N_s-1} + 2}{\Omega_{N_s-1}^3} U^n +
!<                \frac{(\Omega_{N_s-1} + 1)^2(\Omega_{N_s-1} - 2)}{\Omega_{N_s-1}^3} U^{n+N_s-1} +
!<                \frac{\Omega_{N_s-1} + 1}{\Omega_{N_s-1}^2} \Delta t^{n+N_s} R(U^n) +
!<                \frac{(\Omega_{N_s-1} + 1)^2}{\Omega_{N_s-1}^2} \Delta t^{n+N_s} R(U^{n+N_s-1}) $$
!<
!<where \(N_s\) is the number of previous steps considered and
!<
!< $$ \Omega_s = \sum_{i=1}^s { \omega_i }\quad 1 \leq s \leq N_s $$
!< $$ \omega_i = \frac{\Delta t^{n + s}}{\Delta t^{n + N_s}} $$
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit.
!<
!<#### Bibliography
!< [1] *Strong Stability Preserving Explicit Linear Multistep Methods with Variable Step Size*, Y. Hadjmichael, D. Ketcheson,
!< L. Loczi, A. Nemeth, 2016, SIAM, Vol. 54, N. 5, pp. 2799-2832.

use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrand_object, only : integrand_object
use foodie_integrator_multistep_object, only : integrator_multistep_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_lmm_ssp_vss

character(len=99), parameter :: class_name_='lmm_ssp_vss'                                       !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:5)=[trim(class_name_)//'_steps_2_order_2', &
                                                         trim(class_name_)//'_steps_3_order_2', &
                                                         trim(class_name_)//'_steps_3_order_3', &
                                                         trim(class_name_)//'_steps_4_order_3', &
                                                         trim(class_name_)//'_steps_5_order_3'] !< List of supported schemes.

logical, parameter :: has_fast_mode_=.true.  !< Flag to check if integrator provides *fast mode* integrate.

type, extends(integrator_multistep_object) :: integrator_lmm_ssp_vss
  !< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property and
  !< variable stepsize (VSS), from 2nd to 3rd order accurate.
  !<
  !< @note The integrator must be created or initialized before used.
  private
  procedure(integrate_interface),      pointer :: integrate_ => integrate_order_2           !< Integrate integrand field.
  procedure(integrate_fast_interface), pointer :: integrate_fast_ => integrate_order_2_fast !< Integrate integrand field, fast.
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
    ! private methods
    procedure, pass(self), private :: integrate_order_2      !< Integrate integrand field by 2nd order formula.
    procedure, pass(self), private :: integrate_order_3      !< Integrate integrand field by 3rd order formula.
    procedure, pass(self), private :: integrate_order_2_fast !< Integrate integrand field by 2nd order formula, fast mode.
    procedure, pass(self), private :: integrate_order_3_fast !< Integrate integrand field by 3rd order formula, fast mode.
endtype integrator_lmm_ssp_vss

abstract interface
  !< Abstract interfaces of [[integrator_lmm_ssp_vss]] methods.
  subroutine integrate_interface(self, U, previous, Dt, t)
  !< Integrate field with LMM-SSP class scheme.
  import :: integrand_object, integrator_lmm_ssp_vss, R_P
  class(integrator_lmm_ssp_vss), intent(inout) :: self         !< Integrator.
  class(integrand_object),       intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),       intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(in)    :: Dt           !< Time step.
  real(R_P),                     intent(in)    :: t            !< Time.
  endsubroutine integrate_interface

  subroutine integrate_fast_interface(self, U, previous, buffer, Dt, t)
  !< Integrate field with LMM-SSP class scheme.
  import :: integrand_object, integrator_lmm_ssp_vss, R_P
  class(integrator_lmm_ssp_vss), intent(inout) :: self         !< Integrator.
  class(integrand_object),       intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),       intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  class(integrand_object),       intent(inout) :: buffer       !< Temporary buffer for doing fast operation.
  real(R_P),                     intent(in)    :: Dt           !< Time step.
  real(R_P),                     intent(in)    :: t            !< Time.
  endsubroutine integrate_fast_interface
endinterface

contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_lmm_ssp_vss), intent(in) :: self       !< Integrator.
  character(len=99)                         :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_lmm_ssp_vss), intent(in)           :: self             !< Integrator.
  character(*),                  intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                       :: desc             !< Description.
  character(len=:), allocatable                       :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                         :: NL=new_line('a') !< New line character.
  integer(I_P)                                        :: s                !< Counter.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = ''
  desc = desc//prefix_//'Strong Stability preserving Linear-Multistep-Methods Variable Stepsize class'//NL
  desc = desc//prefix_//'  Supported schemes:'//NL
  do s=lbound(supported_schemes_, dim=1), ubound(supported_schemes_, dim=1) - 1
    desc = desc//prefix_//'    + '//supported_schemes_(s)//NL
  enddo
  desc = desc//prefix_//'    + '//supported_schemes_(ubound(supported_schemes_, dim=1))
  endfunction description

  elemental function has_fast_mode(self)
  !< Return .true. if the integrator class has *fast mode* integrate.
  class(integrator_lmm_ssp_vss), intent(in) :: self          !< Integrator.
  logical                                   :: has_fast_mode !< Inquire result.

  has_fast_mode = has_fast_mode_
  endfunction has_fast_mode

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_lmm_ssp_vss), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),      intent(in)    :: rhs !< Right hand side.

  call lhs%assign_multistep(rhs=rhs)
  select type(rhs)
  class is (integrator_lmm_ssp_vss)
    if (associated(rhs%integrate_)) lhs%integrate_ => rhs%integrate_
  endselect
  endsubroutine integr_assign_integr

  subroutine integrate(self, U, Dt, t)
  !< Integrate field with LMM-SSP class scheme.
  class(integrator_lmm_ssp_vss), intent(inout) :: self !< Integrator.
  class(integrand_object),       intent(inout) :: U    !< Field to be integrated.
  real(R_P),                     intent(in)    :: Dt   !< Time step.
  real(R_P),                     intent(in)    :: t    !< Time.

  call self%integrate_(U=U, previous=self%previous, Dt=Dt, t=t)
  endsubroutine integrate

  subroutine integrate_fast(self, U, Dt, t)
  !< Integrate field with LMM-SSP class scheme, fast mode.
  class(integrator_lmm_ssp_vss), intent(inout) :: self !< Integrator.
  class(integrand_object),       intent(inout) :: U    !< Field to be integrated.
  real(R_P),                     intent(in)    :: Dt   !< Time step.
  real(R_P),                     intent(in)    :: t    !< Time.

  call self%integrate_fast_(U=U, previous=self%previous, buffer=self%buffer, Dt=Dt, t=t)
  endsubroutine integrate_fast

  subroutine integrate_ub(self, U, previous, Dt, t)
  !< Integrate field with LMM-SSP class scheme, unbuffered.
  class(integrator_lmm_ssp_vss), intent(inout) :: self         !< Integrator.
  class(integrand_object),       intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),       intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(in)    :: Dt           !< Time step.
  real(R_P),                     intent(in)    :: t            !< Time.

  call self%integrate_(U=U, previous=previous, Dt=Dt, t=t)
  endsubroutine integrate_ub

  subroutine integrate_ub_fast(self, U, previous, Dt, t)
  !< Integrate field with LMM-SSP class scheme, unbuffered, fast mode.
  class(integrator_lmm_ssp_vss), intent(inout) :: self         !< Integrator.
  class(integrand_object),       intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),       intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(in)    :: Dt           !< Time step.
  real(R_P),                     intent(in)    :: t            !< Time.

  call self%integrate_fast_(U=U, previous=previous, buffer=self%buffer, Dt=Dt, t=t)
  endsubroutine integrate_ub_fast

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_lmm_ssp_vss), intent(in) :: self         !< Integrator.
  character(*),                  intent(in) :: scheme       !< Selected scheme.
  logical                                   :: is_supported !< Inquire result.
  integer(I_P)                              :: s            !< Counter.

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
  class(integrator_lmm_ssp_vss), intent(in) :: self       !< Integrator.
  character(len=99), allocatable            :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_lmm_ssp_vss), intent(inout) :: self !< Integrator.

  call self%destroy_multistep
  self%integrate_ => integrate_order_2
  endsubroutine destroy

  subroutine initialize(self, scheme, autoupdate, U, stop_on_fail)
  !< Create the actual LMM-SSP-VSS integrator.
  !<
  !< @note If the integrator is initialized with a bad (unsupported) number of required time steps the initialization fails and
  !< the integrator error status is updated consistently for external-provided errors handling.
  class(integrator_lmm_ssp_vss), intent(inout)        :: self         !< Integrator.
  character(*),                  intent(in)           :: scheme       !< Selected scheme.
  logical,                       intent(in), optional :: autoupdate   !< Enable cyclic autoupdate of previous time steps.
  class(integrand_object),       intent(in), optional :: U            !< Integrand molding prototype.
  logical,                       intent(in), optional :: stop_on_fail !< Stop execution if initialization fail.

  if (self%is_supported(scheme=scheme)) then
    call self%destroy
    select case(trim(adjustl(scheme)))
    case('lmm_ssp_vss_steps_2_order_2')
      self%steps = 2
      self%integrate_ => integrate_order_2
      self%integrate_fast_ => integrate_order_2_fast
    case('lmm_ssp_vss_steps_3_order_2')
      self%steps = 3
      self%integrate_ => integrate_order_2
      self%integrate_fast_ => integrate_order_2_fast
    case('lmm_ssp_vss_steps_3_order_3')
      self%steps = 3
      self%integrate_ => integrate_order_3
      self%integrate_fast_ => integrate_order_3_fast
    case('lmm_ssp_vss_steps_4_order_3')
      self%steps = 4
      self%integrate_ => integrate_order_3
      self%integrate_fast_ => integrate_order_3_fast
    case('lmm_ssp_vss_steps_5_order_3')
      self%steps = 5
      self%integrate_ => integrate_order_3
      self%integrate_fast_ => integrate_order_3_fast
    endselect
    self%autoupdate = .true. ; if (present(autoupdate)) self%autoupdate = autoupdate
    self%registers = self%steps
    if (present(U)) call self%allocate_integrand_members(U=U)
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=stop_on_fail)
  endif
  endsubroutine initialize

  ! private methods
  subroutine integrate_order_2(self, U, previous, Dt, t)
  !< Integrate field with LMM-SSP-VSS 2nd order class scheme.
  class(integrator_lmm_ssp_vss), intent(inout) :: self         !< Integrator.
  class(integrand_object),       intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),       intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(in)    :: Dt           !< Time step.
  real(R_P),                     intent(in)    :: t            !< Time.
  real(R_P)                                    :: omega_       !< Omega coefficient.
  real(R_P)                                    :: omega_sq     !< Square of omega coefficient.

  omega_ = omega(Dt=self%Dt, s=self%steps-1)
  omega_sq = omega_ * omega_
  U = (previous(1) * (1._R_P / omega_sq)) + (previous(self%steps) * ((omega_sq - 1._R_P) / omega_sq)) + &
      (previous(self%steps)%t(t=self%t(self%steps)) * (self%Dt(self%steps) * (omega_ + 1._R_P) / omega_))
  if (self%autoupdate) call self%update_previous(U=U, previous=previous, Dt=Dt, t=t, previous_Dt=self%Dt, previous_t=self%t)
  endsubroutine integrate_order_2

  subroutine integrate_order_3(self, U, previous, Dt, t)
  !< Integrate field with LMM-SSP-VSS 3rd order class scheme.
  class(integrator_lmm_ssp_vss), intent(inout) :: self         !< Integrator.
  class(integrand_object),       intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),       intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(in)    :: Dt           !< Time step.
  real(R_P),                     intent(in)    :: t            !< Time.
  real(R_P)                                    :: omega_       !< Omega coefficient.

  omega_= omega(Dt=self%Dt, s=self%steps-1)
  U = (previous(1) * ((3._R_P * omega_ + 2._R_P) / omega_ ** 3 )) +                           &
      (previous(self%steps) * (((omega_ + 1._R_P) ** 2) * (omega_ - 2._R_P) / omega_ ** 3)) + &
      (previous(1)%t(t=self%t(1)) * (self%Dt(self%steps) * (omega_ + 1._R_P) / omega_ ** 2)) +          &
      (previous(self%steps)%t(t=self%t(self%steps)) * (self%Dt(self%steps) * (omega_ + 1._R_P) ** 2 / omega_ ** 2))
  if (self%autoupdate) call self%update_previous(U=U, previous=previous, Dt=Dt, t=t, previous_Dt=self%Dt, previous_t=self%t)
  endsubroutine integrate_order_3

  subroutine integrate_order_2_fast(self, U, previous, buffer, Dt, t)
  !< Integrate field with LMM-SSP-VSS 2nd order class scheme.
  class(integrator_lmm_ssp_vss), intent(inout) :: self         !< Integrator.
  class(integrand_object),       intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),       intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  class(integrand_object),       intent(inout) :: buffer       !< Temporary buffer for doing fast operation.
  real(R_P),                     intent(in)    :: Dt           !< Time step.
  real(R_P),                     intent(in)    :: t            !< Time.
  real(R_P)                                    :: omega_       !< Omega coefficient.
  real(R_P)                                    :: omega_sq     !< Square of omega coefficient.

  omega_ = omega(Dt=self%Dt, s=self%steps-1)
  omega_sq = omega_ * omega_
  call U%multiply_fast(lhs=previous(1), rhs=1._R_P / omega_sq)
  call buffer%multiply_fast(lhs=previous(self%steps), rhs=(omega_sq - 1._R_P) / omega_sq)
  call U%add_fast(lhs=U, rhs=buffer)
  buffer = previous(self%steps)
  call buffer%t_fast(t=self%t(self%steps))
  call buffer%multiply_fast(lhs=buffer, rhs=self%Dt(self%steps) * (omega_ + 1._R_P) / omega_)
  call U%add_fast(lhs=U, rhs=buffer)
  if (self%autoupdate) call self%update_previous(U=U, previous=previous, Dt=Dt, t=t, previous_Dt=self%Dt, previous_t=self%t)
  endsubroutine integrate_order_2_fast

  subroutine integrate_order_3_fast(self, U, previous, buffer, Dt, t)
  !< Integrate field with LMM-SSP-VSS 3rd order class scheme, fast mode.
  class(integrator_lmm_ssp_vss), intent(inout) :: self         !< Integrator.
  class(integrand_object),       intent(inout) :: U            !< Field to be integrated.
  class(integrand_object),       intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  class(integrand_object),       intent(inout) :: buffer       !< Temporary buffer for doing fast operation.
  real(R_P),                     intent(in)    :: Dt           !< Time step.
  real(R_P),                     intent(in)    :: t            !< Time.
  real(R_P)                                    :: omega_       !< Omega coefficient.

  omega_= omega(Dt=self%Dt, s=self%steps-1)

  call U%multiply_fast(lhs=previous(1), rhs=(3._R_P * omega_ + 2._R_P) / (omega_ ** 3))
  call buffer%multiply_fast(lhs=previous(self%steps), rhs=(((omega_ + 1._R_P) ** 2) * (omega_ - 2._R_P) / (omega_ ** 3)))
  call U%add_fast(lhs=U, rhs=buffer)
  buffer = previous(1)
  call buffer%t_fast(t=self%t(1))
  call buffer%multiply_fast(lhs=buffer, rhs=self%Dt(self%steps) * (omega_ + 1._R_P) / (omega_ ** 2))
  call U%add_fast(lhs=U, rhs=buffer)
  buffer = previous(self%steps)
  call buffer%t_fast(t=self%t(self%steps))
  call buffer%multiply_fast(lhs=buffer, rhs=(self%Dt(self%steps) * (omega_ + 1._R_P) ** 2 / (omega_ ** 2)))
  call U%add_fast(lhs=U, rhs=buffer)
  if (self%autoupdate) call self%update_previous(U=U, previous=previous, Dt=Dt, t=t, previous_Dt=self%Dt, previous_t=self%t)
  endsubroutine integrate_order_3_fast

  ! private non TBP
  pure function dt_ratio(Dt, s) result(ratio)
  !< Return `Dt(n+s)/Dt(n+Ns)` ratio.
  real(R_P),    intent(in) :: Dt(:) !< Time steps.
  integer(I_P), intent(in) :: s     !< Step index.
  real(R_P)                :: ratio !< Time steps ratio.

  ratio = Dt(s)/Dt(ubound(Dt, dim=1))
  endfunction dt_ratio

  pure function omega(Dt, s)
  !< Return `omega=sum(dt_ratio(i)), i=1, s`.
  real(R_P),    intent(in) :: Dt(:) !< Time steps.
  integer(I_P), intent(in) :: s     !< Step index.
  real(R_P)                :: omega !< Omega sum.
  integer(I_P)             :: i     !< Counter.

  omega = 0._R_P
  do i=1, s
    omega = omega + dt_Ratio(Dt=Dt, s=i)
  enddo
  endfunction omega
endmodule foodie_integrator_lmm_ssp_vss
