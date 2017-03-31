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

use foodie_adt_integrand, only : integrand
use foodie_error_codes, only : ERROR_INTEGRATOR_INIT_FAIL
use foodie_kinds, only : I_P, R_P
use foodie_integrator_object, only : integrator_object
use foodie_utils, only : is_admissible

implicit none
private
public :: integrator_lmm_ssp_vss

character(len=99), parameter :: supported_steps='2-5'        !< List of supported steps number. Valid format is `1-2,4,9-23...`.
character(len=99), parameter :: supported_orders='2-3'       !< List of supported accuracy orders. Valid format is `1-2,4,9-23...`.
integer(I_P),      parameter :: supported_steps_(1:2)=[2,5]  !< Supported steps range.
integer(I_P),      parameter :: supported_orders_(1:2)=[2,3] !< Supported orders range.

type, extends(integrator_object) :: integrator_lmm_ssp_vss
  !< FOODIE integrator: provide an explicit class of Linear Multi-step Methods (LLM) with Strong Stability Preserving property and
  !< variable stepsize (VSS), from 2nd to 3rd order accurate.
  !<
  !< @note The integrator must be created or initialized before used.
  private
  integer(I_P) :: steps=0 !< Number of time steps.
  integer(I_P) :: order=0 !< Order of accuracy.
  procedure(integrate_interface), pointer :: integrate_ => integrate_order_2 !< Integrate integrand field.
  contains
    ! deferred methods
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    ! public methods
    procedure, pass(self) :: destroy         !< Destroy the integrator.
    procedure, pass(self) :: init            !< Initialize (create) the integrator.
    procedure, pass(self) :: integrate       !< Integrate integrand field.
    procedure, nopass     :: is_supported    !< Check if the queried number of steps is supported or not.
    procedure, nopass     :: min_steps       !< Return the minimum number of steps supported.
    procedure, nopass     :: max_steps       !< Return the maximum number of steps supported.
    procedure, pass(self) :: update_previous !< Cyclic update previous time steps.
    ! private methods
    procedure, pass(self), private :: integrate_order_2 !< Integrate integrand field by 2nd order formula.
endtype integrator_lmm_ssp_vss

abstract interface
  !< Abstract interfaces of [[integrator_lmm_ssp_vss]] methods.
  subroutine integrate_interface(self, U, previous, Dt, t, autoupdate)
  !< Integrate field with LMM-SSP class scheme.
  import :: integrand, integrator_lmm_ssp_vss, R_P
  class(integrator_lmm_ssp_vss), intent(in)    :: self         !< Integrator.
  class(integrand),              intent(inout) :: U            !< Field to be integrated.
  class(integrand),              intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(inout) :: Dt(:)        !< Time steps.
  real(R_P),                     intent(in)    :: t(:)         !< Times.
  logical, optional,             intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
  endsubroutine integrate_interface
endinterface

contains
  ! deferred methods
  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_lmm_ssp_vss), intent(in)           :: self             !< Integrator.
  character(*),                  intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                       :: desc             !< Description.
  character(len=:), allocatable                       :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                         :: NL=new_line('a') !< New line character.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = ''
  desc = desc//prefix_//'Strong Stability preserving Linear-Multistep-Methods Variable Stepsize class'//NL
  desc = desc//prefix_//'  Supported steps numbers: ['//trim(adjustl(supported_steps))//']'//NL
  desc = desc//prefix_//'  Supported orders:        ['//trim(adjustl(supported_orders))//']'//NL
  desc = desc//prefix_//'  Not all steps/order combinations are supported'
  endfunction description

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_lmm_ssp_vss), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),      intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  select type(rhs)
  class is (integrator_lmm_ssp_vss)
                                    lhs%steps      =  rhs%steps
                                    lhs%order      =  rhs%order
    if (associated(rhs%integrate_)) lhs%integrate_ => rhs%integrate_
  endselect
  endsubroutine integr_assign_integr

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_lmm_ssp_vss), intent(inout) :: self !< Integrator.

  call self%destroy_abstract
  self%steps = 0
  self%order = 0
  self%integrate_ => integrate_order_2
  endsubroutine destroy

  subroutine init(self, steps, order, stop_on_fail)
  !< Create the actual LMM-SSP-VSS integrator.
  !<
  !< @note If the integrator is initialized with a bad (unsupported) number of required time steps the initialization fails and
  !< the integrator error status is updated consistently for external-provided errors handling.
  class(integrator_lmm_ssp_vss), intent(inout)        :: self         !< Integrator.
  integer(I_P),                  intent(in)           :: steps        !< Number of time steps used.
  integer(I_P),                  intent(in)           :: order        !< Order of accuracy.
  logical,                       intent(in), optional :: stop_on_fail !< Stop execution if initialization fail.

  if (self%is_supported(steps=steps, order=order)) then
    call self%destroy
    self%steps = steps
    self%order = order
    select case(order)
    case(2) ! LMM-SSP-VSS(steps,2)
      self%integrate_ => integrate_order_2
    case(3) ! LMM-SSP-VSS(steps,3)
      self%integrate_ => integrate_order_3
    endselect
  else
    call self%trigger_error(error=ERROR_INTEGRATOR_INIT_FAIL,                                            &
                            error_message='bad (unsupported) scheme'//new_line('a')//self%description(), &
                            is_severe=stop_on_fail)
  endif
  endsubroutine init

  subroutine integrate(self, U, previous, Dt, t, autoupdate)
  !< Integrate field with LMM-SSP class scheme.
  class(integrator_lmm_ssp_vss), intent(in)    :: self         !< Integrator.
  class(integrand),              intent(inout) :: U            !< Field to be integrated.
  class(integrand),              intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(inout) :: Dt(:)        !< Time steps.
  real(R_P),                     intent(in)    :: t(:)         !< Times.
  logical, optional,             intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.

  call self%integrate_(U=U, previous=previous, Dt=Dt, t=t, autoupdate=autoupdate)
  endsubroutine integrate

  elemental function is_supported(steps, order)
  !< Check if the queried number of steps is supported or not.
  integer(I_P), intent(in) :: steps        !< Number of time steps used.
  integer(I_P), intent(in) :: order        !< Order of accuracy.
  logical                  :: is_supported !< Is true if the steps number is in *supported_steps*.

  is_supported = is_admissible(n=steps, adm_range=trim(supported_steps))
  if (is_supported) is_supported = is_admissible(n=order, adm_range=trim(supported_orders))
  if (is_supported) then
     select case(order)
     case(2)
       is_supported = steps >= 2
     case(3)
       is_supported = steps >= 4
     endselect
  endif
  endfunction is_supported

  pure function min_steps()
  !< Return the minimum number of steps supported.
  integer(I_P) :: min_steps !< Minimum number of steps supported.

  min_steps = supported_steps_(1)
  endfunction min_steps

  pure function max_steps()
  !< Return the maximum number of steps supported.
  integer(I_P) :: max_steps !< Maximum number of steps supported.

  max_steps = supported_steps_(1)
  endfunction max_steps

  subroutine update_previous(self, U, previous, Dt)
  !< Cyclic update previous time steps.
  class(integrator_lmm_ssp_vss), intent(in)    :: self         !< Integrator.
  class(integrand),              intent(in)    :: U            !< Field to be integrated.
  class(integrand),              intent(inout) :: previous(1:) !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(inout) :: Dt(:)        !< Time steps.
  integer(I_P)                                 :: s            !< Steps counter.

  do s=1, self%steps - 1
    previous(s) = previous(s + 1)
    Dt(s) = Dt(s + 1)
  enddo
  previous(self%steps) = U
  endsubroutine update_previous

  ! private methods
  subroutine integrate_order_2(self, U, previous, Dt, t, autoupdate)
  !< Integrate field with LMM-SSP-VSS 2nd order class scheme.
  class(integrator_lmm_ssp_vss), intent(in)    :: self          !< Integrator.
  class(integrand),              intent(inout) :: U             !< Field to be integrated.
  class(integrand),              intent(inout) :: previous(1:)  !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(inout) :: Dt(:)         !< Time steps.
  real(R_P),                     intent(in)    :: t(:)          !< Times.
  logical, optional,             intent(in)    :: autoupdate    !< Perform cyclic autoupdate of previous time steps.
  logical                                      :: autoupdate_   !< Perform cyclic autoupdate of previous time steps, dummy var.
  integer(I_P)                                 :: s             !< Steps counter.
  real(R_P)                                    :: omega_        !< Omega coefficient.
  real(R_P)                                    :: omega_sq      !< Square of omega coefficient.

  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  omega_= omega(Dt=Dt, s=self%steps-1)
  omega_sq = omega_ * omega_
  U = previous(1) * (1._R_P / omega_sq) + previous(self%steps) * ((omega_sq - 1._R_P) / omega_sq) + &
      previous(self%steps)%t(t=t(self%steps)) * (Dt(self%steps) * (omega_ + 1._R_P) / omega_)
  if (autoupdate_) call self%update_previous(U=U, previous=previous, Dt=Dt)
  endsubroutine integrate_order_2

  subroutine integrate_order_3(self, U, previous, Dt, t, autoupdate)
  !< Integrate field with LMM-SSP-VSS 3rd order class scheme.
  class(integrator_lmm_ssp_vss), intent(in)    :: self          !< Integrator.
  class(integrand),              intent(inout) :: U             !< Field to be integrated.
  class(integrand),              intent(inout) :: previous(1:)  !< Previous time steps solutions of integrand field.
  real(R_P),                     intent(inout) :: Dt(:)         !< Time steps.
  real(R_P),                     intent(in)    :: t(:)          !< Times.
  logical, optional,             intent(in)    :: autoupdate    !< Perform cyclic autoupdate of previous time steps.
  logical                                      :: autoupdate_   !< Perform cyclic autoupdate of previous time steps, dummy var.
  integer(I_P)                                 :: s             !< Steps counter.
  real(R_P)                                    :: omega_        !< Omega coefficient.

  autoupdate_ = .true. ; if (present(autoupdate)) autoupdate_ = autoupdate
  omega_= omega(Dt=Dt, s=self%steps-1)
  U = previous(1) * ((3._R_P * omega_ + 2._R_P) / omega_ ** 3 ) +                           &
      previous(self%steps) * (((omega_ + 1._R_P) ** 2) * (omega_ - 2._R_P) / omega_ ** 3) + &
      previous(1)%t(t=t(1)) * (Dt(self%steps) * (omega_ + 1._R_P) / omega_ ** 2) +          &
      previous(self%steps)%t(t=t(self%steps)) * (Dt(self%steps) * (omega_ + 1._R_P) ** 2 / omega_ ** 2)
  if (autoupdate_) call self%update_previous(U=U, previous=previous, Dt=Dt)
  endsubroutine integrate_order_3

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
