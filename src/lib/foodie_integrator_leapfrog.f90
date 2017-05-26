!< FOODIE integrator: provide an explicit class of leapfrog multi-step schemes, 2nd order accurate.

module foodie_integrator_leapfrog
!< FOODIE integrator: provide an explicit class of leapfrog multi-step schemes, 2nd order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the leapfrog class scheme implemented (see [3]) is:
!<
!< $$ U^{n+2} = U^{n} + 2\Delta t \cdot R(t^{n+1}, U^{n+1}) $$
!<
!< Optionally, the Robert-Asselin-Williams (RAW) filter (see [3]) is applied to the computed integration steps:
!< $$ \Delta = \frac{\nu}{2}(U^{n} - 2 U^{n+1} + U^{n+2}) $$
!< $$ U^{n+1} = U^{n+1} + \Delta * \alpha $$
!< $$ U^{n+2} = U^{n+2} + \Delta * (\alpha-1) $$
!< Note that for \(\alpha=1\) the filter reverts back to the standard Robert-Asselin scheme.
!< The filter coefficients should be taken as \(\nu \in (0,1]\) and \(\alpha \in (0.5,1]\). The default values are
!<
!<  + \(\nu=0.01\)
!<  + \(\alpha=0.53\)
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit. The filter coefficients \(\nu,\,\alpha \) define the actual scheme.
!<
!<#### Bibliography
!<
!< [1] *The integration of a low order spectral form of the primitive meteorological equations*, Robert, A. J., J. Meteor. Soc.
!< Japan,vol. 44, pages 237--245, 1966.
!<
!< [2] *Frequency filter for time integrations*, Asselin, R., Monthly Weather Review, vol. 100, pages 487--490, 1972.
!<
!< [3] *The RAW filter: An improvement to the Robert–Asselin filter in semi-implicit integrations*, Williams, P.D., Monthly
!< Weather Review, vol. 139(6), pages 1996--2007, June 2011.

use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrand_object, only : integrand_object
use foodie_integrator_multistep_object, only : integrator_multistep_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_leapfrog

character(len=99), parameter :: class_name_='leapfrog'                              !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:2)=[trim(class_name_)//'    ', &
                                                         trim(class_name_)//'_raw'] !< List of supported schemes.

logical, parameter :: has_fast_mode_=.true.  !< Flag to check if integrator provides *fast mode* integrate.

type, extends(integrator_multistep_object) :: integrator_leapfrog
  !< FOODIE integrator: provide an explicit class of leapfrog multi-step schemes, 2nd order accurate.
  !<
  !< @note The integrator must be initialized before used.
  !<
  !< @note The time steps `Dt(1:steps)` passed to the integrate methods must be identical: this integrator supports only
  !< fixed time steps.
  private
  real(R_P)                            :: nu          !< Robert-Asselin filter coefficient.
  real(R_P)                            :: alpha       !< Robert-Asselin-Williams filter coefficient.
  logical                              :: is_filtered !< Flag to check if the integration if RAW filtered.
  class(integrand_object), allocatable :: filter      !< Filter field displacement.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(self) :: has_fast_mode        !< Return .true. if the integrator class has *fast mode* integrate.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: integrate            !< Integrate integrand field.
    procedure, pass(self) :: integrate_fast       !< Integrate integrand field, fast mode.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy    !< Destroy the integrator.
    procedure, pass(self) :: initialize !< Initialize (create) the integrator.
    ! overridden public methods
    procedure, pass(self) :: allocate_integrand_members !< Allocate integrand members.
endtype integrator_leapfrog

contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_leapfrog), intent(in) :: self       !< Integrator.
  character(len=99)                      :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_leapfrog), intent(in)           :: self    !< Integrator.
  character(*),               intent(in), optional :: prefix  !< Prefixing string.
  character(len=:), allocatable                    :: desc    !< Description.
  character(len=:), allocatable                    :: prefix_ !< Prefixing string, local variable.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = desc//prefix_//'Explicit leapfrog multi-step 2nd order scheme'
  endfunction description

  elemental function has_fast_mode(self)
  !< Return .true. if the integrator class has *fast mode* integrate.
  class(integrator_leapfrog), intent(in) :: self          !< Integrator.
  logical                                :: has_fast_mode !< Inquire result.

  has_fast_mode = has_fast_mode_
  endfunction has_fast_mode

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_leapfrog), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),   intent(in)    :: rhs !< Right hand side.

  call lhs%assign_multistep(rhs=rhs)
  select type(rhs)
  class is(integrator_leapfrog)
     lhs%nu = rhs%nu
     lhs%alpha = rhs%alpha
     lhs%is_filtered = rhs%is_filtered
     if (allocated(lhs%filter)) deallocate(lhs%filter)
     if (allocated(rhs%filter)) then
        allocate(lhs%filter, mold=rhs%filter)
        lhs%filter = rhs%filter
     endif
  endselect
  endsubroutine integr_assign_integr

  subroutine integrate(self, U, Dt, t)
  !< Integrate field with leapfrog class scheme.
  class(integrator_leapfrog), intent(inout) :: self !< Integrator.
  class(integrand_object),    intent(inout) :: U    !< Field to be integrated.
  real(R_P),                  intent(in)    :: Dt   !< Time step.
  real(R_P),                  intent(in)    :: t    !< Time.

  U = self%previous(1) + (self%previous(2)%t(t=t) * (Dt * 2._R_P))
  if (self%is_filtered) then
    self%filter = (self%previous(1) - (self%previous(2) * 2._R_P) + U) * self%nu * 0.5_R_P
    self%previous(2) = self%previous(2) + (self%filter * self%alpha)
    U = U + (self%filter * (self%alpha - 1._R_P))
  endif
  if (self%autoupdate) call self%update_previous(U=U, previous=self%previous)
  endsubroutine integrate

  subroutine integrate_fast(self, U, Dt, t)
  !< Integrate field with leapfrog class scheme, fast mode.
  class(integrator_leapfrog), intent(inout) :: self !< Integrator.
  class(integrand_object),    intent(inout) :: U    !< Field to be integrated.
  real(R_P),                  intent(in)    :: Dt   !< Time step.
  real(R_P),                  intent(in)    :: t    !< Time.

  self%buffer = self%previous(2)
  call self%buffer%t_fast(t=t)
  call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * 2._R_P)
  call U%add_fast(lhs=self%previous(1), rhs=self%buffer)
  if (self%is_filtered) then
    call self%buffer%multiply_fast(lhs=self%previous(2), rhs=2._R_P)
    call self%buffer%subtract_fast(lhs=self%previous(1), rhs=self%buffer)
    call self%buffer%add_fast(lhs=self%buffer, rhs=U)
    call self%filter%multiply_fast(lhs=self%buffer, rhs=self%nu * 0.5_R_P)

    call self%buffer%multiply_fast(lhs=self%filter, rhs=self%alpha)
    call self%previous(2)%add_fast(lhs=self%previous(2), rhs=self%buffer)

    call self%buffer%multiply_fast(lhs=self%filter, rhs=self%alpha - 1._R_P)
    call U%add_fast(lhs=U, rhs=self%buffer)
  endif
  if (self%autoupdate) call self%update_previous(U=U, previous=self%previous)
  endsubroutine integrate_fast

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_leapfrog), intent(in) :: self         !< Integrator.
  character(*),               intent(in) :: scheme       !< Selected scheme.
  logical                                :: is_supported !< Inquire result.
  integer(I_P)                           :: s            !< Counter.

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
  class(integrator_leapfrog), intent(in) :: self       !< Integrator.
  character(len=99), allocatable         :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_leapfrog), intent(INOUT) :: self !< Integrator.

  call self%destroy_multistep
  self%nu = 0._R_P
  self%alpha = 0._R_P
  self%is_filtered = .false.
  if (allocated(self%filter)) deallocate(self%filter)
  endsubroutine destroy

  subroutine initialize(self, scheme, nu, alpha, autoupdate, U, stop_on_fail)
  !< Create the actual leapfrog integrator: initialize the filter coefficient.
  class(integrator_leapfrog), intent(inout)        :: self         !< Integrator.
  character(*),               intent(in)           :: scheme       !< Selected scheme.
  real(R_P),                  intent(in), optional :: nu           !< Williams-Robert-Asselin filter coefficient.
  real(R_P),                  intent(in), optional :: alpha        !< Robert-Asselin filter coefficient.
  logical,                    intent(in), optional :: autoupdate   !< Enable cyclic autoupdate of previous time steps.
  class(integrand_object),    intent(in), optional :: U            !< Integrand molding prototype.
  logical,                    intent(in), optional :: stop_on_fail !< Stop execution if initialization fail.

  if (self%is_supported(scheme=scheme)) then
    call self%destroy
    select case(trim(adjustl(scheme)))
    case('leapfrog_raw')
       self%nu = 0.01_R_P ; if (present(nu)) self%nu = nu
       self%alpha = 0.53_R_P ; if (present(alpha)) self%alpha = alpha
       self%is_filtered = .true.
    endselect
    self%autoupdate = .true. ; if (present(autoupdate)) self%autoupdate = autoupdate
    self%steps = 2
    self%registers = self%steps
    if (present(U)) call self%allocate_integrand_members(U=U)
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=stop_on_fail)
  endif
  endsubroutine initialize

   ! overridden public methods
   pure subroutine allocate_integrand_members(self, U)
   !< Allocate members of interpolator being of [[integrand_object]] class.
   !<
   !< @note It is assumed that the integrator has been properly initialized before calling this method.
   class(integrator_leapfrog), intent(inout) :: self !< Integrator.
   class(integrand_object),    intent(in)    :: U    !< Integrand.
   integer(I_P)                              :: s    !< Counter.

   if (self%is_multistep() .and. self%registers > 0) then
      if (allocated(self%Dt)) deallocate(self%Dt)
      allocate(self%Dt(1:self%registers)) ; self%Dt = 0._R_P
      if (allocated(self%t)) deallocate(self%t)
      allocate(self%t(1:self%registers)) ; self%t = 0._R_P
      if (allocated(self%previous)) deallocate(self%previous)
      allocate(self%previous(1:self%registers), mold=U)
      do s=1, self%registers
         self%previous(s) = U
      enddo
   endif
   if (self%has_fast_mode()) then
      if (allocated(self%buffer)) deallocate(self%buffer)
      allocate(self%buffer, mold=U)
      self%buffer = U
   endif
   if (self%is_filtered) then
      if (allocated(self%filter)) deallocate(self%filter)
      allocate(self%filter, mold=U)
      self%filter = U
   endif
   endsubroutine allocate_integrand_members
endmodule foodie_integrator_leapfrog
