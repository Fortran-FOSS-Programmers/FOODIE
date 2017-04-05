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

use foodie_adt_integrand, only : integrand
use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_kinds, only : I_P, R_P
use foodie_integrator_object, only : integrator_object

implicit none
private
public :: integrator_leapfrog

character(len=99), parameter :: class_name_='leapfrog'                              !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:2)=[trim(class_name_)//'    ', &
                                                         trim(class_name_)//'_raw'] !< List of supported schemes.

type, extends(integrator_object) :: integrator_leapfrog
  !< FOODIE integrator: provide an explicit class of leapfrog multi-step schemes, 2nd order accurate.
  !<
  !< @note The integrator could be used without initialialization (initialize the time filter coefficients) if the defulat values
  !< are suitable for the problem.
  private
  integer(I_P), public :: steps=2        !< Number of time steps.
  real(R_P)            :: nu=0.01_R_P    !< Robert-Asselin filter coefficient.
  real(R_P)            :: alpha=0.53_R_P !< Robert-Asselin-Williams filter coefficient.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy    !< Destroy the integrator.
    procedure, pass(self) :: initialize !< Initialize (create) the integrator.
    procedure, pass(self) :: integrate  !< Integrate integrand field.
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

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_leapfrog), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),   intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  select type(rhs)
  class is(integrator_leapfrog)
    lhs%steps = rhs%steps
    lhs%nu    = rhs%nu
    lhs%alpha = rhs%alpha
  endselect
  endsubroutine integr_assign_integr

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

  call self%destroy_abstract
  self%steps = 2
  self%nu = 0.01_R_P
  self%alpha = 0.53_R_P
  endsubroutine destroy

  subroutine initialize(self, scheme, nu, alpha)
  !< Create the actual leapfrog integrator: initialize the filter coefficient.
  class(integrator_leapfrog), intent(inout)        :: self   !< Integrator.
  character(*),               intent(in)           :: scheme !< Selected scheme.
  real(R_P),                  intent(in), optional :: nu     !< Williams-Robert-Asselin filter coefficient.
  real(R_P),                  intent(in), optional :: alpha  !< Robert-Asselin filter coefficient.

  if (self%is_supported(scheme=scheme)) then
    call self%destroy
    if (present(nu)) self%nu = nu
    if (present(alpha)) self%alpha = alpha
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=.true.)
  endif
  endsubroutine initialize

  subroutine integrate(self, U, previous, Dt, t, filter)
  !< Integrate field with leapfrog class scheme.
  class(integrator_leapfrog), intent(in)    :: self          !< Integrator.
  class(integrand),           intent(inout) :: U             !< Field to be integrated.
  class(integrand),           intent(inout) :: previous(1:2) !< Previous time steps solutions of integrand field.
  real(R_P),                  intent(in)    :: Dt            !< Time step.
  real(R_P),                  intent(in)    :: t             !< Time.
  class(integrand), optional, intent(inout) :: filter        !< Filter field displacement.

  U = previous(1) + previous(2)%t(t=t) * (Dt * 2._R_P)
  if (present(filter)) then
    filter = (previous(1) - previous(2) * 2._R_P + U) * self%nu * 0.5_R_P
    previous(2) = previous(2) + filter * self%alpha
    U = U + filter * (self%alpha - 1._R_P)
  endif
  previous(1) = previous(2)
  previous(2) = U
  endsubroutine integrate
endmodule foodie_integrator_leapfrog
