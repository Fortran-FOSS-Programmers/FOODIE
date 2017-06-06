!< FOODIE integrator: provide explicit Euler scheme, it being 1st order accurate.

module foodie_integrator_euler_explicit
!< FOODIE integrator: provide explicit Euler scheme, it being 1st order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the Forward Explicit Euler scheme implemented is:
!<
!< $$ U^{n+1} = U^n +\Delta t R(t, U^n) $$
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.

use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrand_object, only : integrand_object
use foodie_integrator_multistage_object, only : integrator_multistage_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_euler_explicit

character(len=99), parameter :: class_name_='euler_explicit'                !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:1)=[trim(class_name_)] !< List of supported schemes.

logical, parameter :: has_fast_mode_=.true. !< Flag to check if integrator provides *fast mode* integrate.

type, extends(integrator_multistage_object) :: integrator_euler_explicit
  !< FOODIE integrator: provide explicit Euler scheme, it being 1st order accurate.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(self) :: has_fast_mode        !< Return .true. if the integrator class has *fast mode* integrate.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: integrate            !< Integrate integrand field.
    procedure, pass(self) :: integrate_fast       !< Integrate integrand field, fast mode.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy    !< Destroy the integrator.
    procedure, pass(self) :: initialize !< Initialize (create) the integrator.
endtype integrator_euler_explicit
contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_euler_explicit), intent(in) :: self       !< Integrator.
  character(len=99)                            :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  elemental function has_fast_mode(self)
  !< Return .true. if the integrator class has *fast mode* integrate.
  class(integrator_euler_explicit), intent(in) :: self          !< Integrator.
  logical                                      :: has_fast_mode !< Inquire result.

  has_fast_mode = has_fast_mode_
  endfunction has_fast_mode

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_euler_explicit), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),         intent(in)    :: rhs !< Right hand side.

  call lhs%assign_multistage(rhs=rhs)
  endsubroutine integr_assign_integr

  subroutine integrate(self, U, Dt, t, new_Dt)
  !< Integrate field with explicit Euler scheme, 1st order.
  class(integrator_euler_explicit), intent(inout) :: self   !< Integrator.
  class(integrand_object),          intent(inout) :: U      !< Field to be integrated.
  real(R_P),                        intent(in)    :: Dt     !< Time step.
  real(R_P),                        intent(in)    :: t      !< Time.
  real(R_P), optional,              intent(out)   :: new_Dt !< New adapted time step.

  U = U + (U%t(t=t) * Dt)
  if (present(new_Dt)) new_Dt = Dt
  endsubroutine integrate

  subroutine integrate_fast(self, U, Dt, t, new_Dt)
  !< Integrate field with explicit Euler scheme, 1st order, fast mode.
  class(integrator_euler_explicit), intent(inout) :: self   !< Integrator.
  class(integrand_object),          intent(inout) :: U      !< Field to be integrated.
  real(R_P),                        intent(in)    :: Dt     !< Time step.
  real(R_P),                        intent(in)    :: t      !< Time.
  real(R_P), optional,              intent(out)   :: new_Dt !< New adapted time step.

  self%buffer = U
  call self%buffer%t_fast(t=t)
  call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt)
  call U%add_fast(lhs=U, rhs=self%buffer)
  if (present(new_Dt)) new_Dt = Dt
  endsubroutine integrate_fast

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_euler_explicit), intent(in) :: self         !< Integrator.
  character(*),                     intent(in) :: scheme       !< Selected scheme.
  logical                                      :: is_supported !< Inquire result.
  integer(I_P)                                 :: s            !< Counter.

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
  class(integrator_euler_explicit), intent(in) :: self       !< Integrator.
  character(len=99), allocatable               :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_euler_explicit), intent(inout) :: self !< Integrator.

  call self%destroy_multistage
  endsubroutine destroy

   subroutine initialize(self, scheme, U, stop_on_fail)
   !< Create the actual RK integrator: initialize the Butcher' table coefficients.
   class(integrator_euler_explicit), intent(inout)        :: self         !< Integrator.
   character(*),                     intent(in)           :: scheme       !< Selected scheme.
   class(integrand_object),          intent(in), optional :: U            !< Integrand molding prototype.
   logical,                          intent(in), optional :: stop_on_fail !< Stop execution if initialization fail.

   if (self%is_supported(scheme=scheme)) then
      call self%destroy
      self%description_ = trim(adjustl(scheme))
      self%stages = 0
      self%registers = self%stages
      if (present(U)) call self%allocate_integrand_members(U=U)
   else
      call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                              error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                              is_severe=stop_on_fail)
   endif
   endsubroutine initialize
endmodule foodie_integrator_euler_explicit
