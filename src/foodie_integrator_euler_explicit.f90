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

use foodie_integrand_object, only : integrand_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_euler_explicit

character(len=99), parameter :: class_name_='euler_explicit'                !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:1)=[trim(class_name_)] !< List of supported schemes.

type, extends(integrator_object) :: integrator_euler_explicit
  !< FOODIE integrator: provide explicit Euler scheme, it being 1st order accurate.
  !<
  !< @note The integrator can be used directly without any initialization.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(self) :: description          !< Return pretty-printed object description.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy   !< Destroy the integrator.
    procedure, nopass     :: integrate !< Integrate integrand field.
endtype integrator_euler_explicit
contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_euler_explicit), intent(in) :: self       !< Integrator.
  character(len=99)                            :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  class(integrator_euler_explicit), intent(in)           :: self    !< Integrator.
  character(*),                     intent(in), optional :: prefix  !< Prefixing string.
  character(len=:), allocatable                          :: desc    !< Description.
  character(len=:), allocatable                          :: prefix_ !< Prefixing string, local variable.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = prefix_//'Euler, Explicit (1 step/stage) 1st order scheme'
  endfunction description

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_euler_explicit), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),         intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  endsubroutine integr_assign_integr

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

  call self%destroy_abstract
  endsubroutine destroy

  subroutine integrate(U, Dt, t)
  !< Integrate field with explicit Euler scheme, 1st order.
  class(integrand_object), intent(inout) :: U  !< Field to be integrated.
  real(R_P),               intent(in)    :: Dt !< Time step.
  real(R_P), optional,     intent(in)    :: t  !< Time.

  U = U + (U%t(t=t) * Dt)
  endsubroutine integrate
endmodule foodie_integrator_euler_explicit
