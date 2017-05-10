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

logical, parameter :: has_fast_mode_=.true.  !< Flag to check if integrator provides *fast mode* integrate.
logical, parameter :: is_multistage_=.false. !< Flag to check if integrator is multistage.
logical, parameter :: is_multistep_=.false.  !< Flag to check if integrator is multistep.

type, extends(integrator_object) :: integrator_euler_explicit
  !< FOODIE integrator: provide explicit Euler scheme, it being 1st order accurate.
  !<
  !< @note The integrator can be used directly without any initialization.
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
    procedure, pass(self) :: destroy        !< Destroy the integrator.
    procedure, nopass     :: integrate      !< Integrate integrand field.
    procedure, nopass     :: integrate_fast !< Integrate integrand field, fast mode.
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

  call lhs%assign_abstract(rhs=rhs)
  endsubroutine integr_assign_integr

  elemental function is_multistage(self)
  !< Return .true. for multistage integrator.
  class(integrator_euler_explicit), intent(in) :: self          !< Integrator.
  logical                                      :: is_multistage !< Inquire result.

  is_multistage = is_multistage_
  endfunction is_multistage

  elemental function is_multistep(self)
  !< Return .true. for multistage integrator.
  class(integrator_euler_explicit), intent(in) :: self         !< Integrator.
  logical                                      :: is_multistep !< Inquire result.

  is_multistep = is_multistep_
  endfunction is_multistep

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

  elemental function stages_number(self)
  !< Return number of stages used.
  class(integrator_euler_explicit), intent(in) :: self          !< Integrator.
  integer(I_P)                                 :: stages_number !< Number of stages used.

  stages_number = 1
  endfunction stages_number

  elemental function steps_number(self)
  !< Return number of steps used.
  class(integrator_euler_explicit), intent(in) :: self         !< Integrator.
  integer(I_P)                                 :: steps_number !< Number of steps used.

  steps_number = 1
  endfunction steps_number

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

  subroutine integrate_fast(U, buffer, Dt, t)
  !< Integrate field with explicit Euler scheme, 1st order, fast mode.
  class(integrand_object), intent(inout) :: U      !< Field to be integrated.
  class(integrand_object), intent(inout) :: buffer !< Temporary buffer for doing fast operation.
  real(R_P),               intent(in)    :: Dt     !< Time step.
  real(R_P), optional,     intent(in)    :: t      !< Time.

  buffer = U
  call buffer%t_fast(t=t)
  call buffer%multiply_fast(lhs=buffer, rhs=Dt)
  call U%add_fast(lhs=U, rhs=buffer)
  endsubroutine integrate_fast
endmodule foodie_integrator_euler_explicit
