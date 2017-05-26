!< Define the abstract type [[integrator_object]] of FOODIE ODE integrators.

module foodie_integrator_object
!< Define the abstract type [[integrator_object]] of FOODIE ODE integrators.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P

implicit none
private
public :: integrator_object

type, abstract :: integrator_object
  !< Abstract type of FOODIE ODE integrators.
  character(len=:), allocatable :: description_  !< Informative description of the integrator.
  integer(I_P)                  :: error=0       !< Error status code.
  character(len=:), allocatable :: error_message !< Error message, hopefully meaningful.
  contains
    ! public methods
    procedure, pass(lhs)  :: assign_abstract  !< Assign ony members of abstract [[integrator_object]] type.
    procedure, pass(self) :: check_error      !< Check for error occurrencies.
    procedure, pass(self) :: description      !< Return informative integrator description.
    procedure, pass(self) :: destroy_abstract !< Destroy only members of abstract [[integrator_object]] type.
    procedure, pass(self) :: trigger_error    !< Trigger an error.
    ! deferred methods
    procedure(class_name_interface),        pass(self), deferred :: class_name           !< Return the class name of schemes.
    procedure(has_fast_mode_interface),     pass(self), deferred :: has_fast_mode        !< Return .true. if the integrator class
                                                                                         !< has *fast mode* integrate.
    procedure(assignment_interface),        pass(lhs),  deferred :: integr_assign_integr !< Operator `=`.
    procedure(is_multistagestep_interface), pass(self), deferred :: is_multistage        !< Return .true. for multistage integrator.
    procedure(is_multistagestep_interface), pass(self), deferred :: is_multistep         !< Return .true. for multistep integrator.
    procedure(is_supported_interface),      pass(self), deferred :: is_supported         !< Return .true. if the integrator class
                                                                                         !< support the given scheme.
    procedure(stagesteps_number_interface), pass(self), deferred :: stages_number        !< Return number of stages used.
    procedure(stagesteps_number_interface), pass(self), deferred :: steps_number         !< Return number of steps used.
    procedure(supported_schemes_interface), pass(self), deferred :: supported_schemes    !< Return the list of supported schemes.
    ! operators
    generic :: assignment(=) => integr_assign_integr !< Overload `=`.
endtype integrator_object

abstract interface
  !< Abstract interfaces of deferred methods of [[integrator_object]].
  pure function class_name_interface(self) result(class_name)
  !< Return the class name of schemes.
  import :: integrator_object
  class(integrator_object), intent(in) :: self       !< Integrator.
  character(len=99)                    :: class_name !< Class name.
  endfunction class_name_interface

  pure function description_interface(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  import :: integrator_object
  class(integrator_object), intent(in)           :: self   !< Integrator.
  character(*),             intent(in), optional :: prefix !< Prefixing string.
  character(len=:), allocatable                  :: desc   !< Description.
  endfunction description_interface

  pure subroutine assignment_interface(lhs, rhs)
  !< Operator `=`.
  import :: integrator_object
  class(integrator_object), intent(inout) :: lhs !< Left hand side.
  class(integrator_object), intent(in)    :: rhs !< Right hand side.
  endsubroutine assignment_interface

  elemental function has_fast_mode_interface(self) result(has_fast_mode)
  !< Return .true. if the integrator class has *fast mode* integrate.
  import :: integrator_object
  class(integrator_object), intent(in) :: self          !< Integrator.
  logical                              :: has_fast_mode !< Inquire result.
  endfunction has_fast_mode_interface

  elemental function is_multistagestep_interface(self) result(is_multistagestep)
  !< Return .true. for multistage or multistep integrator.
  import :: integrator_object
  class(integrator_object), intent(in) :: self              !< Integrator.
  logical                              :: is_multistagestep !< Inquire result.
  endfunction is_multistagestep_interface

  elemental function is_supported_interface(self, scheme) result(is_supported)
  !< Return .true. if the integrator class support the given scheme.
  import :: integrator_object
  class(integrator_object), intent(in) :: self         !< Integrator.
  character(*),             intent(in) :: scheme       !< Queried scheme.
  logical                              :: is_supported !< Inquire result.
  endfunction is_supported_interface

  elemental function stagesteps_number_interface(self) result(stagesteps_number)
  !< Return number of stages/steps used.
  import :: integrator_object, I_P
  class(integrator_object), intent(in) :: self              !< Integrator.
  integer(I_P)                         :: stagesteps_number !< Inquire result.
  endfunction stagesteps_number_interface

  pure function supported_schemes_interface(self) result(schemes)
  !< Return the list of supported schemes.
  import :: integrator_object
  class(integrator_object), intent(in) :: self       !< Integrator.
  character(len=99), allocatable       :: schemes(:) !< Queried scheme.
  endfunction supported_schemes_interface
endinterface

contains
  ! public methods
  pure subroutine assign_abstract(lhs, rhs)
  !< Assign ony members of abstract [[integrator_object]] type.
  class(integrator_object), intent(inout) :: lhs !< Left hand side.
  class(integrator_object), intent(in)    :: rhs !< Right hand side.

  if (allocated(rhs%description_ )) lhs%description_  = rhs%description_
                                    lhs%error         = rhs%error
  if (allocated(rhs%error_message)) lhs%error_message = rhs%error_message
  endsubroutine assign_abstract

  subroutine check_error(self, is_severe)
  !< Check for error occurencies.
  !<
  !< If `is_severe=.true.` an stop is called.
  class(integrator_object), intent(in)           :: self      !< Integrator.
  logical,                  intent(in), optional :: is_severe !< Flag to activate severe faliure, namely errors trigger a stop.

  if (self%error /= 0) then
    if (allocated(self%error_message)) then
      write(stderr, '(A)')'error: '//self%error_message
    else
      write(stderr, '(A)')'error: an obscure error occurred!'
    endif
    write(stderr, '(A,I4)')'error code: ', self%error
    if (present(is_severe)) then
      if (is_severe) stop
    endif
  endif
  endsubroutine check_error

   pure function description(self, prefix) result(desc)
   !< Return informative integrator description.
   class(integrator_object), intent(in)           :: self    !< Integrator.
   character(*),             intent(in), optional :: prefix  !< Prefixing string.
   character(len=:), allocatable                  :: desc    !< Description.
   character(len=:), allocatable                  :: prefix_ !< Prefixing string, local variable.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   if (allocated(self%description_)) then
      desc = prefix//self%description_
   else
      desc = prefix//self%class_name()
   endif
   endfunction description

  elemental subroutine destroy_abstract(self)
  !< Destroy only members of abstract [[integrator_object]] type.
  class(integrator_object), intent(inout) :: self !< Integrator.

  if (allocated(self%description_)) deallocate(self%description_)
  self%error = 0
  if (allocated(self%error_message)) deallocate(self%error_message)
  endsubroutine destroy_abstract

  subroutine trigger_error(self, error, error_message, is_severe)
  !< Check for error occurencies.
  !<
  !< If `is_severe=.true.` an stop is called.
  class(integrator_object), intent(inout)        :: self          !< Integrator.
  integer(I_P),             intent(in)           :: error         !< Error status code.
  character(len=*),         intent(in), optional :: error_message !< Error message, hopefully meaningful.
  logical,                  intent(in), optional :: is_severe     !< Flag to activate severe faliure, namely errors trigger a stop.

                              self%error         = error
  if (present(error_message)) self%error_message = error_message
  call self%check_error(is_severe=is_severe)
  endsubroutine trigger_error
endmodule foodie_integrator_object
