!< Define the abstract type *integrator* of FOODIE ODE integrators.

module foodie_integrator_object
!< Define the abstract type *integrator* of FOODIE ODE integrators.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foodie_kinds, only : I_P

implicit none
private
public :: integrator_object

type, abstract :: integrator_object
  !< Abstract type of FOODIE ODE integrators.
  integer(I_P)                  :: error=0       !< Error status code.
  character(len=:), allocatable :: error_message !< Error message, hopefully meaningful.
  contains
    ! public methods
    procedure, pass(lhs)  :: assign_abstract  !< Assign ony members of abstract [[integrator_object]] type.
    procedure, pass(self) :: check_error      !< Check for error occurrencies.
    procedure, pass(self) :: destroy_abstract !< Destroy only members of abstract [[integrator_object]] type.
    procedure, pass(self) :: trigger_error    !< Trigger an error.
    ! deferred methods
    procedure(description_interface), pass(self), deferred :: description          !< Return pretty-printed object description.
    procedure(assignment_interface),  pass(lhs),  deferred :: integr_assign_integr !< Operator `=`.
    ! operators
    generic :: assignment(=) => integr_assign_integr !< Overload `=`.
endtype integrator_object

abstract interface
  !< Abstract interfaces of deferred methods of [[interpolator_object]].
  pure function description_interface(self, prefix) result(desc)
  !< Return a pretty-formatted object description.
  import :: integrator_object
  class(integrator_object), intent(in)           :: self   !< Interpolator.
  character(*),             intent(in), optional :: prefix !< Prefixing string.
  character(len=:), allocatable                  :: desc   !< Description.
  endfunction description_interface

  pure subroutine assignment_interface(lhs, rhs)
  !< Operator `=`.
  import :: integrator_object
  class(integrator_object), intent(inout) :: lhs !< Left hand side.
  class(integrator_object), intent(in)    :: rhs !< Right hand side.
  endsubroutine assignment_interface
endinterface

contains
  ! public methods
  pure subroutine assign_abstract(lhs, rhs)
  !< Assign ony members of abstract [[integrator_object]] type.
  class(integrator_object), intent(inout) :: lhs !< Left hand side.
  class(integrator_object), intent(in)    :: rhs !< Right hand side.

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

  elemental subroutine destroy_abstract(self)
  !< Destroy only members of abstract [[integrator_object]] type.
  class(integrator_object), intent(inout) :: self !< Integrator.

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
