!< Define the abstract type [[integrator_multistep_implicit_object]] of FOODIE ODE integrators.

module foodie_integrator_multistep_implicit_object
!< Define the abstract type [[integrator_multistep_implicit_object]] of FOODIE ODE integrators.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foodie_integrand_object, only : integrand_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_multistep_implicit_object

type, extends(integrator_object), abstract :: integrator_multistep_implicit_object
   !< Abstract type of FOODIE ODE integrators of the multistep-implicit family.
   integer(I_P)                         :: registers   !< Number of registers used for stages.
   integer(I_P)                         :: steps       !< Number of time steps.
   class(integrand_object), allocatable :: previous(:) !< Previous steps.
   class(integrand_object), allocatable :: buffer      !< Buffer used for fast integration.
   contains
      ! deferred methods
      procedure(integrate_interface),         pass(self), deferred :: integrate         !< Integrate integrand field.
      procedure(integrate_fast_interface),    pass(self), deferred :: integrate_fast    !< Integrate integrand field, fast mode.
      procedure(integrate_ub_interface),      pass(self), deferred :: integrate_ub      !< Integrate integrand field, unbuffered.
      procedure(integrate_ub_fast_interface), pass(self), deferred :: integrate_ub_fast !< Integrate integrand field, fast mode,
                                                                                        !< unbuffered.
      ! implemented deferred methods of parent
      procedure, pass(self) :: is_multistage !< Return .true. for multistage integrator.
      procedure, pass(self) :: is_multistep  !< Return .true. for multistep integrator.
      procedure, pass(self) :: stages_number !< Return number of stages used.
      procedure, pass(self) :: steps_number  !< Return number of steps used.
      ! public methods
      procedure, pass(self) :: allocate_integrand_members !< Allocate integrand members.
      procedure, pass(self) :: destroy_multistep          !< Destroy the integrator.
      procedure, pass(self) :: update_previous            !< Cyclic update previous time steps.
endtype integrator_multistep_implicit_object

abstract interface
   !< Abstract interfaces of deferred methods of [[integrator_multistep_implicit_object]].
   subroutine integrate_interface(self, U, Dt, t, iterations, autoupdate)
   !< Integrate integrand field.
   !<
   !< @note This method uses integrand previous-steps-buffer stored inside integrator.
   import :: integrand_object, integrator_multistep_implicit_object, I_P, R_P
   class(integrator_multistep_implicit_object), intent(inout) :: self       !< Integrator.
   class(integrand_object),                     intent(inout) :: U          !< Integrand.
   real(R_P),                                   intent(in)    :: Dt(1:)     !< Time steps.
   real(R_P),                                   intent(in)    :: t(1:)      !< Times.
   integer(I_P), optional,                      intent(in)    :: iterations !< Fixed point iterations.
   logical,      optional,                      intent(in)    :: autoupdate !< Perform cyclic autoupdate of previous time steps.
   endsubroutine integrate_interface

   subroutine integrate_fast_interface(self, U, Dt, t, iterations, autoupdate)
   !< Integrate integrand field, fast mode.
   !<
   !< @note This method uses integrand previous-steps-buffer stored inside integrator.
   import :: integrand_object, integrator_multistep_implicit_object, I_P, R_P
   class(integrator_multistep_implicit_object), intent(inout) :: self       !< Integrator.
   class(integrand_object),                     intent(inout) :: U          !< Field to be integrated.
   real(R_P),                                   intent(in)    :: Dt(1:)     !< Time steps.
   real(R_P),                                   intent(in)    :: t(1:)      !< Times.
   integer(I_P), optional,                      intent(in)    :: iterations !< Fixed point iterations.
   logical, optional,                           intent(in)    :: autoupdate !< Perform cyclic autoupdate of previous time steps.
   endsubroutine integrate_fast_interface

   subroutine integrate_ub_interface(self, U, previous, Dt, t, iterations, autoupdate)
   !< Integrate integrand field, unbuffered.
   import :: integrand_object, integrator_multistep_implicit_object, I_P, R_P
   class(integrator_multistep_implicit_object), intent(inout) :: self         !< Integrator.
   class(integrand_object),                     intent(inout) :: U            !< Integrand.
   class(integrand_object),                     intent(inout) :: previous(1:) !< Integrand.
   real(R_P),                                   intent(in)    :: Dt(1:)       !< Time steps.
   real(R_P),                                   intent(in)    :: t(1:)        !< Times.
   integer(I_P), optional,                      intent(in)    :: iterations   !< Fixed point iterations.
   logical, optional,                           intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
   endsubroutine integrate_ub_interface

   subroutine integrate_ub_fast_interface(self, U, previous, Dt, t, iterations, autoupdate)
   !< Integrate integrand field, unbuffered, fast mode.
   import :: integrand_object, integrator_multistep_implicit_object, I_P, R_P
   class(integrator_multistep_implicit_object), intent(inout) :: self         !< Integrator.
   class(integrand_object),                     intent(inout) :: U            !< Field to be integrated.
   class(integrand_object),                     intent(inout) :: previous(1:) !< Integrand.
   real(R_P),                                   intent(in)    :: Dt(1:)       !< Time steps.
   real(R_P),                                   intent(in)    :: t(1:)        !< Times.
   integer(I_P), optional,                      intent(in)    :: iterations   !< Fixed point iterations.
   logical, optional,                           intent(in)    :: autoupdate   !< Perform cyclic autoupdate of previous time steps.
   endsubroutine integrate_ub_fast_interface
endinterface

contains
   ! deferred methods
   elemental function is_multistage(self)
   !< Return .true. for multistage integrator.
   class(integrator_multistep_implicit_object), intent(in) :: self          !< Integrator.
   logical                                                 :: is_multistage !< Inquire result.

   is_multistage = .false.
   endfunction is_multistage

   elemental function is_multistep(self)
   !< Return .true. for multistage integrator.
   class(integrator_multistep_implicit_object), intent(in) :: self         !< Integrator.
   logical                                                 :: is_multistep !< Inquire result.

   is_multistep = .true.
   endfunction is_multistep

   elemental function stages_number(self)
   !< Return number of stages used.
   class(integrator_multistep_implicit_object), intent(in) :: self          !< Integrator.
   integer(I_P)                                            :: stages_number !< Number of stages used.

   stages_number = 0
   endfunction stages_number

   elemental function steps_number(self)
   !< Return number of steps used.
   class(integrator_multistep_implicit_object), intent(in) :: self         !< Integrator.
   integer(I_P)                                            :: steps_number !< Number of steps used.

   steps_number = self%steps
   endfunction steps_number

   ! public methods
   pure subroutine allocate_integrand_members(self, U)
   !< Allocate members of interpolator being of [[integrand_object]] class.
   !<
   !< @note It is assumed that the integrator has been properly initialized before calling this method.
   class(integrator_multistep_implicit_object), intent(inout) :: self !< Integrator.
   class(integrand_object),                     intent(in)    :: U    !< Integrand.
   integer(I_P)                                               :: s    !< Counter.

   if (self%is_multistep() .and. self%registers > 0) then
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
   endsubroutine allocate_integrand_members

   elemental subroutine destroy_multistep(self)
   !< Destroy the integrator.
   class(integrator_multistep_implicit_object), intent(inout) :: self !< Integrator.

   call self%destroy_abstract
   self%registers = 0
   self%steps = -1
   if (allocated(self%previous)) deallocate(self%previous)
   if (allocated(self%buffer)) deallocate(self%buffer)
   endsubroutine destroy_multistep

   subroutine update_previous(self, U, previous, is_like_explicit)
   !< Cyclic update previous time steps.
   class(integrator_multistep_implicit_object), intent(in)    :: self              !< Integrator.
   class(integrand_object),                     intent(in)    :: U                 !< Field to be integrated.
   class(integrand_object),                     intent(inout) :: previous(1:)      !< Previous time steps solutions of integrand.
   logical, optional,                           intent(in)    :: is_like_explicit  !< Use explicit-like update.
   logical                                                    :: is_like_explicit_ !< Use explicit-like update, local variable.
   integer(I_P)                                               :: s                 !< Steps counter.

   is_like_explicit_ = .false. ; if (present(is_like_explicit)) is_like_explicit_ = is_like_explicit
   if (is_like_explicit_) then
      do s=1, self%steps - 1
        previous(s) = previous(s + 1)
      enddo
      previous(self%steps) = U
   else
      if (self%steps > 0) then
        do s=0, self%steps - 2
          previous(s + 1) = previous(s + 2)
        enddo
        previous(self%steps) = U
      endif
   endif
   endsubroutine update_previous
endmodule foodie_integrator_multistep_implicit_object
