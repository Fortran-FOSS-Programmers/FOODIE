#include "preprocessor_macros.h"
!< Define the abstract type [[integrator_multistage_multistep_object]] of FOODIE ODE integrators.

module foodie_integrator_multistage_multistep_object
!< Define the abstract type [[integrator_multistage_multistep_object]] of FOODIE ODE integrators.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foodie_integrand_object, only : integrand_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_multistage_multistep_object

type, extends(integrator_object), abstract :: integrator_multistage_multistep_object
   !< Abstract type of FOODIE ODE integrators of the multistage/multistep family.
   integer(I_P)                         :: registers_stages !< Number of registers used for stages.
   integer(I_P)                         :: registers_steps  !< Number of registers used for steps.
   integer(I_P)                         :: stages           !< Number of stages.
   integer(I_P)                         :: steps            !< Number of time steps.
   logical                              :: autoupdate       !< Perform cyclic autoupdate of previous time steps buffers.
   integer(I_P)                         :: iterations       !< Implicit iterations.
   real(R_P),               allocatable :: Dt(:)            !< Previous time steps.
   real(R_P),               allocatable :: t(:)             !< Previous times.
   class(integrand_object), allocatable :: previous(:)      !< Previous steps.
   class(integrand_object), allocatable :: stage(:)         !< Stages.
   class(integrand_object), allocatable :: buffer           !< Buffer used for fast integration.
   contains
      ! deferred methods
      procedure(integrate_interface),      pass(self), deferred :: integrate      !< Integrate integrand field.
      procedure(integrate_fast_interface), pass(self), deferred :: integrate_fast !< Integrate integrand field, fast mode.
      ! implemented deferred methods of parent
      procedure, pass(self) :: is_multistage !< Return .true. for multistage integrator.
      procedure, pass(self) :: is_multistep  !< Return .true. for multistep integrator.
      procedure, pass(self) :: stages_number !< Return number of stages used.
      procedure, pass(self) :: steps_number  !< Return number of steps used.
      ! public methods
      procedure, pass(self) :: allocate_integrand_members   !< Allocate integrand members.
      procedure, pass(lhs)  :: assign_multistage_multistep  !< Assign members of [[integrator_multistage_multistep_object]].
      procedure, pass(self) :: destroy_multistage_multistep !< Destroy the integrator.
      procedure, nopass     :: update_previous              !< Cyclic update previous time steps.
endtype integrator_multistage_multistep_object

abstract interface
   !< Abstract interfaces of deferred methods of [[integrator_multistage_multistep_object]].
   subroutine integrate_interface(self, U, Dt, t)
   !< Integrate integrand field.
   import :: integrand_object, integrator_multistage_multistep_object, R_P
   class(integrator_multistage_multistep_object), intent(inout) :: self !< Integrator.
   class(integrand_object),                       intent(inout) :: U    !< Integrand.
   real(R_P),                                     intent(in)    :: Dt   !< Time step.
   real(R_P),                                     intent(in)    :: t    !< Time.
   endsubroutine integrate_interface

   subroutine integrate_fast_interface(self, U, Dt, t)
   !< Integrate integrand field, fast mode.
   import :: integrand_object, integrator_multistage_multistep_object, R_P
   class(integrator_multistage_multistep_object), intent(inout) :: self !< Integrator.
   class(integrand_object),                       intent(inout) :: U    !< Field to be integrated.
   real(R_P),                                     intent(in)    :: Dt   !< Time step.
   real(R_P),                                     intent(in)    :: t    !< Time.
   endsubroutine integrate_fast_interface
endinterface

contains
   ! deferred methods
   elemental function is_multistage(self)
   !< Return .true. for multistage integrator.
   class(integrator_multistage_multistep_object), intent(in) :: self          !< Integrator.
   logical                                                   :: is_multistage !< Inquire result.

   is_multistage = .true.
   endfunction is_multistage

   elemental function is_multistep(self)
   !< Return .true. for multistage integrator.
   class(integrator_multistage_multistep_object), intent(in) :: self         !< Integrator.
   logical                                                   :: is_multistep !< Inquire result.

   is_multistep = .true.
   endfunction is_multistep

   elemental function stages_number(self)
   !< Return number of stages used.
   class(integrator_multistage_multistep_object), intent(in) :: self          !< Integrator.
   integer(I_P)                                              :: stages_number !< Number of stages used.

   stages_number = self%stages
   endfunction stages_number

   elemental function steps_number(self)
   !< Return number of steps used.
   class(integrator_multistage_multistep_object), intent(in) :: self         !< Integrator.
   integer(I_P)                                              :: steps_number !< Number of steps used.

   steps_number = self%steps
   endfunction steps_number

   ! public methods
   subroutine allocate_integrand_members(self, U)
   !< Allocate members of interpolator being of [[integrand_object]] class.
   !<
   !< @note It is assumed that the integrator has been properly initialized before calling this method.
   class(integrator_multistage_multistep_object), intent(inout) :: self !< Integrator.
   class(integrand_object),                       intent(in)    :: U    !< Integrand.
   integer(I_P)                                                 :: s    !< Counter.

   if (self%is_multistage() .and. self%registers_stages > 0) then
      if (allocated(self%stage)) deallocate(self%stage)
      allocate(self%stage(1:self%registers_stages), mold=U)
      do s=1, self%registers_stages
         self%stage(s) = U
      enddo
   endif
   if (self%is_multistep() .and. self%registers_steps > 0) then
      if (allocated(self%Dt)) deallocate(self%Dt)
      allocate(self%Dt(1:self%registers_steps)) ; self%Dt = 0._R_P
      if (allocated(self%t)) deallocate(self%t)
      allocate(self%t(1:self%registers_steps)) ; self%t = 0._R_P
      if (allocated(self%previous)) deallocate(self%previous)
      allocate(self%previous(1:self%registers_steps), mold=U)
      do s=1, self%registers_steps
         self%previous(s) = U
      enddo
   endif
   if (self%has_fast_mode()) then
      if (allocated(self%buffer)) deallocate(self%buffer)
      allocate(self%buffer, mold=U)
      self%buffer = U
   endif
   endsubroutine allocate_integrand_members

   _PURE_ subroutine assign_multistage_multistep(lhs, rhs)
   !< Assign members of [[integrator_multistage_multistep_object]] and parents.
   class(integrator_multistage_multistep_object), intent(inout) :: lhs !< Left hand side.
   class(integrator_object),                      intent(in)    :: rhs !< Right hand side.
   integer(I_P)                                                 :: s   !< Counter.

   call lhs%assign_abstract(rhs=rhs)
   select type(rhs)
   class is (integrator_multistage_multistep_object)
     lhs%registers_stages = rhs%registers_stages
     lhs%registers_steps = rhs%registers_steps
     lhs%stages = rhs%stages
     lhs%steps = rhs%steps
     lhs%autoupdate = rhs%autoupdate
     lhs%iterations = rhs%iterations
     if (allocated(lhs%Dt)) deallocate(lhs%Dt)
     if (allocated(rhs%Dt)) lhs%Dt = rhs%Dt
     if (allocated(lhs%t)) deallocate(lhs%t)
     if (allocated(rhs%t)) lhs%t = rhs%t
     if (allocated(lhs%previous)) deallocate(lhs%previous)
     if (allocated(rhs%previous)) then
        allocate(lhs%previous(1:lhs%registers_steps), mold=rhs%previous)
        do s=1, lhs%registers_steps
           lhs%previous(s) = rhs%previous(s)
        enddo
     endif
     if (allocated(lhs%stage)) deallocate(lhs%stage)
     if (allocated(rhs%stage)) then
        allocate(lhs%stage(1:lhs%registers_stages), mold=rhs%stage)
        do s=1, lhs%registers_stages
           lhs%stage(s) = rhs%stage(s)
        enddo
     endif
     if (allocated(lhs%buffer)) deallocate(lhs%buffer)
     if (allocated(rhs%buffer)) then
        allocate(lhs%buffer, mold=rhs%buffer)
        lhs%buffer = rhs%buffer
     endif
   endselect
   endsubroutine assign_multistage_multistep

   elemental subroutine destroy_multistage_multistep(self)
   !< Destroy the integrator.
   class(integrator_multistage_multistep_object), intent(inout) :: self !< Integrator.

   call self%destroy_abstract
   self%registers_stages = 0
   self%registers_steps = 0
   self%stages = 0
   self%steps = -1
   self%autoupdate = .false.
   self%iterations = 0
   if (allocated(self%Dt)) deallocate(self%Dt)
   if (allocated(self%t)) deallocate(self%t)
   if (allocated(self%previous)) deallocate(self%previous)
   if (allocated(self%stage)) deallocate(self%stage)
   if (allocated(self%buffer)) deallocate(self%buffer)
   endsubroutine destroy_multistage_multistep

   subroutine update_previous(U, previous, Dt, t, previous_Dt, previous_t)
   !< Cyclic update previous time steps.
   class(integrand_object), intent(in)              :: U               !< Field to be integrated.
   class(integrand_object), intent(inout)           :: previous(1:)    !< Previous time steps solutions of integrand.
   real(R_P),               intent(in),    optional :: Dt              !< Time step.
   real(R_P),               intent(in),    optional :: t               !< Time.
   real(R_P),               intent(inout), optional :: previous_Dt(1:) !< Time step.
   real(R_P),               intent(inout), optional :: previous_t(1:)  !< Time.
   integer(I_P)                                     :: last_step       !< Last step.
   integer(I_P)                                     :: s               !< Steps counter.

   last_step = size(previous, dim=1)
   do s=1, last_step - 1
     previous(s) = previous(s + 1)
     if (present(previous_Dt)) previous_Dt(s) = previous_Dt(s + 1)
     if (present(previous_t)) previous_t(s) = previous_t(s + 1)
   enddo
   previous(last_step) = U
   if (present(previous_Dt)) previous_Dt(last_step) = Dt
   if (present(previous_t)) previous_t(last_step) = t + Dt
   endsubroutine update_previous
endmodule foodie_integrator_multistage_multistep_object
