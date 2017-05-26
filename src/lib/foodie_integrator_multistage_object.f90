!< Define the abstract type [[integrator_multistage_object]] of FOODIE ODE integrators.

module foodie_integrator_multistage_object
!< Define the abstract type [[integrator_multistage_object]] of FOODIE ODE integrators.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foodie_integrand_object, only : integrand_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_multistage_object

type, extends(integrator_object), abstract :: integrator_multistage_object
   !< Abstract type of FOODIE ODE integrators of the multistage family.
   integer(I_P)                         :: registers !< Number of registers used for stages.
   integer(I_P)                         :: stages    !< Number of stages.
   class(integrand_object), allocatable :: stage(:)  !< Stages.
   class(integrand_object), allocatable :: buffer    !< Buffer used for fast integration.
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
      procedure, pass(self) :: allocate_integrand_members !< Allocate integrand members.
      procedure, pass(lhs)  :: assign_multistage          !< Assign members of [[integrator_multistage_object]] and parents.
      procedure, pass(self) :: destroy_multistage         !< Destroy the integrator.
endtype integrator_multistage_object

abstract interface
   !< Abstract interfaces of deferred methods of [[integrator_multistage_object]].
   subroutine integrate_interface(self, U, Dt, t, new_Dt)
   !< Integrate integrand field.
   import :: integrand_object, integrator_multistage_object, R_P
   class(integrator_multistage_object), intent(inout)         :: self   !< Integrator.
   class(integrand_object),             intent(inout)         :: U      !< Integrand.
   real(R_P),                           intent(in)            :: Dt     !< Time step.
   real(R_P),                           intent(in)            :: t      !< Time.
   real(R_P),                           intent(out), optional :: new_Dt !< New adapted time step.
   endsubroutine integrate_interface

   subroutine integrate_fast_interface(self, U, Dt, t, new_Dt)
   !< Integrate integrand field, fast mode.
   import :: integrand_object, integrator_multistage_object, R_P
   class(integrator_multistage_object), intent(inout)         :: self   !< Integrator.
   class(integrand_object),             intent(inout)         :: U      !< Field to be integrated.
   real(R_P),                           intent(in)            :: Dt     !< Time step.
   real(R_P),                           intent(in)            :: t      !< Time.
   real(R_P),                           intent(out), optional :: new_Dt !< New adapted time step.
   endsubroutine integrate_fast_interface
endinterface

contains
   ! deferred methods
   elemental function is_multistage(self)
   !< Return .true. for multistage integrator.
   class(integrator_multistage_object), intent(in) :: self          !< Integrator.
   logical                                         :: is_multistage !< Inquire result.

   is_multistage = .true.
   endfunction is_multistage

   elemental function is_multistep(self)
   !< Return .true. for multistage integrator.
   class(integrator_multistage_object), intent(in) :: self         !< Integrator.
   logical                                         :: is_multistep !< Inquire result.

   is_multistep = .false.
   endfunction is_multistep

   elemental function stages_number(self)
   !< Return number of stages used.
   class(integrator_multistage_object), intent(in) :: self          !< Integrator.
   integer(I_P)                                    :: stages_number !< Number of stages used.

   stages_number = self%stages
   endfunction stages_number

   elemental function steps_number(self)
   !< Return number of steps used.
   class(integrator_multistage_object), intent(in) :: self         !< Integrator.
   integer(I_P)                                    :: steps_number !< Number of steps used.

   steps_number = 0
   endfunction steps_number

   ! public methods
   pure subroutine allocate_integrand_members(self, U)
   !< Allocate members of interpolator being of [[integrand_object]] class.
   !<
   !< @note It is assumed that the integrator has been properly initialized before calling this method.
   class(integrator_multistage_object), intent(inout) :: self !< Integrator.
   class(integrand_object),             intent(in)    :: U    !< Integrand.
   integer(I_P)                                       :: s    !< Counter.

   if (self%is_multistage() .and. self%registers > 0) then
      if (allocated(self%stage)) deallocate(self%stage)
      allocate(self%stage(1:self%registers), mold=U)
      do s=1, self%registers
         self%stage(s) = U
      enddo
   endif
   if (self%has_fast_mode()) then
      if (allocated(self%buffer)) deallocate(self%buffer)
      allocate(self%buffer, mold=U)
      self%buffer = U
   endif
   endsubroutine allocate_integrand_members

   pure subroutine assign_multistage(lhs, rhs)
   !< Assign members of [[integrator_multistage_object]] and parents.
   class(integrator_multistage_object), intent(inout) :: lhs !< Left hand side.
   class(integrator_object),            intent(in)    :: rhs !< Right hand side.
   integer(I_P)                                       :: s   !< Counter.

   call lhs%assign_abstract(rhs=rhs)
   select type(rhs)
   class is (integrator_multistage_object)
     lhs%registers = rhs%registers
     lhs%stages = rhs%stages
     if (allocated(lhs%stage)) deallocate(lhs%stage)
     if (allocated(rhs%stage)) then
        allocate(lhs%stage(1:lhs%registers), mold=rhs%stage)
        do s=1, lhs%registers
           lhs%stage(s) = rhs%stage(s)
        enddo
     endif
     if (allocated(lhs%buffer)) deallocate(lhs%buffer)
     if (allocated(rhs%buffer)) then
        allocate(lhs%buffer, mold=rhs%buffer)
        lhs%buffer = rhs%buffer
     endif
   endselect
   endsubroutine assign_multistage

   elemental subroutine destroy_multistage(self)
   !< Destroy the integrator.
   class(integrator_multistage_object), intent(inout) :: self !< Integrator.

   call self%destroy_abstract
   self%registers = 0
   self%stages = 0
   if (allocated(self%stage)) deallocate(self%stage)
   if (allocated(self%buffer)) deallocate(self%buffer)
   endsubroutine destroy_multistage
endmodule foodie_integrator_multistage_object
