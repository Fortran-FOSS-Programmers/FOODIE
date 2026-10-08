module my_integrand
use foodie, only : integrand_object
use penf, only : I_P, R_P

implicit none
private
public :: integrand_lcce

type, extends(integrand_object) :: integrand_lcce
   !< The linear constant coefficient equation field, U_t = a*U + b.
   real(R_P) :: a=0._R_P  !< *a* constant.
   real(R_P) :: b=0._R_P  !< *b* constant.
   real(R_P) :: U=0._R_P  !< Integrand (state) variable.
   contains
      ! integrand_object deferred methods
      procedure, pass(self), public :: description         !< Return an informative description of the integrand.
      procedure, pass(self), public :: integrand_dimension !< Return integrand dimension.
      procedure, pass(self), public :: t => dU_dt          !< Time derivative, residuals.
      ! operators
      procedure, pass(lhs), public :: local_error !<`||integrand_lcce - integrand_lcce||` operator.
      ! +
      procedure, pass(lhs), public :: integrand_add_integrand !< `+` operator.
      procedure, pass(lhs), public :: integrand_add_real      !< `+ real` operator.
      procedure, pass(rhs), public :: real_add_integrand      !< `real +` operator.
      ! *
      procedure, pass(lhs), public :: integrand_multiply_integrand   !< `*` operator.
      procedure, pass(lhs), public :: integrand_multiply_real        !< `* real` operator.
      procedure, pass(rhs), public :: real_multiply_integrand        !< `real *` operator.
      procedure, pass(lhs), public :: integrand_multiply_real_scalar !< `* real_scalar` operator.
      procedure, pass(rhs), public :: real_scalar_multiply_integrand !< `real_scalar *` operator.
      ! -
      procedure, pass(lhs), public :: integrand_sub_integrand !< `-` operator.
      procedure, pass(lhs), public :: integrand_sub_real      !< `- real` operator.
      procedure, pass(rhs), public :: real_sub_integrand      !< `real -` operator.
      ! =
      procedure, pass(lhs), public :: assign_integrand !< `=` operator.
      procedure, pass(lhs), public :: assign_real      !< `= real` operator.
endtype integrand_lcce

contains
   pure function description(self, prefix) result(desc)
   !< Return informative integrator description.
   class(integrand_lcce), intent(in)           :: self    !< Integrand.
   character(*),          intent(in), optional :: prefix  !< Prefixing string.
   character(len=:), allocatable               :: desc    !< Description.
   character(len=:), allocatable               :: prefix_ !< Prefixing string, local variable.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = prefix_//'linear_constant_coefficients_eq'
   endfunction description

   pure function integrand_dimension(self)
   !< return integrand dimension.
   class(integrand_lcce), intent(in) :: self                !< integrand.
   integer(I_P)                      :: integrand_dimension !< integrand dimension.

   integrand_dimension = 1
   endfunction integrand_dimension

   pure function dU_dt(self, t) result(dState_dt)
   !< Time derivative of field.
   class(integrand_lcce), intent(in)           :: self         !< Integrand.
   real(R_P),             intent(in), optional :: t            !< Time.
   real(R_P), allocatable                      :: dState_dt(:) !< Integrand time derivative.

   dState_dt = [self%a * self%U + self%b]
   endfunction dU_dt

   pure function local_error(lhs, rhs) result(error)
   !< Estimate local truncation error between 2 lcce approximations.
   class(integrand_lcce),   intent(in) :: lhs   !< Left hand side.
   class(integrand_object), intent(in) :: rhs   !< Right hand side.
   real(R_P)                           :: error !< Error estimation.

   select type(rhs)
   class is(integrand_lcce)
      if (lhs%U /= 0._R_P) then
         error = sqrt(((lhs%U - rhs%U) ** 2) / lhs%U **2)
      elseif (rhs%U /= 0._R_P) then
         error = sqrt(((lhs%U - rhs%U) ** 2) / rhs%U **2)
      else
         error = 0._R_P
      endif
   endselect
   endfunction local_error

   pure function integrand_add_integrand(lhs, rhs) result(opr)
   !< `+` operator.
   class(integrand_lcce),   intent(in) :: lhs    !< Left hand side.
   class(integrand_object), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable              :: opr(:) !< Operator result.

   select type(rhs)
   class is(integrand_lcce)
     opr = [lhs%U + rhs%U]
   endselect
   endfunction integrand_add_integrand

   pure function integrand_add_real(lhs, rhs) result(opr)
   !< `+ real` operator.
   class(integrand_lcce), intent(in) :: lhs     !< Left hand side.
   real(R_P),             intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable            :: opr(:)  !< Operator result.

   opr = [lhs%U + rhs(1)]
   endfunction integrand_add_real

   pure function real_add_integrand(lhs, rhs) result(opr)
   !< `real +` operator.
   real(R_P),             intent(in) :: lhs(1:) !< Left hand side.
   class(integrand_lcce), intent(in) :: rhs     !< Left hand side.
   real(R_P), allocatable            :: opr(:)  !< Operator result.

   opr = [lhs(1) + rhs%U]
   endfunction real_add_integrand

   pure function integrand_multiply_integrand(lhs, rhs) result(opr)
   !< `*` operator.
   class(integrand_lcce),   intent(in) :: lhs    !< Left hand side.
   class(integrand_object), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable              :: opr(:) !< Operator result.

   select type(rhs)
   class is(integrand_lcce)
     opr = [lhs%U * rhs%U]
   endselect
   endfunction integrand_multiply_integrand

   pure function integrand_multiply_real(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(integrand_lcce), intent(in) :: lhs     !< Left hand side.
   real(R_P),             intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable            :: opr(:)  !< Operator result.

   opr = [lhs%U * rhs(1)]
   endfunction integrand_multiply_real

   pure function real_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   class(integrand_lcce), intent(in) :: rhs     !< Right hand side.
   real(R_P),             intent(in) :: lhs(1:) !< Left hand side.
   real(R_P), allocatable            :: opr(:)  !< Operator result.

   opr = [lhs(1) * rhs%U]
   endfunction real_multiply_integrand

   pure function integrand_multiply_real_scalar(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(integrand_lcce), intent(in) :: lhs    !< Left hand side.
   real(R_P),             intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable            :: opr(:) !< Operator result.

   opr = [lhs%U * rhs]
   endfunction integrand_multiply_real_scalar

   pure function real_scalar_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   real(R_P),             intent(in) :: lhs    !< Left hand side.
   class(integrand_lcce), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable            :: opr(:) !< Operator result.

   opr = [lhs * rhs%U]
   endfunction real_scalar_multiply_integrand

   pure function integrand_sub_integrand(lhs, rhs) result(opr)
   !< `-` operator.
   class(integrand_lcce),   intent(in) :: lhs    !< Left hand side.
   class(integrand_object), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable              :: opr(:) !< Operator result.

   select type(rhs)
   class is(integrand_lcce)
     opr = [lhs%U - rhs%U]
   endselect
   endfunction integrand_sub_integrand

   pure function integrand_sub_real(lhs, rhs) result(opr)
   !< `- real` operator.
   class(integrand_lcce), intent(in) :: lhs     !< Left hand side.
   real(R_P),             intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable            :: opr(:)  !< Operator result.

   opr = [lhs%U - rhs(1)]
   endfunction integrand_sub_real

   pure function real_sub_integrand(lhs, rhs) result(opr)
   !< `real -` operator.
   real(R_P),             intent(in) :: lhs(1:) !< Left hand side.
   class(integrand_lcce), intent(in) :: rhs     !< Left hand side.
   real(R_P), allocatable            :: opr(:)  !< Operator result.

   opr = [lhs(1) - rhs%U]
   endfunction real_sub_integrand

   pure subroutine assign_integrand(lhs, rhs)
   !< `=` operator.
   class(integrand_lcce),   intent(inout) :: lhs !< Left hand side.
   class(integrand_object), intent(in)    :: rhs !< Right hand side.

   select type(rhs)
   class is(integrand_lcce)
      lhs%a  = rhs%a
      lhs%b  = rhs%b
      lhs%U  = rhs%U
   endselect
   endsubroutine assign_integrand

   pure subroutine assign_real(lhs, rhs)
   !< `= real` operator.
   class(integrand_lcce), intent(inout) :: lhs     !< Left hand side.
   real(R_P),             intent(in)    :: rhs(1:) !< Right hand side.

   lhs%U = rhs(1)
   endsubroutine assign_real
endmodule my_integrand

program quickstart
use foodie, only : foodie_integrator_factory, integrator_object, integrator_multistage_object, integrator_multistep_object
use my_integrand, only : integrand_lcce
use penf, only : I_P, R_P

implicit none
type(integrand_lcce)                  :: integrand  !< The integrand.
class(integrator_object), allocatable :: integrator !< The integrator.
real(R_P), parameter                  :: Dt=0.01_R_P !< Time step.
real(R_P)                             :: time        !< Time.
integer(I_P)                          :: step        !< Time steps counter.

integrand%a = -1._R_P ; integrand%b = 0._R_P ; integrand%U = 1._R_P
time = 0._R_P

call foodie_integrator_factory(scheme='runge_kutta_ssp_stages_3_order_3', integrator=integrator, U=integrand)
select type(integrator)
class is(integrator_multistage_object)
   do step=1, 100
      call integrator%integrate(U=integrand, Dt=Dt, t=time)
      time = time + Dt
   enddo
endselect
print '(A,ES23.15,A,ES23.15)', 'U(1) = ', integrand%U, ', exact ', exp(-1._R_P)

! multistep: start-up with a one-step scheme
integrand%U = 1._R_P ; time = 0._R_P
block
class(integrator_object), allocatable :: starter
call foodie_integrator_factory(scheme='runge_kutta_ssp_stages_3_order_3', integrator=starter, U=integrand)
call foodie_integrator_factory(scheme='adams_bashforth_3', integrator=integrator, U=integrand)
select type(integrator)
class is(integrator_multistep_object)
   do step=1, 100
      if (step <= integrator%steps_number()) then
         select type(starter)
         class is(integrator_multistage_object)
            call starter%integrate(U=integrand, Dt=Dt, t=time)
         endselect
         time = time + Dt
         integrator%previous(step) = integrand
         integrator%Dt(step) = Dt
         integrator%t(step) = time
      else
         call integrator%integrate(U=integrand, Dt=Dt, t=time)
         time = time + Dt
      endif
   enddo
endselect
endblock
print '(A,ES23.15,A,ES23.15)', 'U(1) = ', integrand%U, ', exact ', exp(-1._R_P)
endprogram quickstart
