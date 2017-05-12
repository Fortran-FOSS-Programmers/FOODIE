!< Define [[integrand_oscillation]], the Oscillation test field that is a concrete extension of the abstract integrand type.

module foodie_test_integrand_oscillation
!< Define [[integrand_oscillation]], the Oscillation test field that is a concrete extension of the abstract integrand type.

use foodie, only : integrand_object
use penf, only : R_P, I_P

implicit none
private
public :: integrand_oscillation

type, extends(integrand_object) :: integrand_oscillation
  !< The oscillation equations field.
  !<
  !< It is a FOODIE integrand class concrete extension.
  !<
  !<### Oscillation ODEs system
  !<The (inertial) Oscillation equations system is a non linear system of pure ODEs and it can be written as:
  !<
  !<$$\begin{matrix}
  !< U_t = R(U)  \\
  !< U = \begin{bmatrix}
  !< v_1 \\
  !< v_2
  !< \end{bmatrix}\;\;\;
  !< R(U) = \begin{bmatrix}
  !< -f v_2 \\
  !< f v_1
  !< \end{bmatrix}
  !<\end{matrix}$$
  !<
  !<The frequency *f* is constant and it is here selected as *f=10^-4*.
  !<
  !<In the space *v1-v2* the path of the oscillation must be a circle, but for the frequency selected some ODE solvers are not
  !<stable leading to a wrong path, see [1].
  !<
  !<#### Bibliography
  !<
  !<[1] *Numerical Methods for Fluid Dynamics With Applications to Geophysics*, Dale R. Durran, Springer, 2010.
  !<
  !<#### State variables organization
  !< State variables are organized as an array (rank 1) of reals of *n=2* elements.
  private
  real(R_P) :: f=0._R_P                 !< Oscillation frequency (Hz).
  real(R_P) :: U(1:2)=[0._R_P, 0._R_P]  !< Integrand (state) variables.
  real(R_P) :: U0(1:2)=[0._R_P, 0._R_P] !< Initial state.
  contains
     ! auxiliary methods
     procedure, pass(self), public :: initialize     !< Initialize integrand.
     procedure, pass(self), public :: exact_solution !< Return exact solution.
     procedure, pass(self), public :: output         !< Extract integrand state field.
     ! public deferred methods
     procedure, pass(self), public :: integrand_dimension !< Return integrand dimension.
     procedure, pass(self), public :: t => dU_dt          !< Time derivative, residuals.
     ! operators
     procedure, pass(lhs), public :: local_error !<`||integrand_oscillation - integrand_oscillation||` operator.
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
     ! override fast operators
     procedure, pass(self), public :: t_fast                              !< Time derivative, residuals, fast mode.
     procedure, pass(opr),  public :: integrand_add_integrand_fast        !< `+` fast operator.
     procedure, pass(opr),  public :: integrand_multiply_integrand_fast   !< `*` fast operator.
     procedure, pass(opr),  public :: integrand_multiply_real_scalar_fast !< `* real_scalar` fast operator.
     procedure, pass(opr),  public :: integrand_subtract_integrand_fast   !< `-` fast operator.
endtype integrand_oscillation

contains
   ! auxiliary methods
   pure subroutine initialize(self, U0, frequency)
   !< Initialize integrand.
   class(integrand_oscillation), intent(inout) :: self      !< Integrand.
   real(R_P),                    intent(in)    :: U0(1:2)   !< Initial state.
   real(R_P),                    intent(in)    :: frequency !< Frequency of oscillation.

   self%U = U0
   self%f = frequency
   self%U0 = U0
   endsubroutine initialize

   pure function exact_solution(self, t) result(exact)
   !< Return exact solution.
   class(integrand_oscillation), intent(in) :: self       !< Integrand.
   real(R_P),                    intent(in) :: t          !< Time.
   real(R_P)                                :: exact(1:2) !< Exact solution.

   exact(1) = self%U0(1) * cos(self%f * t) - self%U0(2) * sin(self%f * t)
   exact(2) = self%U0(1) * sin(self%f * t) + self%U0(2) * cos(self%f * t)
   endfunction exact_solution

   pure function output(self) result(state)
   !< Extract integrand state field.
   class(integrand_oscillation), intent(in) :: self       !< Integrand.
   real(R_P)                                :: state(1:2) !< State.

   state = self%U
   endfunction output

   ! deferred methods
   pure function integrand_dimension(self)
   !< return integrand dimension.
   class(integrand_oscillation), intent(in) :: self                !< integrand.
   integer(i_p)                             :: integrand_dimension !< integrand dimension.

   integrand_dimension = size(self%U, dim=1)
   endfunction integrand_dimension

   pure function dU_dt(self, t) result(dState_dt)
   !< Time derivative of integrand_oscillation field.
   class(integrand_oscillation), intent(in)           :: self         !< Integrand.
   real(R_P),                    intent(in), optional :: t            !< Time.
   real(R_P), allocatable                             :: dState_dt(:) !< Integrand time derivative.

   dState_dt = [-self%f * self%U(2), &
                 self%f * self%U(1)]
   endfunction dU_dt

   pure function local_error(lhs, rhs) result(error)
   !< Estimate local truncation error between 2 oscillation approximations.
   !<
   !< The estimation is done by norm L2 of U:
   !<
   !< $$ error = \sqrt{ \sum_i{ \frac{(lhs\%U_i - rhs\%U_i)^2}{lhs\%U_i^2} }} $$
   class(integrand_oscillation), intent(in) :: lhs   !< Left hand side.
   class(integrand_object),      intent(in) :: rhs   !< Right hand side.
   real(R_P)                                :: error !< Error estimation.
   integer(I_P)                             :: i     !< Space counter.

   select type(rhs)
   class is(integrand_oscillation)
      error = 0._R_P
      do i=1, size(lhs%U, dim=1)
         error = error + (lhs%U(i) - rhs%U(i)) ** 2 / lhs%U(i) ** 2
      enddo
      error = sqrt(error)
   endselect
   endfunction local_error

   ! +
   pure function integrand_add_integrand(lhs, rhs) result(opr)
   !< `+` operator.
   class(integrand_oscillation), intent(in) :: lhs    !< Left hand side.
   class(integrand_object),      intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                   :: opr(:) !< Operator result.

   select type(rhs)
   class is(integrand_oscillation)
      opr = lhs%U + rhs%U
   endselect
   endfunction integrand_add_integrand

   pure function integrand_add_real(lhs, rhs) result(opr)
   !< `+ real` operator.
   class(integrand_oscillation), intent(in) :: lhs     !< Left hand side.
   real(R_P),                    intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable                   :: opr(:)  !< Operator result.

   opr = lhs%U + rhs
   endfunction integrand_add_real

   pure function real_add_integrand(lhs, rhs) result(opr)
   !< `real +` operator.
   real(R_P),                    intent(in) :: lhs(1:) !< Left hand side.
   class(integrand_oscillation), intent(in) :: rhs     !< Left hand side.
   real(R_P), allocatable                   :: opr(:)  !< Operator result.

   opr = lhs + rhs%U
   endfunction real_add_integrand

   ! *
   pure function integrand_multiply_integrand(lhs, rhs) result(opr)
   !< `*` operator.
   class(integrand_oscillation), intent(in) :: lhs    !< Left hand side.
   class(integrand_object),      intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                   :: opr(:) !< Operator result.

   select type(rhs)
   class is(integrand_oscillation)
      opr = lhs%U * rhs%U
   endselect
   endfunction integrand_multiply_integrand

   pure function integrand_multiply_real(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(integrand_oscillation), intent(in) :: lhs     !< Left hand side.
   real(R_P),                    intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable                   :: opr(:)  !< Operator result.

   opr = lhs%U * rhs
   endfunction integrand_multiply_real

   pure function real_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   class(integrand_oscillation), intent(in) :: rhs     !< Right hand side.
   real(R_P),                    intent(in) :: lhs(1:) !< Left hand side.
   real(R_P), allocatable                   :: opr(:)  !< Operator result.

   opr = lhs * rhs%U
   endfunction real_multiply_integrand

   pure function integrand_multiply_real_scalar(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(integrand_oscillation), intent(in) :: lhs    !< Left hand side.
   real(R_P),                    intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                   :: opr(:) !< Operator result.

   opr = lhs%U * rhs
   endfunction integrand_multiply_real_scalar

   pure function real_scalar_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   real(R_P),                    intent(in) :: lhs    !< Left hand side.
   class(integrand_oscillation), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                   :: opr(:) !< Operator result.

   opr = lhs * rhs%U
   endfunction real_scalar_multiply_integrand

   ! -
   pure function integrand_sub_integrand(lhs, rhs) result(opr)
   !< `-` operator.
   class(integrand_oscillation), intent(in) :: lhs    !< Left hand side.
   class(integrand_object),      intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                   :: opr(:) !< Operator result.

   select type(rhs)
   class is(integrand_oscillation)
      opr = lhs%U - rhs%U
   endselect
   endfunction integrand_sub_integrand

   pure function integrand_sub_real(lhs, rhs) result(opr)
   !< `- real` operator.
   class(integrand_oscillation), intent(in) :: lhs     !< Left hand side.
   real(R_P),                    intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable                   :: opr(:)  !< Operator result.

   opr = lhs%U - rhs
   endfunction integrand_sub_real

   pure function real_sub_integrand(lhs, rhs) result(opr)
   !< `real -` operator.
   class(integrand_oscillation), intent(in) :: rhs     !< Left hand side.
   real(R_P),                    intent(in) :: lhs(1:) !< Left hand side.
   real(R_P), allocatable                   :: opr(:)  !< Operator result.

   opr = lhs - rhs%U
   endfunction real_sub_integrand

   ! =
   pure subroutine assign_integrand(lhs, rhs)
   !< `=` operator.
   class(integrand_oscillation), intent(inout) :: lhs !< Left hand side.
   class(integrand_object),      intent(in)    :: rhs !< Right hand side.

   select type(rhs)
   class is(integrand_oscillation)
      lhs%U = rhs%U
      lhs%f = rhs%f
      lhs%U0 = rhs%U0
   endselect
   endsubroutine assign_integrand

   pure subroutine assign_real(lhs, rhs)
   !< `= real` operator.
   class(integrand_oscillation), intent(inout) :: lhs     !< Left hand side.
   real(R_P),                    intent(in)    :: rhs(1:) !< Right hand side.

   lhs%U = rhs
   endsubroutine assign_real

   ! fast operators
   ! time derivative
   subroutine t_fast(self, t)
   !< Time derivative function of integrand class, i.e. the residuals function. Fast mode acting directly on self.
   class(integrand_oscillation), intent(inout)        :: self   !< Oscillation field.
   real(R_P),                    intent(in), optional :: t      !< Time.
   real(R_P)                                          :: buffer !< Temporary buffer.

   buffer = self%U(1)
   self%U(1) = - self%f * self%U(2)
   self%U(2) =   self%f * buffer
   endsubroutine t_fast

   ! +
   pure subroutine integrand_add_integrand_fast(opr, lhs, rhs)
   !< `+` fast operator.
   class(integrand_oscillation), intent(inout) :: opr !< Operator result.
   class(integrand_object),      intent(in)    :: lhs !< Left hand side.
   class(integrand_object),      intent(in)    :: rhs !< Right hand side.

   select type(lhs)
   class is(integrand_oscillation)
      select type(rhs)
      class is(integrand_oscillation)
         opr%U = lhs%U + rhs%U
      endselect
   endselect
   endsubroutine integrand_add_integrand_fast

   ! *
   pure subroutine integrand_multiply_integrand_fast(opr, lhs, rhs)
   !< `*` fast operator.
   class(integrand_oscillation), intent(inout) :: opr !< Operator result.
   class(integrand_object),      intent(in)    :: lhs !< Left hand side.
   class(integrand_object),      intent(in)    :: rhs !< Right hand side.

   select type(lhs)
   class is(integrand_oscillation)
      select type(rhs)
      class is(integrand_oscillation)
         opr%U = lhs%U * rhs%U
      endselect
   endselect
   endsubroutine integrand_multiply_integrand_fast

   pure subroutine integrand_multiply_real_scalar_fast(opr, lhs, rhs)
   !< `* real_scalar` fast operator.
   class(integrand_oscillation), intent(inout) :: opr !< Operator result.
   class(integrand_object),      intent(in)    :: lhs !< Left hand side.
   real(R_P),                    intent(in)    :: rhs !< Right hand side.

   select type(lhs)
   class is(integrand_oscillation)
      opr%U = lhs%U * rhs
   endselect
   endsubroutine integrand_multiply_real_scalar_fast

   ! -
   pure subroutine integrand_subtract_integrand_fast(opr, lhs, rhs)
   !< `-` fast operator.
   class(integrand_oscillation), intent(inout) :: opr !< Operator result.
   class(integrand_object),      intent(in)    :: lhs !< Left hand side.
   class(integrand_object),      intent(in)    :: rhs !< Right hand side.

   select type(lhs)
   class is(integrand_oscillation)
      select type(rhs)
      class is(integrand_oscillation)
         opr%U = lhs%U - rhs%U
      endselect
   endselect
   endsubroutine integrand_subtract_integrand_fast
endmodule foodie_test_integrand_oscillation
