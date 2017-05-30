!< Define [[integrand_oscillation]], the Oscillation test field that is a concrete extension of the abstract integrand type.

module foodie_test_integrand_oscillation
!< Define [[integrand_oscillation]], the Oscillation test field that is a concrete extension of the abstract integrand type.

use flap, only : command_line_interface
use foodie, only : integrand_object
use foodie_test_integrand_tester_object, only : integrand_tester_object
use penf, only : FR_P, R_P, I_P, str

implicit none
private
public :: integrand_oscillation

type, extends(integrand_tester_object) :: integrand_oscillation
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
   real(R_P) :: U0(1:2)=[0._R_P, 0._R_P] !< Initial state.
   real(R_P) :: U(1:2)=[0._R_P, 0._R_P]  !< Integrand (state) variables.
   contains
      ! auxiliary methods
      procedure, pass(self), public :: amplitude_phase !< Return amplitude and phase of the oscillation.
      procedure, pass(self), public :: output          !< Extract integrand state field.
      ! integrand_tester_object deferred methods
      procedure, pass(self), public :: description    !< Return an informative description of the test.
      procedure, pass(self), public :: error          !< Return error.
      procedure, pass(self), public :: exact_solution !< Return exact solution.
      procedure, pass(self), public :: export_tecplot !< Export integrand to Tecplot file.
      procedure, pass(self), public :: initialize     !< Initialize field.
      procedure, pass(self), public :: parse_cli      !< Initialize from command line interface.
      procedure, nopass,     public :: set_cli        !< Set command line interface.
      ! integrand_object deferred methods
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
   function amplitude_phase(self) result(ap)
   !< Compute amplitude and phase of the oscillation.
   class(integrand_oscillation), intent(in) :: self    !< Advection field.
   real(R_P)                                :: ap(1:2) !< Amplitude and phase.

   ap(1) = sqrt(self%U(1)**2 + self%U(2)**2)
   ap(2) = atan(-self%U(1) / self%U(2))
   endfunction amplitude_phase

   pure function output(self) result(state)
   !< Extract integrand state field.
   class(integrand_oscillation), intent(in) :: self       !< Integrand.
   real(R_P)                                :: state(1:2) !< State.

   state = self%U
   endfunction output

   ! integrand_tester_object deferred methods
   pure function description(self, prefix) result(desc)
   !< Return informative integrator description.
   class(integrand_oscillation), intent(in)           :: self    !< Integrand.
   character(*),                 intent(in), optional :: prefix  !< Prefixing string.
   character(len=:), allocatable                      :: desc    !< Description.
   character(len=:), allocatable                      :: prefix_ !< Prefixing string, local variable.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = prefix//'oscillation'
   endfunction description

   pure function error(self, t, t0, U0)
   !< Return error.
   class(integrand_oscillation), intent(in)           :: self     !< Integrand.
   real(R_P),                    intent(in)           :: t        !< Time.
   real(R_P),                    intent(in), optional :: t0       !< Initial time.
   class(integrand_object),      intent(in), optional :: U0       !< Initial conditions.
   real(R_P), allocatable                             :: error(:) !< Error.

   allocate(error(1:2))
   error = abs(self%U - self%exact_solution(t=t))
   endfunction error

   pure function exact_solution(self, t, t0, U0) result(exact)
   !< Return exact solution.
   class(integrand_oscillation), intent(in)           :: self     !< Integrand.
   real(R_P),                    intent(in)           :: t        !< Time.
   real(R_P),                    intent(in), optional :: t0       !< Initial time.
   class(integrand_object),      intent(in), optional :: U0       !< Initial conditions.
   real(R_P), allocatable                             :: exact(:) !< Exact solution.

   exact = [self%U0(1) * cos(self%f * t) - self%U0(2) * sin(self%f * t), &
            self%U0(1) * sin(self%f * t) + self%U0(2) * cos(self%f * t)]
   endfunction exact_solution

   subroutine export_tecplot(self, file_name, t, scheme, close_file)
   !< Export integrand to Tecplot file.
   class(integrand_oscillation), intent(in)           :: self            !< Advection field.
   character(*),                 intent(in), optional :: file_name       !< File name.
   real(R_P),                    intent(in), optional :: t               !< Time.
   character(*),                 intent(in), optional :: scheme          !< Scheme used to integrate integrand.
   logical,                      intent(in), optional :: close_file      !< Flag for closing file.
   logical, save                                      :: is_open=.false. !< Flag for checking if file is open.
   integer(I_P), save                                 :: file_unit       !< File unit.

   if (present(close_file)) then
      if (close_file .and. is_open) then
         close(unit=file_unit)
         is_open = .false.
      endif
   else
      if (present(file_name)) then
         if (is_open) close(unit=file_unit)
         open(newunit=file_unit, file=trim(adjustl(file_name)))
         is_open = .true.
         write(unit=file_unit, fmt='(A)') 'VARIABLES="t" "x" "y" "amplitude" "phase"'
      endif
      if (present(t) .and. present(scheme) .and. is_open) then
         write(unit=file_unit, fmt='(A)') 'ZONE T="'//trim(adjustl(scheme))//'"'
         write(unit=file_unit, fmt='(5('//FR_P//',1X))') t, self%U, self%amplitude_phase()
      elseif (present(t) .and. is_open) then
         write(unit=file_unit, fmt='(5('//FR_P//',1X))') t, self%U, self%amplitude_phase()
      endif
   endif
   endsubroutine export_tecplot

   pure subroutine initialize(self, Dt)
   !< Initialize integrand.
   !<
   !< Intentionally empty, all is done in `parse_cli` method.
   class(integrand_oscillation), intent(inout) :: self !< Integrand.
   real(R_P),                    intent(in)    :: Dt   !< Time step.
   endsubroutine initialize

   subroutine parse_cli(self, cli)
   !< Initialize from command line interface.
   class(integrand_oscillation), intent(inout) :: self !< Advection field.
   type(command_line_interface), intent(inout) :: cli  !< Command line interface handler.

   call cli%get(group='oscillation', switch='-f', val=self%f, error=cli%error) ; if (cli%error/=0) stop
   call cli%get(group='oscillation', switch='-U0', val=self%U0, error=cli%error) ; if (cli%error/=0) stop
   self%U = self%U0
   endsubroutine parse_cli

   subroutine set_cli(cli)
   !< Set command line interface.
   type(command_line_interface), intent(inout) :: cli !< Command line interface handler.

   call cli%add_group(description='oscillation test settings', group='oscillation')
   call cli%add(group='oscillation', switch='--frequency', switch_ab='-f', help='frequency', required=.false., def='1e-4', &
                act='store')
   call cli%add(group='oscillation', switch='--U0', switch_ab='-U0', nargs='2', help='initial state', required=.false., &
                def='0.0 1.0', act='store')
   endsubroutine set_cli

   ! integrand_object deferred methods
   pure function integrand_dimension(self)
   !< return integrand dimension.
   class(integrand_oscillation), intent(in) :: self                !< integrand.
   integer(I_P)                             :: integrand_dimension !< integrand dimension.

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
      lhs%f = rhs%f
      lhs%U0 = rhs%U0
      lhs%U = rhs%U
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
