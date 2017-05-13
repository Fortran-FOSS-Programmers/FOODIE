!< Define [[integrand_ladvection]], the 1D linear advection PDE test field that is a concrete extension of the
!< abstract integrand type.

module foodie_test_integrand_ladvection
!< Define [[integrand_ladvection]], the 1D linear advection PDE test field that is a concrete extension of the
!< abstract integrand type.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foodie, only : integrand_object
use penf, only : I_P, R_P
use wenoof, only : interpolator_object, wenoof_create

implicit none
private
public :: integrand_ladvection

type, extends(integrand_object) :: integrand_ladvection
   !< 1D linear advection field.
   !<
   !< It is a FOODIE integrand class concrete extension.
   !<
   !<### 1D linear advection field
   !< The 1D linear advection equation is a conservation law that reads as
   !<$$
   !<\begin{matrix}
   !<u_t = R(u)  \Leftrightarrow u_t = F(u)_x \\
   !<F(u) = a * u
   !<$$
   !< where `a` is scalar constant coefficient. The PDE must completed with the proper initial and boundary conditions.
   !<
   !<#### Numerical grid organization
   !< The finite volume, Godunov's like approach is employed. The conservative variables (and the primitive ones) are co-located at
   !< the cell center. The cell and (inter)faces numeration is as follow.
   !<```
   !<                cell            (inter)faces
   !<                 |                   |
   !<                 v                   v
   !<     |-------|-------|-.....-|-------|-------|-------|-------|-.....-|-------|-------|-------|-.....-|-------|-------|
   !<     | 1-Ng  | 2-Ng  | ..... |  -1   |   0   |   1   |  2    | ..... |  Ni   | Ni+1  | Ni+1  | ..... |Ni+Ng-1| Ni+Ng |
   !<     |-------|-------|-.....-|-------|-------|-------|-------|-.....-|-------|-------|-------|-.....-|-------|-------|
   !<    0-Ng                             -1      0       1       2      Ni-1     Ni                                    Ni+Ng
   !<```
   !< Where *Ni* are the finite volumes (cells) used for discretizing the domain and *Ng* are the ghost cells used for imposing the
   !< left and right boundary conditions (for a total of *2Ng* cells).
   integer(I_P)                            :: weno_order=0           !< WENO reconstruction order.
   real(R_P)                               :: weno_eps=10._R_P**(-6) !< WENO epsilon to avoid division by zero, default value.
   integer(I_P)                            :: Ni=0                   !< Space dimension.
   integer(I_P)                            :: Ng=0                   !< Ghost cells number.
   real(R_P)                               :: Dx=0._R_P              !< Space step.
   real(R_P)                               :: a=0._R_P               !< Advection coefficient.
   real(R_P), allocatable                  :: u(:)                   !< Integrand (state) variable.
   character(:), allocatable               :: BC_L                   !< Left boundary condition type.
   character(:), allocatable               :: BC_R                   !< Right boundary condition type.
   class(interpolator_object), allocatable :: interpolator           !< WENO interpolator.
   contains
      ! auxiliary methods
      procedure, pass(self), public :: destroy          !< Destroy field.
      procedure, pass(self), public :: dt => compute_dt !< Compute the current time step, by means of CFL condition.
      procedure, pass(self), public :: exact_solution   !< Return exact solution.
      procedure, pass(self), public :: initialize       !< Initialize integrand.
      procedure, pass(self), public :: output           !< Extract integrand state field.
      ! public deferred methods
      procedure, pass(self), public :: integrand_dimension !< Return integrand dimension.
      procedure, pass(self), public :: t => dU_dt          !< Time derivative, residuals.
      ! operators
      procedure, pass(lhs), public :: local_error !<`||integrand_ladvection - integrand_ladvection||` operator.
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
      ! private methods
      procedure, pass(self), private :: impose_boundary_conditions !< Impose boundary conditions.
      procedure, pass(self), private :: reconstruct_interfaces     !< Reconstruct interface states.
endtype integrand_ladvection

contains
   ! auxiliary methods
   pure subroutine destroy(self)
   !< Destroy field.
   class(integrand_ladvection), intent(inout) :: self !< Advection field.

   self%weno_order = 0
   self%Ni = 0
   self%Ng = 0
   self%Dx = 0._R_P
   if (allocated(self%u)) deallocate(self%u)
   if (allocated(self%BC_L)) deallocate(self%BC_L)
   if (allocated(self%BC_R)) deallocate(self%BC_R)
   if (allocated(self%interpolator)) deallocate(self%interpolator)
   endsubroutine destroy

   pure function compute_dt(self, steps_max, t_max, t, CFL) result(Dt)
   !< Compute the current time step by means of CFL condition.
   class(integrand_ladvection), intent(in) :: self      !< Advection field.
   integer(I_P),                intent(in) :: steps_max !< Maximun number of time steps.
   real(R_P),                   intent(in) :: t_max     !< Maximum integration time.
   real(R_P),                   intent(in) :: t         !< Time.
   real(R_P),                   intent(in) :: CFL       !< CFL value.
   real(R_P)                               :: Dt        !< Time step.

   associate(Ni=>self%Ni, Dx=>self%Dx)
      Dt = Dx * CFL / abs(self%a)
      if (steps_max <= 0 .and. t_max > 0._R_P) then
         if ((t + Dt) > t_max) Dt = t_max - t
      endif
   endassociate
   endfunction compute_dt

   pure function exact_solution(self, u0, t) result(exact)
   !< Return exact solution.
   class(integrand_ladvection), intent(in) :: self     !< Integrand.
   real(R_P),                   intent(in) :: u0(1:)   !< Initial state
   real(R_P),                   intent(in) :: t        !< Time.
   real(R_P), allocatable                  :: exact(:) !< Exact solution.

   endfunction exact_solution

   subroutine initialize(self, Ni, Dx, BC_L, BC_R, initial_state, w_scheme, a, weno_order)
   !< Initialize field.
   class(integrand_ladvection), intent(inout)        :: self              !< Advection field.
   integer(I_P),                intent(in)           :: Ni                !< Space dimension.
   real(R_P),                   intent(in)           :: Dx                !< Space step.
   character(*),                intent(in)           :: BC_L              !< Left boundary condition type.
   character(*),                intent(in)           :: BC_R              !< Right boundary condition type.
   real(R_P),                   intent(in)           :: initial_state(1:) !< Initial state.
   character(*),                intent(in)           :: w_scheme          !< WENO scheme.
   real(R_P),                   intent(in), optional :: a                 !< Advection coefficient.
   integer(I_P),                intent(in), optional :: weno_order        !< WENO reconstruction order.
   integer(I_P)                                      :: i                 !< Space couner.

   call self%destroy
   self%a = 1._R_P ; if (present(a)) self%a = a
   self%weno_order = 1 ; if (present(weno_order)) self%weno_order = weno_order
   self%Ni = Ni
   self%Ng = (self%weno_order + 1) / 2
   self%Dx = Dx
   if (allocated(self%u)) deallocate(self%u) ; allocate(self%u(1-self%Ng:self%Ni+self%Ng))
   do i=1, Ni
      self%u(i) = initial_state(i)
   enddo
   self%BC_L = BC_L
   self%BC_R = BC_R

   if (self%weno_order>1) call wenoof_create(interpolator_type=trim(adjustl(w_scheme)), S=self%Ng, interpolator=self%interpolator, &
                                             eps=self%weno_eps)
   endsubroutine initialize

   pure function output(self) result(state)
   !< Output the advection field state.
   class(integrand_ladvection), intent(in) :: self     !< Advection field.
   real(R_P), allocatable                  :: state(:) !< Advection state

   state = self%u(1:self%Ni)
   endfunction output

   ! ADT integrand deferred methods
   function dU_dt(self, t) result(dState_dt)
   !< Time derivative of advection field, the residuals function.
   class(integrand_ladvection), intent(in)           :: self                         !< Advection field.
   real(R_P),                   intent(in), optional :: t                            !< Time.
   real(R_P), allocatable                            :: dState_dt(:)                 !< Advection field time derivative.
   real(R_P)                                         :: u(1-self%Ng:self%Ni+self%Ng) !< Conservative variable.
   real(R_P)                                         :: ur(1:2,0:self%Ni+1)          !< Reconstructed conservative variable.
   real(R_P)                                         :: f(0:self%Ni)                 !< Flux of conservative variable.
   integer(I_P)                                      :: i                            !< Counter.

   do i=1, self%Ni
      U(i) = self%U(i)
   enddo
   call self%impose_boundary_conditions(u=u)
   call self%reconstruct_interfaces(conservative=u, r_conservative=ur)
   do i=0, self%Ni
      call solve_riemann_problem(state_left=ur(2, i), state_right=ur(1, i+1), flux=f(i))
   enddo
   allocate(dState_dt(1:self%Ni))
   do i=1, self%Ni
       dState_dt(i) = (f(i - 1) - f(i)) / self%Dx
   enddo
   contains
      subroutine solve_riemann_problem(state_left, state_right, flux)
      !< Solver Riemann problem of linear advection by upwinding.
      real(R_P), intent(in)  :: state_left  !< Left state.
      real(R_P), intent(in)  :: state_right !< right state.
      real(R_P), intent(out) :: flux        !< Flux of conservative variable.

      if (self%a > 0._R_P) then
         flux = self%a * state_left
      else
         flux = self%a * state_right
      endif
      endsubroutine solve_riemann_problem
   endfunction dU_dt

   pure function integrand_dimension(self)
   !< return integrand dimension.
   class(integrand_ladvection), intent(in) :: self                !< integrand.
   integer(I_P)                            :: integrand_dimension !< integrand dimension.

   integrand_dimension = self%Ni
   endfunction integrand_dimension

   function local_error(lhs, rhs) result(error)
   !< Estimate local truncation error between 2 advection approximations.
   !<
   !< The estimation is done by norm L2 of U:
   !<
   !< $$ error = \sqrt{ \sum_i{\sum_i{ \frac{(lhs\%u_i - rhs\%u_i)^2}{lhs\%u_i^2} }} } $$
   class(integrand_ladvection), intent(in) :: lhs   !< Left hand side.
   class(integrand),               intent(in) :: rhs   !< Right hand side.
   real(R_P)                                  :: error !< Error estimation.
   integer(I_P)                               :: i     !< Space counter.

   select type(rhs)
   class is (integrand_ladvection)
      error = 0._R_P
      do i=1, lhs%Ni
         error = error + (lhs%u(i) - rhs%u(i)) ** 2 / lhs%u(i) ** 2
      enddo
      error = sqrt(error)
   endselect
   endfunction advection_local_error

   ! +
   pure function integrand_add_integrand(lhs, rhs) result(opr)
   !< `+` operator.
   class(integrand_ladvection), intent(in) :: lhs    !< Left hand side.
   class(integrand_object),     intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                  :: opr(:) !< Operator result.

   select type(rhs)
   class is(integrand_ladvection)
     opr = lhs%u + rhs%u
   endselect
   endfunction integrand_add_integrand

   pure function integrand_add_real(lhs, rhs) result(opr)
   !< `+ real` operator.
   class(integrand_ladvection), intent(in) :: lhs     !< Left hand side.
   real(R_P),                   intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable                  :: opr(:)  !< Operator result.

   opr = lhs%u + rhs
   endfunction integrand_add_real

   pure function real_add_integrand(lhs, rhs) result(opr)
   !< `real +` operator.
   real(R_P),                   intent(in) :: lhs(1:) !< Left hand side.
   class(integrand_ladvection), intent(in) :: rhs     !< Left hand side.
   real(R_P), allocatable                  :: opr(:)  !< Operator result.

   opr = lhs + rhs%U
   endfunction real_add_integrand

   ! *
   pure function integrand_multiply_integrand(lhs, rhs) result(opr)
   !< `*` operator.
   class(integrand_ladvection), intent(in) :: lhs    !< Left hand side.
   class(integrand_object),     intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                  :: opr(:) !< Operator result.

   select type(rhs)
   class is(integrand_ladvection)
     opr = lhs%U * rhs%U
   endselect
   endfunction integrand_multiply_integrand

   pure function integrand_multiply_real(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(integrand_ladvection), intent(in) :: lhs     !< Left hand side.
   real(R_P),                   intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable                  :: opr(:)  !< Operator result.

   opr = lhs%U * rhs
   endfunction integrand_multiply_real

   pure function real_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   class(integrand_ladvection), intent(in) :: rhs     !< Right hand side.
   real(R_P),                   intent(in) :: lhs(1:) !< Left hand side.
   real(R_P), allocatable                  :: opr(:)  !< Operator result.

   opr = lhs * rhs%U
   endfunction real_multiply_integrand

   pure function integrand_multiply_real_scalar(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(integrand_ladvection), intent(in) :: lhs    !< Left hand side.
   real(R_P),                   intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                  :: opr(:) !< Operator result.

   opr = lhs%U * rhs
   endfunction integrand_multiply_real_scalar

   pure function real_scalar_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   real(R_P),                   intent(in) :: lhs    !< Left hand side.
   class(integrand_ladvection), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                  :: opr(:) !< Operator result.

   opr = lhs * rhs%U
   endfunction real_scalar_multiply_integrand

   ! -
   pure function integrand_sub_integrand(lhs, rhs) result(opr)
   !< `-` operator.
   class(integrand_ladvection), intent(in) :: lhs    !< Left hand side.
   class(integrand_object),     intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable                  :: opr(:) !< Operator result.

   select type(rhs)
   class is(integrand_ladvection)
     opr = lhs%U - rhs%U
   endselect
   endfunction integrand_sub_integrand

   pure function integrand_sub_real(lhs, rhs) result(opr)
   !< `- real` operator.
   class(integrand_ladvection), intent(in) :: lhs     !< Left hand side.
   real(R_P),                   intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable                  :: opr(:)  !< Operator result.

   opr = lhs%U - rhs
   endfunction integrand_sub_real

   pure function real_sub_integrand(lhs, rhs) result(opr)
   !< `real -` operator.
   real(R_P),                   intent(in) :: lhs(1:) !< Left hand side.
   class(integrand_ladvection), intent(in) :: rhs     !< Left hand side.
   real(R_P), allocatable                  :: opr(:)  !< Operator result.

   opr = lhs - rhs%U
   endfunction real_sub_integrand

   pure subroutine assign_integrand(lhs, rhs)
   !< `=` operator.
   class(integrand_ladvection), intent(inout) :: lhs !< Left hand side.
   class(integrand),            intent(in)    :: rhs !< Right hand side.
   integer(I_P)                               :: i   !< Counter.

   select type(rhs)
   class is(integrand_ladvection)
      lhs%weno_order = rhs%weno_order
      lhs%Ni         = rhs%Ni
      lhs%Ng         = rhs%Ng
      lhs%Dx         = rhs%Dx
      lhs%a          = rhs%a
      if (allocated(rhs%u)) then
         if (allocated(lhs%u)) deallocate(lhs%u) ; allocate(lhs%u(1:lhs%Ni))
         select type(rhs)
         class is(integrand_ladvection)
            if (allocated(rhs%U)) then
               do i=1, lhs%Ni
                  lhs%u(i) = rhs%u(i)
               enddo
            endif
         endselect
      endif
      if (allocated(rhs%BC_L)) lhs%BC_L = rhs%BC_L
      if (allocated(rhs%BC_R)) lhs%BC_R = rhs%BC_R
      if (allocated(rhs%interpolator)) then
         if (allocated(lhs%interpolator)) deallocate(lhs%interpolator)
         allocate(lhs%interpolator, mold=rhs%interpolator)
         lhs%interpolator = rhs%interpolator
      endif
   endselect
   endsubroutine assign_advection

   pure subroutine assign_real(lhs, rhs)
   !< Assign one real to an advection field.
   class(integrand_ladvection), intent(inout) :: lhs     !< Left hand side.
   real(R_P),                   intent(in)    :: rhs(1:) !< Right hand side.

   lhs%u = rhs
   endsubroutine assign_real

   ! private methods
   pure subroutine impose_boundary_conditions(self, u)
   !< Impose boundary conditions.
   class(integrand_ladvection), intent(in)    :: self          !< Advection field.
   real(R_P),                   intent(inout) :: u(1-self%Ng:) !< Conservative variables.
   integer(I_P)                               :: i             !< Space counter.

   select case(trim(adjustl(self%BC_L)))
      case('TRA') ! trasmissive (non reflective) BC
         do i=1-self%Ng, 0
            u(i) = u(-i+1)
         enddo
      case('REF') ! reflective BC
         do i=1-self%Ng, 0
            u(i) = - u(-i+1)
         enddo
      case('PER') ! periodic BC
         do i=1-self%Ng, 0
            u(i) = u(self%Ni+i)
         enddo
   endselect

   select case(trim(adjustl(self%BC_R)))
      case('TRA') ! trasmissive (non reflective) BC
         do i=self%Ni+1, self%Ni+self%Ng
            u(i) = u(self%Ni-(i-self%Ni-1))
         enddo
      case('REF') ! reflective BC
         do i=self%Ni+1, self%Ni+self%Ng
            u(i) = - u(self%Ni-(i-self%Ni-1))
         enddo
      case('PER') ! periodic BC
         do i=self%Ni+1, self%Ni+self%Ng
            u(i) = u(i-self%Ni)
         enddo
   endselect
   endsubroutine impose_boundary_conditions

   subroutine reconstruct_interfaces(self, conservative, r_conservative)
   !< Reconstruct interfaces states.
   class(integrand_ladvection), intent(in)    :: self                         !< Advection field.
   real(R_P),                   intent(in)    :: conservative(1-self%Ng:)     !< Conservative variables.
   real(R_P),                   intent(inout) :: r_conservative(1:, 0:)       !< Reconstructed conservative vars.
   real(R_P), allocatable                     :: U(:)                         !< Serialized conservative variables.
   real(R_P)                                  :: C(1:2, 1-self%Ng:-1+self%Ng) !< Stencils.
   real(R_P)                                  :: CR(1:2)                      !< Reconstrcuted intrafaces.
   integer(I_P)                               :: i                            !< Counter.
   integer(I_P)                               :: j                            !< Counter.
   integer(I_P)                               :: f                            !< Counter.

   select case(self%weno_order)
   case(1) ! 1st order piecewise constant reconstruction
      do i=0, self%Ni+1
         r_conservative(1, i) = conservative(i)
         r_conservative(2, i) = r_conservative(1, i)
      enddo
   case(3, 5, 7, 9, 11, 13, 15, 17) ! 3rd-17th order WENO reconstruction
      do i=0, self%Ni+1
         do j=i+1-self%Ng, i-1+self%Ng
            do f=1, 2
               C(f, j-i) = conservative(j)
            enddo
         enddo
         call self%interpolator%interpolate(stencil=C(:, :), interpolation=CR(:))
         do f=1, 2
            r_conservative(f, i)  = CR(f)
         enddo
      enddo
   endselect
   endsubroutine reconstruct_interfaces
endmodule foodie_test_integrand_ladvection

   subroutine save_results(results, output, scheme, frequency, U0, save_exact_solution, solution)
   !< Save results (and plots).
   logical,      intent(in)      :: results             !< Flag for activating results saving.
   character(*), intent(in)      :: output              !< Output files basename coming from CLI.
   character(*), intent(in)      :: scheme              !< Selected scheme: must be defined into *solvers*.
   real(R_P),    intent(in)      :: frequency           !< Oscillation frequency.
   real(R_P),    intent(in)      :: U0(1:)              !< Initial state.
   logical,      intent(in)      :: save_exact_solution !< Flag for saving exact solution.
   real(R_P),    intent(in)      :: solution(0:,0:)     !< Solution at each time step.
   character(len=:), allocatable :: title               !< Output files title.
   character(len=:), allocatable :: basename            !< Output files basename.
   integer(I_P)                  :: rawfile             !< Raw file unit for saving results.
   type(integrand_ladvection)    :: oscillator          !< Oscillation field.
   integer(I_P)                  :: s                   !< Counter.

   basename = trim(adjustl(output))//'-'//trim(strz(ubound(solution, dim=2), 10))//'-time_steps-'//trim(adjustl(scheme))
   title = 'oscillation equations integration, solver='//trim(adjustl(scheme))
   if (results) then
     open(newunit=rawfile, file=basename//'.dat')
     write(rawfile, '(A)')'TITLE="'//title//'"'
     write(rawfile, '(A)')'VARIABLES="t" "x" "y" "amplitude" "phase"'
     write(rawfile, '(A)')'ZONE T="'//trim(adjustl(scheme))//'"'
     do s=0, ubound(solution, dim=2)
       write(rawfile, '(5('//FR_P//',1X))')solution(:, s), amplitude_phase(solution(1:2, s))
     enddo
     close(rawfile)
   endif
   if (save_exact_solution) then
      call oscillator%initialize(U0=U0, frequency=frequency)
      basename = trim(adjustl(output))//'-'//trim(strz(ubound(solution, dim=2), 10))//'-time_steps-exact_solution'
      title = 'linear constant coefficients equation integration, solver=exact solution'
      open(newunit=rawfile, file=basename//'.dat')
      write(rawfile, '(A)')'TITLE="'//title//'"'
      write(rawfile, '(A)')'VARIABLES="t" "u"'
      write(rawfile, '(A)')'ZONE T="exact solution"'
      do s=0, ubound(solution, dim=2)
         write(rawfile, '(5('//FR_P//',1X))')solution(0, s), oscillator%exact_solution(t=solution(0, s)), &
                                             amplitude_phase(oscillator%exact_solution(t=solution(0, s)))
      enddo
      close(rawfile)
   endif
   contains
      function amplitude_phase(sol) result(ap)
      !< Compute amplitude and phase of the solution provided in X-Y domain.
      real(R_P), intent(in) :: sol(1:) !< Solution in X-Y domain.
      real(R_P)             :: ap(1:2) !< Amplitude and phase solution.

      ap(1) = sqrt(sol(1)**2 + sol(2)**2)
      ap(2) = atan(-sol(1) / sol(2))
      endfunction amplitude_phase
   endsubroutine save_results

   subroutine square_wave_initial_state(initial_state, BC_L, BC_R, Dx)
   real(R_P),    intent(inout) :: initial_state(1:)               !< Initial state of primitive variables.
   character(*), intent(inout) :: BC_L                            !< Left boundary condition.
   character(*), intent(inout) :: BC_R                            !< Rigth boundary condition.
   real(R_P),    intent(in)    :: Dx                              !< Space step.
   real(R_P)                   :: x(1:size(initial_state, dim=1)) !< Cell center x-abscissa values.
   integer(I_P)                :: i                               !< Space counter.

   BC_L = 'PER'
   BC_R = 'PER'
   do i=1, size(x, dim=1)
      x(i) = Dx * i - 0.5_R_P * Dx
      if     (x(i) < 0.25_R_P) then
         initial_state(i) = 0._R_P
      elseif (0.25_R_P <= x(i) .and. x(i) < 0.75_R_P) then
         initial_state(i) = 1._R_P
      else
         initial_state(i) = 0._R_P
      endif
   enddo
   endsubroutine square_wave_initial_state

program wenoof_test_linear_advection
!< Define [[integrand_ladvection]], the 1D linear advection PDE test field that is a concrete extension of the
!< abstract integrand type.

use flap, only : command_line_interface
use foodie, only : integrator_runge_kutta_lssp, integrator_runge_kutta_ssp
use foodie_test_integrand_ladvection, only : integrand_ladvection
use penf, only : cton, FR_P, I_P, R_P, str

implicit none
character(len=99)                          :: s_scheme                   !< Space scheme.
character(len=99)                          :: t_scheme                   !< Time scheme.
integer(I_P)                               :: weno_order                 !< WENO reconstruction order.
real(R_P)                                  :: weno_eps                   !< WENO epsilon value.
integer(I_P)                               :: stages                     !< Number of stages.
type(integrator_runge_kutta_lssp)          :: rk_integrator              !< Runge-Kutta integrator.
type(integrator_runge_kutta_ssp)           :: rk_integrator_s            !< Runge-Kutta integrator.
type(integrand_ladvection), allocatable :: rk_stage(:)                !< Runge-Kutta stages.
real(R_P)                                  :: dt                         !< Time step.
real(R_P)                                  :: t                          !< Time.
integer(I_P)                               :: step                       !< Time steps counter.
type(integrand_ladvection)              :: domain                     !< Domain of Advection equations.
real(R_P)                                  :: CFL                        !< CFL value.
character(3)                               :: BC_L                       !< Left boundary condition type.
character(3)                               :: BC_R                       !< Right boundary condition type.
integer(I_P)                               :: Ni                         !< Number of grid cells.
real(R_P)                                  :: Dx                         !< Space step discretization.
real(R_P)                                  :: a                          !< Advection coefficient.
real(R_P), allocatable                     :: x(:)                       !< Cell center x-abscissa values.
integer(I_P)                               :: steps_max                  !< Maximum number of time steps.
real(R_P)                                  :: t_max                      !< Maximum integration time.
logical                                    :: results                    !< Flag for activating results saving.
logical                                    :: time_serie                 !< Flag for activating time serie-results saving.
logical                                    :: verbose                    !< Flag for activating more verbose output.
real(R_P), parameter                       :: pi = 4._R_P * atan(1._R_P) !< Pi greek.

call parse_command_line_interface
if (verbose) print "(A)", 'Solve 1D linear advection equation " u_t + (a*u)_x =0" with a='//trim(str(a))
call initialize
call save_time_serie(filename='linear_advection-'//                                                      &
                              trim(adjustl(s_scheme))//'-'//trim(str(weno_order, no_sign=.true.))//'-'// &
                              trim(adjustl(t_scheme))//'-'//trim(str(stages, no_sign=.true.))//'-'//     &
                              'Ni_'//trim(str(Ni, no_sign=.true.))//'.dat', t=t)
step = 0
time_loop: do
   step = step + 1
   dt = domain%dt(steps_max=steps_max, t_max=t_max, t=t, CFL=CFL)
   if (trim(adjustl(t_scheme))=='runge_kutta_ssp_stages_5_order_4') then
      call rk_integrator_s%integrate(U=domain, stage=rk_stage, dt=dt, t=t)
   else
      call rk_integrator%integrate(U=domain, stage=rk_stage, dt=dt, t=t)
   endif
   t = t + dt
   call save_time_serie(t=t)
   if (verbose) print "(A)", 'step = '//str(n=step)//', time step = '//str(n=dt)//', time = '//str(n=t)
   if ((t == t_max).or.(step == steps_max)) exit time_loop
enddo time_loop
call save_time_serie(t=t, finish=.true.)

contains
   subroutine initialize()
   !< Initialize the test.
   real(R_P),   allocatable :: initial_state(:) !< Initial state of primitive variables.
   integer(I_P)             :: i                !< Space counter.

   if (trim(adjustl(t_scheme))=='runge_kutta_ssp_stages_5_order_4') then
      call rk_integrator_s%initialize(scheme=t_scheme, stop_on_fail=.true.)
      if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:rk_integrator_s%stages))
   else
      call rk_integrator%initialize(scheme=t_scheme, stages=stages)
      if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:rk_integrator%stages))
   endif
   t = 0._R_P
   if (allocated(x)) deallocate(x) ; allocate(x(1:Ni))
   if (allocated(initial_state)) deallocate(initial_state) ; allocate(initial_state(1:Ni))
   Dx = 1._R_P / Ni
   call square_wave_initial_state(initial_state=initial_state)
   call domain%initialize(Ni=Ni, Dx=Dx,                &
                          BC_L=BC_L, BC_R=BC_R,        &
                          initial_state=initial_state, &
                          s_scheme=s_scheme,           &
                          a=a,                         &
                          weno_order=weno_order)
   endsubroutine initialize

   subroutine parse_command_line_interface()
   !< Parse Command Line Interface (CLI).
   type(command_line_interface)  :: cli    !< Command line interface handler.
   integer(I_P)                  :: error  !< Error handler.
   character(len=:), allocatable :: buffer !< String buffer.

   call cli%init(description = 'WenOOF test: 1D linear advection equation', &
                 examples    = ["wenoof_test_linear_advection         ",    &
                                "wenoof_test_linear_advection --tserie"])
   call cli%add(switch='-a', help='advection coefficient', required=.false., act='store', def='1.0')
   call cli%add(switch='--Ni', help='number finite volumes used', required=.false., act='store', def='100')
   call cli%add(switch='--steps', help='number time steps performed', required=.false., act='store', def='100')
   call cli%add(switch='--t-max', help='maximum integration time', required=.false., act='store', def='0.')
   call cli%add(switch='--s-scheme', help='space intergation scheme', required=.false., act='store', def='reconstructor-JS', &
     choices='reconstructor-JS,reconstructor-M-JS,reconstructor-M-Z,reconstructor-Z')
   call cli%add(switch='--weno-order', help='WENO order', required=.false., act='store', def='1')
   call cli%add(switch='--weno-eps', help='WENO epsilon parameter', required=.false., act='store', def='0.000001')
   call cli%add(switch='--t-scheme', help='time intergation scheme', required=.false., act='store', &
                def='runge_kutta_lssp_stages_s_order_s_1',                                          &
                choices='runge_kutta_lssp_stages_s_order_s_1,runge_kutta_lssp_stages_s_order_s,runge_kutta_ssp_stages_5_order_4')
   call cli%add(switch='--stages', help='number stages', required=.false., act='store', def='2')
   call cli%add(switch='--cfl', help='CFL value', required=.false., act='store', def='0.8')
   call cli%add(switch='--tserie', switch_ab='-t', help='Save time-serie-result', required=.false., act='store_true', def='.false.')
   call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.')
   call cli%parse(error=error)
   call cli%get(switch='-a',           val=a,          error=error) ; if (error/=0) stop
   call cli%get(switch='--Ni',         val=Ni,         error=error) ; if (error/=0) stop
   call cli%get(switch='--steps',      val=steps_max,  error=error) ; if (error/=0) stop
   call cli%get(switch='--t-max',      val=t_max,      error=error) ; if (error/=0) stop
   call cli%get(switch='--s-scheme',   val=s_scheme,   error=error) ; if (error/=0) stop
   call cli%get(switch='--weno-order', val=weno_order, error=error) ; if (error/=0) stop
   call cli%get(switch='--weno-eps',   val=weno_eps,   error=error) ; if (error/=0) stop
   call cli%get(switch='--t-scheme',   val=t_scheme,   error=error) ; if (error/=0) stop
   call cli%get(switch='--stages',     val=stages,     error=error) ; if (error/=0) stop
   call cli%get(switch='--cfl',        val=CFL,        error=error) ; if (error/=0) stop
   call cli%get(switch='--tserie',     val=time_serie, error=error) ; if (error/=0) stop
   call cli%get(switch='--verbose',    val=verbose,    error=error) ; if (error/=0) stop

   if (t_max > 0._R_P) steps_max = 0
   endsubroutine parse_command_line_interface

   subroutine square_wave_initial_state(initial_state)
   real(R_P),   intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I_P)               :: i                 !< Space counter.

   BC_L = 'PER'
   BC_R = 'PER'
   do i=1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      if     (x(i) < 0.25_R_P) then
         initial_state(i) = 0._R_P
      elseif (0.25_R_P <= x(i) .and. x(i) < 0.75_R_P) then
         initial_state(i) = 1._R_P
      else
         initial_state(i) = 0._R_P
      endif
   enddo
   endsubroutine square_wave_initial_state

   subroutine save_time_serie(t, filename, finish)
   !< Save time-serie results.
   real(R_P),    intent(in)           :: t         !< Current integration time.
   character(*), intent(in), optional :: filename  !< Output filename.
   logical,      intent(in), optional :: finish    !< Flag for triggering the file closing.
   integer(I_P), save                 :: tsfile    !< File unit for saving time serie results.
   integer(I_P)                       :: i         !< Counter.

   if (time_serie) then
      if (present(filename)) then
         open(newunit=tsfile, file=filename)
      endif
      write(tsfile, '(A)')'VARIABLES = "x" "u"'
      write(tsfile, '(A)')'ZONE T="'//str(n=t)//'"'
      do i=1, Ni
         write(tsfile, '(4'//'('//FR_P//',1X))')x(i), domain%u(i)
      enddo
      if (present(finish)) then
         if (finish) close(tsfile)
      endif
   endif
   endsubroutine save_time_serie
endprogram wenoof_test_linear_advection
