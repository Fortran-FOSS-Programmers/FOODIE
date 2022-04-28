!< Define [[integrand_euler_1D]], the 1D Euler PDE system field that is a concrete extension of the abstract integrand type.

module foodie_integrand_euler_1D
!< Define [[integrand_euler_1D]], the 1D Euler PDE system field that is a concrete extension of the abstract integrand type.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use flap, only : command_line_interface
use foodie, only : integrand_object, integrator_object, integrator_multistage_object, foodie_integrator_factory
use penf, only : FR_P, I_P, R_P, str, strz
use wenoof, only : interpolator_object, wenoof_create

implicit none
private
public :: integrand_euler_1D

type, extends(integrand_object) :: integrand_euler_1D
   !< 1D Euler PDE system field.
   !<
   !< It is a FOODIE integrand class concrete extension.
   !<
   !<### 1D Euler PDEs system
   !<
   !< The 1D Euler PDEs system considered is a non linear, hyperbolic (inviscid) system of conservation laws for compressible gas
   !< dynamics, that reads as
   !<$$
   !<\begin{matrix}
   !<U_t = R(U)  \Leftrightarrow U_t = F(U)_x \\
   !<U = \begin{bmatrix}
   !<\rho \\
   !<\rho u \\
   !<\rho E
   !<\end{bmatrix}\;\;\;
   !<F(U) = \begin{bmatrix}
   !<\rho u \\
   !<\rho u^2 + p \\
   !<\rho u H
   !<\end{bmatrix}
   !<\end{matrix}
   !<$$
   !< where \(\rho\) is the density, \(u\) is the velocity, \(p\) the pressure, \(E\) the total internal specific energy and \(H\)
   !< the total specific enthalpy. The PDEs system must completed with the proper initial and boundary conditions. Moreover, an ideal
   !< (thermally and calorically perfect) gas is considered
   !<$$
   !<\begin{matrix}
   !<R = c_p - c_v \\
   !<\gamma = \frac{c_p}{c_v}\\
   !<e = c_v T \\
   !<h = c_p T
   !<\end{matrix}
   !<$$
   !< where *R* is the gas constant, \(c_p\,c_v\) are the specific heats at constant pressure and volume (respectively), *e* is the
   !< internal energy, *h* is the internal enthalpy and *T* is the temperature. The following addition equations of state hold:
   !<$$
   !<\begin{matrix}
   !<T = \frac{p}{\rho R} \\
   !<E = \rho e + \frac{1}{2} \rho u^2 \\
   !<H = \rho h + \frac{1}{2} \rho u^2 \\
   !<a = \sqrt{\frac{\gamma p}{\rho}}
   !<\end{matrix}
   !<$$
   !<
   !<### Multi-fluid Euler PDEs system
   !< An extension of the above Euler system is considered allowing the modelling of a multi-fluid mixture of different gas (with
   !< different physical characteristics). The well known Standard Thermodynamic Model is used to model the gas mixture replacing the
   !< density with the density fraction of each specie composing the mixture. This led to the following system:
   !<$$
   !<\begin{matrix}
   !<U_t = R(U)  \Leftrightarrow U_t = F(U)_x \\
   !<U = \begin{bmatrix}
   !<\rho_s \\
   !<\rho u \\
   !<\rho E
   !<\end{bmatrix}\;\;\;
   !<F(U) = \begin{bmatrix}
   !<\rho_s u \\
   !<\rho u^2 + p \\
   !<\rho u H
   !<\end{bmatrix}\;\;\; for\; s=1,2,...N_s \\
   !<\rho = \sum_{s=1}^{N_s}\rho_s \\
   !<c_p = \sum_{s=1}^{N_S} \frac{\rho_s}{\rho} c_{p,s} \quad  c_v = \sum_{s=1}^{N_S} \frac{\rho_s}{\rho} c_{v,s}
   !<\end{matrix}
   !<$$
   !< where \(N_s\) is the number of initial species composing the gas mixture.
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
   !<
   !<#### Primitive variables organization
   !< Primitive variables are organized as an array of reals which the first index means:
   !<
   !< + 1    : density of species 1    (r1)
   !< + 2    : density of species 2    (r2)
   !< + ...  :
   !< + s    : density of species s-th (rs)
   !< + ...  :
   !< + Ns   : density of species Ns   (rNs)
   !< + Ns+1 : velocity                (u)
   !< + Ns+2 : pressure                (p)
   !< + Ns+3 : density                 (r=sum(rs))
   !< + Ns+4 : specific heats ratio    (g)
   !<
   !<#### Conservative variables organization
   !< Conservative variables are organized as an array (rank 2) of reals which the first index means:
   !<
   !< + 1    : mass conservation of species 1    (r1)
   !< + 2    : mass conservation of species 2    (r2)
   !< + ...  :
   !< + s    : mass conservation of species s-th (rs)
   !< + ...  :
   !< + Ns   : mass conservation of species Ns   (rNs)
   !< + Ns+1 : momentum conservation             (r*u)
   !< + Ns+2 : energy conservation               (r*E)
   !<
   private
   type(command_line_interface)            :: cli              !< Command line interface handler.
   character(99)                           :: w_scheme=''      !< WENO Scheme used.
   integer(I_P)                            :: weno_order=0     !< WENO reconstruction order.
   real(R_P)                               :: weno_eps=0._R_P  !< WENO epsilon to avoid division by zero, default value.
   class(interpolator_object), allocatable :: interpolator     !< WENO interpolator.
   real(R_P)                               :: CFL=0._R_P       !< CFL value.
   integer(I_P)                            :: Ni=0             !< Space dimension.
   integer(I_P)                            :: Ng=0             !< Ghost cells number.
   integer(I_P)                            :: Ns=0             !< Number of initial species.
   integer(I_P)                            :: Nc=0             !< Number of conservative variables, Ns+2.
   integer(I_P)                            :: Np=0             !< Number of primitive variables, Ns+4.
   real(R_P)                               :: length=0._R_P    !< Domain length.
   real(R_P)                               :: dx=0._R_P        !< Space step.
   real(R_P)                               :: Tmax=1._R_P      !< Maximum time.
   integer(I_P)                            :: Nmax=-1          !< Maximum time steps.
   real(R_P), allocatable                  :: U(:,:)           !< Integrand (state) variables, physical domain [1:Nc,1:Ni].
   real(R_P), allocatable                  :: cp0(:)           !< Specific heat cp of initial species [1:Ns].
   real(R_P), allocatable                  :: cv0(:)           !< Specific heat cv of initial species [1:Ns].
   character(99)                           :: i_scheme=''      !< Integrator scheme used.
   class(integrator_object), allocatable   :: integrator       !< The integrator.
   character(99)                           :: initial_state='' !< Initial state.
   character(99)                           :: BC_L=''          !< Left boundary condition type.
   character(99)                           :: BC_R=''          !< Right boundary condition type.
   contains
      ! auxiliary methods
      procedure, pass(self), public :: destroy        !< Destroy field.
      procedure, pass(self), public :: dt             !< Return the current time step by means of CFL condition.
      procedure, pass(self), public :: export_csv     !< Export integrand to CSV file.
      procedure, pass(self), public :: export_tecplot !< Export integrand to Tecplot file.
      procedure, pass(self), public :: initialize     !< Initialize integrand.
      procedure, pass(self), public :: integrate      !< Integrate itself.
      procedure, pass(self), public :: is_completed   !< Return true if the integration is completed.
      procedure, pass(self), public :: output         !< Extract integrand state field.
      ! integrand_object deferred methods
      procedure, pass(self), public :: description         !< Return an informative description of the integrand.
      procedure, pass(self), public :: integrand_dimension !< Return integrand dimension.
      procedure, pass(self), public :: t => dEuler_dt      !< Time derivative, residuals.
      ! operators
      procedure, pass(lhs), public :: local_error !<`||integrand_euler_1D - integrand_euler_1D||` operator.
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
      ! private methods
      procedure, pass(self), private :: primitive2conservative        !< Convert primitive variables to conservative ones.
      procedure, pass(self), private :: conservative2primitive        !< Convert conservative variables to primitive ones.
      procedure, pass(self), private :: impose_boundary_conditions    !< Impose boundary conditions.
      procedure, pass(self), private :: reconstruct_interfaces_states !< Reconstruct interfaces states.
      procedure, pass(self), private :: riemann_solver                !< Solve the Riemann Problem at cell interfaces.
endtype integrand_euler_1D

contains
   ! auxiliary methods
   subroutine destroy(self)
   !< Destroy field.
   class(integrand_euler_1D), intent(inout) :: self  !< Euler field.
   type(integrand_euler_1D)                 :: fresh !< Fresh field to reset self.

   self = fresh
   endsubroutine destroy

   pure function dt(self, t)
   !< Compute the current time step by means of CFL condition.
   class(integrand_euler_1D), intent(in) :: self !< Euler field.
   real(R_P),                 intent(in) :: t    !< Time.
   real(R_P)                             :: dt   !< Time step.
   real(R_P), allocatable                :: P(:) !< Primitive variables.
   real(R_P)                             :: vmax !< Maximum propagation speed of signals.
   integer(I_P)                          :: i    !< Counter.

   associate(Ni=>self%Ni, Ns=>self%Ns, dx=>self%dx, CFL=>self%CFL, Nmax=>self%Nmax, Tmax=>self%Tmax)
      vmax = 0._R_P
      do i=1, Ni
         P    = self%conservative2primitive(self%U(:, i))
         vmax = max(abs(P(Ns+1)) + a(p=P(Ns+2), r=P(Ns+3), g=P(Ns+4)), vmax)
      enddo
      dt = dx * CFL / vmax
      if (Nmax<=0) then
         if ((t + dt) > Tmax) dt = Tmax - t
      endif
   endassociate
   endfunction dt

   subroutine export_csv(self, file_name, t)
   !< Export integrand to CSV file.
   class(integrand_euler_1D), intent(in) :: self      !< Advection field.
   character(*),              intent(in) :: file_name !< File name.
   real(R_P),                 intent(in) :: t         !< Time.
   integer(I_P), save                    :: file_unit !< File unit.
   real(R_P), allocatable                :: P(:)      !< Primitive variables.
   integer(I_P)                          :: i         !< Counter.

   open(newunit=file_unit, file=trim(adjustl(file_name)))
   write(unit=file_unit, fmt='(A)') 'x,u,p,r,t'
   do i=1, self%Ni
      P = self%conservative2primitive(self%U(:, i))
      write(unit=file_unit, fmt='(A)') trim(str(self%dx*(i-0.5_R_P)))//','// &
                                       trim(str(P(self%Ns+1)))//','//        &
                                       trim(str(P(self%Ns+2)))//','//        &
                                       trim(str(P(self%Ns+3)))//','//        &
                                       trim(str(t))
   enddo
   close(unit=file_unit)
   endsubroutine export_csv

   subroutine export_tecplot(self, file_name, t, close_file)
   !< Export integrand to Tecplot file.
   class(integrand_euler_1D), intent(in)           :: self            !< Advection field.
   character(*),              intent(in), optional :: file_name       !< File name.
   real(R_P),                 intent(in), optional :: t               !< Time.
   logical,                   intent(in), optional :: close_file      !< Flag for closing file.
   logical, save                                   :: is_open=.false. !< Flag for checking if file is open.
   integer(I_P), save                              :: file_unit       !< File unit.
   real(R_P), allocatable                          :: P(:)            !< Primitive variables.
   integer(I_P)                                    :: i               !< Counter.

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
         write(unit=file_unit, fmt='(A)') 'VARIABLES="x" "u" "p" "r"'
      endif
      if (present(t) .and. is_open) then
         write(unit=file_unit, fmt='(A)') 'ZONE T="'//str(t)//&
                                          '" I='//trim(str(self%Ni))//', J=1, K=1, ZONETYPE = Ordered, DATAPACKING = POINT'
         do i=1, self%Ni
            P = self%conservative2primitive(self%U(:, i))
            write(unit=file_unit, fmt='(A)') trim(str(self%dx*(i-0.5_R_P)))//' '// &
                                             trim(str(P(self%Ns+1)))//' '//        &
                                             trim(str(P(self%Ns+2)))//' '//        &
                                             trim(str(P(self%Ns+3)))
         enddo
      endif
   endif
   endsubroutine export_tecplot

   subroutine initialize(self)
   !< Initialize integrand.
   class(integrand_euler_1D), intent(inout) :: self !< Euler field.
   integer(I_P)                             :: i    !< Space counter.

   call self%destroy

   call set_cli
   call parse_cli

   self%Nc = self%Ns + 2
   self%Np = self%Ns + 4
   self%Ng = (self%weno_order + 1) / 2
   self%dx = self%length / self%Ni
   if (allocated(self%U)) deallocate(self%U) ; allocate(self%U(1:self%Nc, 1:self%Ni))
   if (allocated(self%cp0)) deallocate(self%cp0) ; allocate(self%cp0(1:self%Ns))
   if (allocated(self%cv0)) deallocate(self%cv0) ; allocate(self%cv0(1:self%Ns))
   if (self%weno_order>1) call wenoof_create(interpolator_type=trim(adjustl(self%w_scheme)), &
                                             S=self%Ng,                                      &
                                             interpolator=self%interpolator,                 &
                                             eps=self%weno_eps)
   call foodie_integrator_factory(scheme=self%i_scheme, integrator=self%integrator, U=self)

   select case(trim(adjustl(self%initial_state)))
   case('sod')
      call set_initial_state_sod
   case('bryson')
      call set_initial_state_bryson
   endselect
   contains
      subroutine parse_cli
      !< Initialize from command line interface.

      call self%cli%get(switch='--w-scheme',      val=self%w_scheme,      error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--weno-order',    val=self%weno_order,    error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--weno-eps',      val=self%weno_eps,      error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--cfl',           val=self%CFL,           error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--length',        val=self%length,        error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--Ni',            val=self%Ni,            error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--Ns',            val=self%Ns,            error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--initial_state', val=self%initial_state, error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--BC_L',          val=self%BC_L,          error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--BC_R',          val=self%BC_R,          error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--i-scheme',      val=self%i_scheme,      error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--Tmax',          val=self%Tmax,          error=self%cli%error) ; if (self%cli%error/=0) stop
      call self%cli%get(switch='--Nmax',          val=self%Nmax,          error=self%cli%error) ; if (self%cli%error/=0) stop
      endsubroutine parse_cli

      subroutine set_cli
      !< Set command line interface.

      call self%cli%init(progname    = 'euler_1D',                                                     &
                         authors     = 'Fortran-FOSS-Programmers',                                     &
                         license     = 'GNU GPLv3',                                                    &
                         description = 'FOODIE library APP, 1D Euler equations integration',           &
                         examples    = ["foodie_euler_1D --help                                     ", &
                                        "foodie_euler_1D                                            ", &
                                        "foodie_euler_1D --initial_state bryson --length 20 --Tmax 6"])
      call self%cli%add(switch='--w-scheme', help='WENO scheme', required=.false., act='store', def='reconstructor-JS', &
                        choices='reconstructor-JS,reconstructor-M-JS,reconstructor-M-Z,reconstructor-Z')
      call self%cli%add(switch='--weno-order', help='WENO order', required=.false., act='store', def='3')
      call self%cli%add(switch='--weno-eps', help='WENO epsilon parameter', required=.false., act='store', def='0.000001')
      call self%cli%add(switch='--cfl', help='CFL value', required=.false., act='store', def='0.8')
      call self%cli%add(switch='--length', help='domain lenth', required=.false., act='store', def='1.0')
      call self%cli%add(switch='--Ni', help='number finite volumes used', required=.false., act='store', def='100')
      call self%cli%add(switch='--Ns', help='number initial species', required=.false., act='store', def='1')
      call self%cli%add(switch='--initial_state', help='initial state', required=.false., act='store', def='sod', &
                        choices='sod,bryson')
      call self%cli%add(switch='--BC_L',help='boundary condition, left', required=.false.,act='store',def='TRA',choices='TRA,REF')
      call self%cli%add(switch='--BC_R',help='boundary condition, right',required=.false.,act='store',def='TRA',choices='TRA,REF')
      call self%cli%add(switch='--i-scheme', help='integrator scheme', required=.false., act='store', &
                        def='runge_kutta_ssp_stages_5_order_4',                                       &
                        choices='runge_kutta_ssp_stages_1_order_1,'//&
                                'runge_kutta_ssp_stages_3_order_3,'//&
                                'runge_kutta_ssp_stages_5_order_4')
      call self%cli%add(switch='--Tmax', help='maximum integration time', required=.false., act='store', def='0.2')
      call self%cli%add(switch='--Nmax', help='maximum number of time steps', required=.false., act='store', def='-1')
      endsubroutine set_cli

      subroutine set_initial_state_sod
      !< Set SOD initial state.

      self%cp0(1) = 1040._R_P
      self%cv0(1) = 743._R_P
      do i=1, self%Ni/2
         self%U(:, i) = self%primitive2conservative([1._R_P, &                 ! rho(s)
                                                     0._R_P, &                 ! u
                                                     1._R_P, &                 ! p
                                                     1._R_P, &                 ! sum(rho(s))
                                                     self%cp0(1)/self%cv0(1)]) ! gamma = cp/cv)
      enddo
      do i=self%Ni/2 + 1, self%Ni
         self%U(:, i) = self%primitive2conservative([0.125_R_P, &              ! rho(s)
                                                     0._R_P,    &              ! u
                                                     0.1_R_P,   &              ! p
                                                     0.125_R_P, &              ! sum(rho(s))
                                                     self%cp0(1)/self%cv0(1)]) ! gamma = cp/cv)
      enddo
      endsubroutine set_initial_state_sod

      subroutine set_initial_state_bryson
      !< Set SOD initial state.

      self%cp0(1) = 1040._R_P
      self%cv0(1) = 743._R_P
      do i=1, self%Ni
         if (i*self%dx<=1._R_P) then
            self%U(:, i) = self%primitive2conservative([5.1432014_R_P,          & ! rho(s)
                                                        2.0451067615658363_R_P, & ! u
                                                        9.04545_R_P,            & ! p
                                                        5.1432014_R_P,          & ! sum(rho(s))
                                                        self%cp0(1)/self%cv0(1)]) ! gamma = cp/cv)
         else
            self%U(:, i) = self%primitive2conservative([1.4_R_P,   &              ! rho(s)
                                                        0._R_P,    &              ! u
                                                        1._R_P,    &              ! p
                                                        1.4_R_P,   &              ! sum(rho(s))
                                                        self%cp0(1)/self%cv0(1)]) ! gamma = cp/cv)
         endif
      enddo
      endsubroutine set_initial_state_bryson
   endsubroutine initialize

   subroutine integrate(self)
   !< Integrate itself.
   class(integrand_euler_1D), intent(inout) :: self !< Euler field.
   real(R_P)                                :: dt   !< Time step.
   real(R_P)                                :: t    !< Current time.
   integer(I_P)                             :: n    !< Current time step.

   associate(integrator=>self%integrator)
      call self%initialize
      print '(A)', 'Integrate '//trim(adjustl(self%description()))
      print '(A)', 'Nmax              : '//trim(str(self%Nmax))
      print '(A)', 'Tmax              : '//trim(str(self%Tmax))
      print '(A)', 'CFL               : '//trim(str(self%CFL))
      print '(A)', 'WENO scheme       : '//trim(self%w_scheme)
      print '(A)', 'Integrator scheme : '//trim(self%i_scheme)
      print '(A)', 'Domain length     : '//trim(str(self%length))
      print '(A)', 'Number of FV (Ni) : '//trim(str(self%Ni))
      print '(A)', 'Number of IS (Ns) : '//trim(str(self%Ns))
      print '(A)', 'BC left           : '//trim(self%BC_L)
      print '(A)', 'BC right          : '//trim(self%BC_R)
      n = 0
      t = 0._R_P
      call self%export_csv(file_name=trim(adjustl(self%description()))//'-'//trim(strz(n,9))//'.csv', t=t)
      do
         n = n + 1
         dt = self%dt(t=t)
         select type(integrator)
         class is(integrator_multistage_object)
            call integrator%integrate(U=self, dt=dt, t=t)
         endselect
         t = t + dt
         print '(A)', 'n, t, dt: '//trim(str(n))//', '//trim(str(t))//', '//trim(str(dt))
         call self%export_csv(file_name=trim(adjustl(self%description()))//'-'//trim(strz(n,9))//'.csv', t=t)
         if (self%is_completed(n=n, t=t)) exit
      enddo
   endassociate
   endsubroutine integrate

   pure function is_completed(self, n, t)
   !< Return true if the integration is completed.
   class(integrand_euler_1D), intent(in) :: self         !< Euler field.
   integer(I_P),              intent(in) :: n            !< Time step.
   real(R_P),                 intent(in) :: t            !< Time.
   logical                               :: is_completed !< Returned status.

   is_completed = .false.
   associate(Nmax=>self%Nmax, Tmax=>self%Tmax)
      if ((Nmax<=0)) then
         if ((t>=Tmax)) then
            is_completed = .true.
         endif
      elseif (n>=Nmax) then
         is_completed = .true.
      endif
   endassociate
   endfunction is_completed

   pure function output(self) result(state)
   !< Output the advection field state.
   class(integrand_euler_1D), intent(in) :: self       !< Euler field.
   real(R_P), allocatable                :: state(:,:) !< Euler field state

   state = self%U
   endfunction output

   ! integrand_tester_object deferred methods
   pure function description(self, prefix) result(desc)
   !< Return informative integrator description.
   class(integrand_euler_1D), intent(in)           :: self    !< Euler field.
   character(*),              intent(in), optional :: prefix  !< Prefixing string.
   character(len=:), allocatable                   :: desc    !< Description.
   character(len=:), allocatable                   :: prefix_ !< Prefixing string, local variable.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = prefix//'euler_1D-'//trim(self%initial_state)
   endfunction description

   function local_error(lhs, rhs) result(error)
   !< Estimate local truncation error between 2 Euler approximations.
   !<
   !< The estimation is done by norm L2 of U:
   !<
   !< $$ error = \sqrt{ \sum_i{\sum_i{ \frac{(lhs\%U_i - rhs\%U_i)^2}{lhs\%U_i^2} }} } $$
   class(integrand_euler_1D), intent(in) :: lhs   !< Left hand side.
   class(integrand_object),   intent(in) :: rhs   !< Right hand side.
   real(R_P)                             :: error !< Error estimation.
   integer(I_P)                          :: i     !< Space counter.
   integer(I_P)                          :: v     !< Variables counter.

   select type(rhs)
   class is (integrand_euler_1D)
      error = 0._R_P
      do i=1, lhs%Ni
         do v=1, lhs%Nc
            error = error + (lhs%U(v, i) - rhs%U(v, i))**2/lhs%U(v, i)**2
         enddo
      enddo
      error = sqrt(error)
   endselect
   endfunction local_error

   ! integrand_object deferred methods
   function dEuler_dt(self, t) result(dState_dt)
   !< Time derivative of Euler field, the residuals function.
   class(integrand_euler_1D), intent(in)           :: self         !< Euler field.
   real(R_P),                 intent(in), optional :: t            !< Time.
   real(R_P), allocatable                          :: dState_dt(:) !< Euler field time derivative.
   real(R_P), allocatable                          :: F(:,:)       !< Fluxes of conservative variables.
   real(R_P), allocatable                          :: P(:,:)       !< Primitive variables.
   real(R_P), allocatable                          :: PR(:,:,:)    !< Left (1)/right (2) (rec.) interface values of prim. var.
   integer(I_P)                                    :: c, i, j      !< Counter.

   allocate(F(1:self%Nc, 0:self%Ni))
   F = 0._R_P
   allocate(P(1:self%Np, 1-self%Ng:self%Ni+self%Ng))
   P = 0._R_P
   allocate(PR(1:self%Np, 1:2, 0:self%Ni+1))
   PR = 0._R_P
   ! compute primitive variables
   do i=1, self%Ni
      P(:, i) = self%conservative2primitive(self%U(:, i))
   enddo
   call self%impose_boundary_conditions(primitive=P)
   call self%reconstruct_interfaces_states(primitive=P, r_primitive=PR)
   ! compute fluxes by solving Rimeann Problems at each interface
   do i=0, self%Ni
      call self%riemann_solver(r1=PR(self%Ns+3, 2, i  ), &
                               u1=PR(self%Ns+1, 2, i  ), &
                               p1=PR(self%Ns+2, 2, i  ), &
                               g1=PR(self%Ns+4, 2, i  ), &
                               r4=PR(self%Ns+3, 1, i+1), &
                               u4=PR(self%Ns+1, 1, i+1), &
                               p4=PR(self%Ns+2, 1, i+1), &
                               g4=PR(self%Ns+4, 1, i+1), &
                               F=F(:, i))
      if (self%Ns>1) then
         if (F(1, i)>0._R_P) then
            F(1:self%Ns, i) = PR(1:self%Ns, 2, i  )/PR(self%Ns+3, 2, i  )*F(1, i)
         else
            F(1:self%Ns, i) = PR(1:self%Ns, 1, i+1)/PR(self%Ns+3, 1, i+1)*F(1, i)
         endif
      endif
   enddo
   ! compute residuals
   allocate(dState_dt(1:self%Nc*self%Ni))
   j = 0
   do i=1, self%Ni
      do c=1, self%Nc
         j = j + 1
         dState_dt(j) = (F(c, i - 1) - F(c, i)) / self%dx
      enddo
   enddo
   endfunction dEuler_dt

   pure function integrand_dimension(self)
   !< return integrand dimension.
   class(integrand_euler_1D), intent(in) :: self                !< Euler field.
   integer(I_P)                          :: integrand_dimension !< Integrand dimension.

   integrand_dimension = self%Nc*self%Ni
   endfunction integrand_dimension

   ! +
   pure function integrand_add_integrand(lhs, rhs) result(opr)
   !< `+` operator.
   class(integrand_euler_1D), intent(in) :: lhs     !< Left hand side.
   class(integrand_object),   intent(in) :: rhs     !< Right hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   select type(rhs)
   class is(integrand_euler_1D)
      allocate(opr(1:lhs%Nc*lhs%Ni))
      j = 0
      do i=1, lhs%Ni
         do c=1, lhs%Nc
            j = j + 1
            opr(j) = lhs%U(c,i) + rhs%U(c,i)
         enddo
      enddo
   endselect
   endfunction integrand_add_integrand

   pure function integrand_add_real(lhs, rhs) result(opr)
   !< `+ real` operator.
   class(integrand_euler_1D), intent(in) :: lhs     !< Left hand side.
   real(R_P),                 intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   allocate(opr(1:lhs%Nc*lhs%Ni))
   j = 0
   do i=1, lhs%Ni
      do c=1, lhs%Nc
         j = j + 1
         opr(j) = lhs%U(c,i) + rhs(j)
      enddo
   enddo
   endfunction integrand_add_real

   pure function real_add_integrand(lhs, rhs) result(opr)
   !< `real +` operator.
   real(R_P),                 intent(in) :: lhs(1:) !< Left hand side.
   class(integrand_euler_1D), intent(in) :: rhs     !< Left hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   allocate(opr(1:rhs%Nc*rhs%Ni))
   j = 0
   do i=1, rhs%Ni
      do c=1, rhs%Nc
         j = j + 1
         opr(j) = lhs(j) + rhs%U(c,i)
      enddo
   enddo
   endfunction real_add_integrand

   ! *
   pure function integrand_multiply_integrand(lhs, rhs) result(opr)
   !< `*` operator.
   class(integrand_euler_1D), intent(in) :: lhs     !< Left hand side.
   class(integrand_object),   intent(in) :: rhs     !< Right hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   select type(rhs)
   class is(integrand_euler_1D)
      allocate(opr(1:lhs%Nc*lhs%Ni))
      j = 0
      do i=1, lhs%Ni
         do c=1, lhs%Nc
            j = j + 1
            opr(j) = lhs%U(c,i) * rhs%U(c,i)
         enddo
      enddo
   endselect
   endfunction integrand_multiply_integrand

   pure function integrand_multiply_real(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(integrand_euler_1D), intent(in) :: lhs     !< Left hand side.
   real(R_P),                 intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   allocate(opr(1:lhs%Nc*lhs%Ni))
   j = 0
   do i=1, lhs%Ni
      do c=1, lhs%Nc
         j = j + 1
         opr(j) = lhs%U(c,i) * rhs(j)
      enddo
   enddo
   endfunction integrand_multiply_real

   pure function real_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   real(R_P),                 intent(in) :: lhs(1:) !< Left hand side.
   class(integrand_euler_1D), intent(in) :: rhs     !< Right hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   allocate(opr(1:rhs%Nc*rhs%Ni))
   j = 0
   do i=1, rhs%Ni
      do c=1, rhs%Nc
         j = j + 1
         opr(j) = lhs(j) * rhs%U(c,i)
      enddo
   enddo
   endfunction real_multiply_integrand

   pure function integrand_multiply_real_scalar(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(integrand_euler_1D), intent(in) :: lhs     !< Left hand side.
   real(R_P),                 intent(in) :: rhs     !< Right hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   allocate(opr(1:lhs%Nc*lhs%Ni))
   j = 0
   do i=1, lhs%Ni
      do c=1, lhs%Nc
         j = j + 1
         opr(j) = lhs%U(c,i) * rhs
      enddo
   enddo
   endfunction integrand_multiply_real_scalar

   pure function real_scalar_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   real(R_P),                 intent(in) :: lhs     !< Left hand side.
   class(integrand_euler_1D), intent(in) :: rhs     !< Right hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   allocate(opr(1:rhs%Nc*rhs%Ni))
   j = 0
   do i=1, rhs%Ni
      do c=1, rhs%Nc
         j = j + 1
         opr(j) =  lhs * rhs%U(c,i)
      enddo
   enddo
   endfunction real_scalar_multiply_integrand

   ! -
   pure function integrand_sub_integrand(lhs, rhs) result(opr)
   !< `-` operator.
   class(integrand_euler_1D), intent(in) :: lhs     !< Left hand side.
   class(integrand_object),   intent(in) :: rhs     !< Right hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   select type(rhs)
   class is(integrand_euler_1D)
      allocate(opr(1:lhs%Nc*lhs%Ni))
      j = 0
      do i=1, lhs%Ni
         do c=1, lhs%Nc
            j = j + 1
            opr(j) = lhs%U(c,i) - rhs%U(c,i)
         enddo
      enddo
   endselect
   endfunction integrand_sub_integrand

   pure function integrand_sub_real(lhs, rhs) result(opr)
   !< `- real` operator.
   class(integrand_euler_1D), intent(in) :: lhs     !< Left hand side.
   real(R_P),                 intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   allocate(opr(1:lhs%Nc*lhs%Ni))
   j = 0
   do i=1, lhs%Ni
      do c=1, lhs%Nc
         j = j + 1
         opr(j) = lhs%U(c,i) - rhs(j)
      enddo
   enddo
   endfunction integrand_sub_real

   pure function real_sub_integrand(lhs, rhs) result(opr)
   !< `real -` operator.
   real(R_P),                 intent(in) :: lhs(1:) !< Left hand side.
   class(integrand_euler_1D), intent(in) :: rhs     !< Left hand side.
   real(R_P), allocatable                :: opr(:)  !< Operator result.
   integer(I_P)                          :: c, i, j !< Counter.

   allocate(opr(1:rhs%Nc*rhs%Ni))
   j = 0
   do i=1, rhs%Ni
      do c=1, rhs%Nc
         j = j + 1
         opr(j) = lhs(j) - rhs%U(c,i)
      enddo
   enddo
   endfunction real_sub_integrand

   ! =
   subroutine assign_integrand(lhs, rhs)
   !< `=` operator.
   class(integrand_euler_1D), intent(inout) :: lhs !< Left hand side.
   class(integrand_object),   intent(in)    :: rhs !< Right hand side.

   select type(rhs)
   class is(integrand_euler_1D)
      lhs%cli        = rhs%cli
      lhs%w_scheme   = rhs%w_scheme
      lhs%weno_order = rhs%weno_order
      lhs%weno_eps   = rhs%weno_eps
      lhs%CFL        = rhs%CFL
      lhs%Ni         = rhs%Ni
      lhs%Ng         = rhs%Ng
      lhs%Ns         = rhs%Ns
      lhs%Nc         = rhs%Nc
      lhs%Np         = rhs%Np
      lhs%length     = rhs%length
      lhs%dx         = rhs%dx
      if (allocated(rhs%U)) then
         lhs%U = rhs%U
      else
         if (allocated(lhs%U)) deallocate(lhs%U)
      endif
      if (allocated(rhs%interpolator)) then
         if (allocated(lhs%interpolator)) then
            call lhs%interpolator%destroy
            deallocate(lhs%interpolator)
         endif
         allocate(lhs%interpolator, mold=rhs%interpolator)
         lhs%interpolator = rhs%interpolator
      else
         if (allocated(lhs%interpolator)) deallocate(lhs%interpolator)
      endif
      if (allocated(rhs%integrator)) then
         if (allocated(lhs%integrator)) deallocate(lhs%integrator)
         allocate(lhs%integrator, mold=rhs%integrator)
         lhs%integrator = rhs%integrator
      else
         if (allocated(lhs%integrator)) deallocate(lhs%integrator)
      endif
      lhs%initial_state = rhs%initial_state
      if (allocated(rhs%cp0))  lhs%cp0  = rhs%cp0
      if (allocated(rhs%cv0))  lhs%cv0  = rhs%cv0
      lhs%BC_L = rhs%BC_L
      lhs%BC_R = rhs%BC_R
   endselect
   endsubroutine assign_integrand

   pure subroutine assign_real(lhs, rhs)
   !< Assign one real to an advection field.
   class(integrand_euler_1D), intent(inout) :: lhs     !< Left hand side.
   real(R_P),                 intent(in)    :: rhs(1:) !< Right hand side.
   integer(I_P)                             :: c, i, j !< Counter.

   j = 0
   do i=1, lhs%Ni
      do c=1, lhs%Nc
         j = j + 1
         lhs%U(c,i) = rhs(j)
      enddo
   enddo
   endsubroutine assign_real

   ! private methods
   pure function primitive2conservative(self, primitive) result(conservative)
   !< Convert primitive variables to conservative variables.
   class(integrand_euler_1D), intent(in) :: self                    !< Euler field.
   real(R_P),                 intent(in) :: primitive(:)            !< Primitive variables.
   real(R_P)                             :: conservative(1:self%Nc) !< Conservative variables.

   associate(Ns=>self%Ns)
      conservative(1:Ns) = primitive(1:Ns)
      conservative(Ns+1) = primitive(Ns + 3) * primitive(Ns + 1)
      conservative(Ns+2) = primitive(Ns + 2) / (primitive(Ns + 4) - 1._R_P) + &
                           0.5_R_P*primitive(Ns + 3) * primitive(Ns + 1) * primitive(Ns + 1)
   endassociate
   endfunction primitive2conservative

   pure function conservative2primitive(self, conservative) result(primitive)
   !< Convert conservative variables to primitive variables.
   class(integrand_euler_1D), intent(in) :: self                 !< Euler field.
   real(R_P),                 intent(in) :: conservative(:)      !< Conservative variables.
   real(R_P)                             :: primitive(1:self%Np) !< Primitive variables.
   real(R_P), allocatable                :: c(:)                 !< Species concentration.

   associate(Ns=>self%Ns, cp0=>self%cp0, cv0=>self%cv0)
      primitive(1:Ns) = conservative(1:Ns)
      primitive(Ns+3) = sum(conservative(1:Ns))
      c = primitive(1:Ns) / primitive(Ns + 3)
      primitive(Ns+4) = dot_product(c, cp0) / dot_product(c, cv0)
      primitive(Ns+1) = conservative(Ns + 1) / primitive(Ns + 3)
      primitive(Ns+2) = (conservative(Ns + 2) - 0.5_R_P * primitive(Ns + 3) * primitive(Ns + 1) * primitive(Ns + 1)) * &
                        (primitive(Ns + 4) - 1._R_P)
   endassociate
   endfunction conservative2primitive

   pure subroutine impose_boundary_conditions(self, primitive)
   !< Impose boundary conditions.
   class(integrand_euler_1D), intent(in)    :: self                                           !< Euler field.
   real(R_P),                 intent(INOUT) :: primitive(1:self%Np,1-self%Ng:self%Ni+self%Ng) !< Primitive variables.
   integer(I_P)                             :: i                                              !< Space counter.

   select case(trim(adjustl(self%BC_L)))
      case('TRA') ! trasmissive (non reflective) BC
        do i=1-self%Ng, 0
           primitive(:, i) = primitive(:, -i+1)
        enddo
      case('REF') ! reflective BC
        do i=1-self%Ng, 0
           primitive(:,           i) =  primitive(:,           -i+1) ! all variables
           primitive(self%Ns + 1, i) = -primitive(self%Ns + 1, -i+1) ! only velocity
        enddo
   endselect

   select case(trim(adjustl(self%BC_R)))
      case('TRA') ! trasmissive (non reflective) BC
        do i=self%Ni+1, self%Ni+self%Ng
           primitive(:, i) = primitive(:, self%Ni-(i-self%Ni-1))
        enddo
      case('REF') ! reflective BC
        do i=self%Ni+1, self%Ni+self%Ng
           primitive(:,           i) =  primitive(:,           self%Ni-(i-self%Ni-1)) ! all variables
           primitive(self%Ns + 1, i) = -primitive(self%Ns + 1, self%Ni-(i-self%Ni-1)) ! only velocity
        enddo
   endselect
   endsubroutine impose_boundary_conditions

   subroutine reconstruct_interfaces_states(self, primitive, r_primitive)
   !< Reconstruct the interfaces states (into primitive variables formulation) by the requested order of accuracy.
   class(integrand_euler_1D), intent(in)    :: self                                            !< Euler field.
   real(R_P),                 intent(in)    :: primitive(1:self%Np, 1-self%Ng:self%Ni+self%Ng) !< Primitive variables.
   real(R_P),                 intent(inout) :: r_primitive(1:self%Np, 1:2, 0:self%Ni+1)        !< Reconstructed primitive variables.
   real(R_P)                                :: C(1:2, 1-self%Ng:-1+self%Ng, 1:self%Ns+2)       !< Pseudo characteristic variables.
   real(R_P)                                :: CR(1:self%Ns+2, 1:2)                            !< Pseudo characteristic reconst. vars.
   real(R_P)                                :: Pm(1:self%Np, 1:2)                              !< Mean of primitive variables.
   real(R_P)                                :: LPm(1:self%Ns+2, 1:self%Ns+2, 1:2)              !< Mean left eigenvectors matrix.
   real(R_P)                                :: RPm(1:self%Ns+2, 1:self%Ns+2, 1:2)              !< Mean right eigenvectors matrix.
   integer(I_P)                             :: i, j, f, v                                      !< Counter.

   select case(self%weno_order)
   case(1) ! 1st order piecewise constant reconstruction
      do i=0, self%Ni+1
         r_primitive(:, 1, i) = primitive(:, i)
         r_primitive(:, 2, i) = r_primitive(:, 1, i)
      enddo
   case(3, 5, 7) ! 3rd, 5th or 7th order WENO reconstruction
      do i=0, self%Ni+1
         ! trasform primitive variables to pseudo charteristic ones
         do f=1, 2
            Pm(:,f) = 0.5_R_P * (primitive(:, i+f-2) + primitive(:, i+f-1))
         enddo
         do f=1, 2
            LPm(:, :, f) = eigen_vect_L(Ns=self%Ns, Np=self%Np, primitive=Pm(:, f))
            RPm(:, :, f) = eigen_vect_R(Ns=self%Ns, Np=self%Np, primitive=Pm(:, f))
         enddo
         do j=i+1-self%Ng, i-1+self%Ng
            do f=1, 2
              do v=1, self%Ns+2
                C(f, j-i, v) = dot_product(LPm(v, 1:self%Ns+2, f), primitive(1:self%Ns+2, j))
              enddo
            enddo
         enddo
         ! compute WENO reconstruction of pseudo charteristic variables
         do v=1, self%Ns+2
            call self%interpolator%interpolate(stencil=C(1:2, 1-self%Ng:-1+self%Ng, v), &
                                               interpolation=CR(v, 1:2))
         enddo
         ! trasform back reconstructed pseudo charteristic variables to primitive ones
         do f=1, 2
            do v=1, self%Ns+2
              r_primitive(v, f, i) = dot_product(RPm(v, 1:self%Ns+2, f), CR(1:self%Ns+2, f))
            enddo
            r_primitive(self%Ns+3, f, i) = sum(r_primitive(1:self%Ns, f, i))
            r_primitive(self%Ns+4, f, i) = dot_product(r_primitive(1:self%Ns, f, i) / r_primitive(self%Ns+3, f, i), self%cp0) / &
                                           dot_product(r_primitive(1:self%Ns, f, i) / r_primitive(self%Ns+3, f, i), self%cv0)
         enddo
      enddo
   endselect
   contains
     pure function eigen_vect_L(Ns, Np, primitive) result(L)
     !< Compute left eigenvectors from primitive variables.
     integer(I_P), intent(in) :: Ns               !< Number of initial species.
     integer(I_P), intent(in) :: Np               !< Number of primitive variables.
     real(R_P),    intent(in) :: primitive(1:Np)  !< Primitive variables.
     real(R_P)                :: L(1:Ns+2,1:Ns+2) !< Left eigenvectors matrix.
     real(R_P)                :: gp               !< g*p.
     real(R_P)                :: gp_a             !< g*p/a.
     integer(I_P)             :: i                !< Counter.
     integer(I_P)             :: s                !< Counter.

     gp   = primitive(Ns+4) * primitive(Ns+2)
     gp_a = gp/a(p=primitive(Ns+2), r=primitive(Ns+3), g=primitive(Ns+4))
     L = 0._R_P
                              L(1,    Ns+1) = -gp_a              ; L(1,    Ns+2) =  1._R_P
     do s=2, Ns+1
        if (primitive(s-1)>0) L(s,    s-1 ) =  gp/primitive(s-1) ; L(s,    Ns+2) = -1._R_P
     enddo
                              L(Ns+2, Ns+1) =  gp_a              ; L(Ns+2, Ns+2) =  1._R_P
     return
     endfunction eigen_vect_L

     pure function eigen_vect_R(Ns, Np, primitive) result(R)
     !< Compute right eigenvectors from primitive variables.
     integer(I_P), intent(in) :: Ns               !< Number of initial species.
     integer(I_P), intent(in) :: Np               !< Number of primitive variables.
     real(R_P),    intent(in) :: primitive(1:Np)  !< Primitive variables.
     real(R_P)                :: R(1:Ns+2,1:Ns+2) !< Right eigenvectors matrix.
     real(R_P)                :: gp               !< g*p.
     real(R_P)                :: ss               !< Speed of sound, sqrt(g*p/r).
     real(R_P)                :: gp_inv           !< 1/(g*p).
     integer(I_P)             :: i                !< Counter.
     integer(I_P)             :: s                !< Counter.

     gp = primitive(Ns+4) * primitive(Ns+2)
     ss = a(p=primitive(Ns+2), r=primitive(Ns+3), g=primitive(Ns+4))
     gp_inv = 1._R_P/gp
     R = 0._R_P
     do s=1, Ns
        R(s,    1) =  0.5_R_P*primitive(s) * gp_inv ; R(s, s+1) = primitive(s) * gp_inv ; R(s,    Ns+2) = R(s, 1)
     enddo
        R(Ns+1, 1) = -0.5_R_P* ss *gp_inv           ;                                     R(Ns+1, Ns+2) = 0.5_R_P* ss * gp_inv
        R(Ns+2, 1) =  0.5_R_P                       ;                                     R(Ns+2, Ns+2) = 0.5_R_P
     return
     endfunction eigen_vect_R
  endsubroutine reconstruct_interfaces_states

   pure subroutine riemann_solver(self, p1, r1, u1, g1, p4, r4, u4, g4, F)
   !< Solve the Riemann problem between the state $1$ and $4$ using the (local) Lax Friedrichs (Rusanov) solver.
   class(integrand_euler_1D), intent(in)  :: self         !< Euler field.
   real(R_P),                 intent(in)  :: p1           !< Pressure of state 1.
   real(R_P),                 intent(in)  :: r1           !< Density of state 1.
   real(R_P),                 intent(in)  :: u1           !< Velocity of state 1.
   real(R_P),                 intent(in)  :: g1           !< Specific heats ratio of state 1.
   real(R_P),                 intent(in)  :: p4           !< Pressure of state 4.
   real(R_P),                 intent(in)  :: r4           !< Density of state 4.
   real(R_P),                 intent(in)  :: u4           !< Velocity of state 4.
   real(R_P),                 intent(in)  :: g4           !< Specific heats ratio of state 4.
   real(R_P),                 intent(out) :: F(1:self%Nc) !< Resulting fluxes.
   real(R_P)                              :: F1(1:3)      !< State 1 fluxes.
   real(R_P)                              :: F4(1:3)      !< State 4 fluxes.
   real(R_P)                              :: u            !< Velocity of the intermediate states.
   real(R_P)                              :: p            !< Pressure of the intermediate states.
   real(R_P)                              :: S1           !< Maximum wave speed of state 1 and 4.
   real(R_P)                              :: S4           !< Maximum wave speed of state 1 and 4.
   real(R_P)                              :: lmax         !< Maximum wave speed estimation.

   ! evaluating the intermediates states 2 and 3 from the known states U1,U4 using the PVRS approximation
   call compute_inter_states(p1 = p1, r1 = r1, u1 = u1, g1 = g1, p4 = p4, r4 = r4, u4 = u4, g4 = g4, p = p, S = u, S1 = S1, S4 = S4)
   ! evalutaing the maximum waves speed
   lmax = max(abs(S1), abs(u), abs(S4))
   ! computing the fluxes of state 1 and 4
   F1 = fluxes(p = p1, r = r1, u = u1, g = g1)
   F4 = fluxes(p = p4, r = r4, u = u4, g = g4)
   ! computing the Lax-Friedrichs fluxes approximation
   F(1)         = 0.5_R_P*(F1(1) + F4(1) - lmax*(r4                        - r1                       ))
   F(self%Ns+1) = 0.5_R_P*(F1(2) + F4(2) - lmax*(r4*u4                     - r1*u1                    ))
   F(self%Ns+2) = 0.5_R_P*(F1(3) + F4(3) - lmax*(r4*E(p=p4,r=r4,u=u4,g=g4) - r1*E(p=p1,r=r1,u=u1,g=g1)))
   contains
      pure function fluxes(p, r, u, g) result(Fc)
      !< 1D Euler fluxes from primitive variables.
      real(R_P), intent(in) :: p       !< Pressure.
      real(R_P), intent(in) :: r       !< Density.
      real(R_P), intent(in) :: u       !< Velocity.
      real(R_P), intent(in) :: g       !< Specific heats ratio.
      real(R_P)             :: Fc(1:3) !< State fluxes.

      Fc(1) = r*u
      Fc(2) = Fc(1)*u + p
      Fc(3) = Fc(1)*H(p=p, r=r, u=u, g=g)
      endfunction fluxes
   endsubroutine riemann_solver

   ! non type-bound procedures
   pure subroutine compute_inter_states(r1, p1, u1, g1, r4, p4, u4, g4, p, S, S1, S4)
   !< Compute inter states (23*-states) from state1 and state4.
   real(R_P), intent(in)  :: r1             !< Density of state 1.
   real(R_P), intent(in)  :: p1             !< Pressure of state 1.
   real(R_P), intent(in)  :: u1             !< Velocity of state 1.
   real(R_P), intent(in)  :: g1             !< Specific heat ratio of state 1.
   real(R_P), intent(in)  :: r4             !< Density of state 4.
   real(R_P), intent(in)  :: p4             !< Pressure of state 4.
   real(R_P), intent(in)  :: u4             !< Velocity of state 4.
   real(R_P), intent(in)  :: g4             !< Specific heat ratio of state 4.
   real(R_P), intent(out) :: p              !< Pressure of the intermediate states.
   real(R_P), intent(out) :: S              !< Contact discontinuity signal velocity.
   real(R_P), intent(out) :: S1             !< Left fastest signal velocity.
   real(R_P), intent(out) :: S4             !< Right fastest signal velocity.
   real(R_P)              :: a1             !< Speed of sound of state 1.
   real(R_P)              :: a4             !< Speed of sound of state 4.
   real(R_P)              :: ram            !< Mean value of rho*a.
   real(R_P), parameter   :: toll=1e-10_R_P !< Tollerance.

   ! evaluation of the intermediate states pressure and velocity
   a1  = sqrt(g1 * p1 / r1)                              ! left speed of sound
   a4  = sqrt(g4 * p4 / r4)                              ! right speed of sound
   ram = 0.5_R_P * (r1 + r4) * 0.5_R_P * (a1 + a4)       ! product of mean density for mean speed of sound
   S   = 0.5_R_P * (u1 + u4) - 0.5_R_P * (p4 - p1) / ram ! evaluation of the contact wave speed (velocity of intermediate states)
   p   = 0.5_R_P * (p1 + p4) - 0.5_R_P * (u4 - u1) * ram ! evaluation of the pressure of the intermediate states
   ! evaluation of the left wave speeds
   if (p<=p1*(1._R_P + toll)) then
      ! rarefaction
      S1 = u1 - a1
   else
      ! shock
      S1 = u1 - a1 * sqrt(1._R_P + (g1 + 1._R_P) / (2._R_P * g1) * (p / p1 - 1._R_P))
   endif
   ! evaluation of the right wave speeds
   if (p<=p4 * (1._R_P + toll)) then
      ! rarefaction
      S4 = u4 + a4
   else
      ! shock
      S4 = u4 + a4 * sqrt(1._R_P + (g4 + 1._R_P) / (2._R_P * g4) * ( p / p4 - 1._R_P))
   endif
   endsubroutine compute_inter_states

  elemental function p(r, a, g) result(pressure)
  !< Compute the pressure for an ideal calorically perfect gas.
  real(R_P), intent(in) :: r        !< Density.
  real(R_P), intent(in) :: a        !< Speed of sound.
  real(R_P), intent(in) :: g        !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P)             :: pressure !< Pressure.

  pressure = r*a*a/g
  endfunction p

  elemental function r(p, a, g) result(density)
  !< Compute the density for an ideal calorically perfect gas.
  real(R_P), intent(in) :: p       !< Pressure.
  real(R_P), intent(in) :: a       !< Speed of sound.
  real(R_P), intent(in) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P)             :: density !< Density.

  density = g*p/(a*a)
  endfunction r

  elemental function a(p, r, g) result(ss)
  !< Compute the speed of sound for an ideal calorically perfect gas.

  real(R_P), intent(in) :: p  !< Pressure.
  real(R_P), intent(in) :: r  !< Density.
  real(R_P), intent(in) :: g  !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P)             :: ss !< Speed of sound.

  ss = sqrt(g*p/r)
  endfunction a

  elemental function E(p, r, u, g) result(energy)
  !< Compute total specific energy (per unit of mass).
  !<$$
  !<  E = \frac{p}{{\left( {\g  - 1} \right)\r }} + \frac{{u^2 }}{2}
  !<$$
  real(R_P), intent(in) :: p      !< Pressure.
  real(R_P), intent(in) :: r      !< Density.
  real(R_P), intent(in) :: u      !< Module of velocity vector.
  real(R_P), intent(in) :: g      !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P)             :: energy !< Total specific energy (per unit of mass).

  energy = p/((g - 1._R_P) * r) + 0.5_R_P * u * u
  endfunction E

  elemental function H(p, r, u, g) result(entalpy)
  !< Compute total specific entalpy (per unit of mass).
  !<$$
  !<  H = \frac{{\g p}}{{\left( {\g  - 1} \right)\r }} + \frac{{u^2 }}{2}
  !<$$
  real(R_P), intent(in) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P), intent(in) :: p       !< Pressure.
  real(R_P), intent(in) :: r       !< Density.
  real(R_P), intent(in) :: u       !< Module of velocity vector.
  real(R_P)             :: entalpy !< Total specific entalpy (per unit of mass).

  entalpy = g * p / ((g - 1._R_P) * r) + 0.5_R_P * u * u
  endfunction H
endmodule foodie_integrand_euler_1D
