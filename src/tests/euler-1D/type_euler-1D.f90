!< Define Euler 1D field that is a concrete extension of the abstract integrand type.
module type_euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define Euler 1D field that is a concrete extension of the abstract integrand type.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P
use foodie, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(integrand) :: euler_1D
  !< Euler 1D PDEs system field.
  !<
  !< It is a FOODiE integrand class.
  !<
  !<### 1D Euler PDEs system
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
  !< Conservative  variables are organized as an array of reals which the first index means:
  !<
  !< + 1    : mass conservation of species 1    (r1)
  !< + 2    : mass conservation of species 2    (r2)
  !< + ...  :
  !< + s    : mass conservation of species s-th (rs)
  !< + ...  :
  !< + Ns   : mass conservation of species Ns   (rNs)
  !< + Ns+1 : momentum conservation             (r*u)
  !< + Ns+2 : energy conservation               (r*E)
  private
  integer(I_P)              :: steps=0         !< Number of time steps stored.
  integer(I_P)              :: ord=0           !< Space accuracy formal order.
  integer(I_P)              :: Ni=0            !< Space dimension.
  integer(I_P)              :: Ng=0            !< Number of ghost cells for boundary conditions handling.
  integer(I_P)              :: Ns=0            !< Number of initial species.
  integer(I_P)              :: Nc=0            !< Number of conservative variables, Ns+2.
  integer(I_P)              :: Np=0            !< Number of primitive variables, Ns+4.
  real(R_P)                 :: Dx=0._R_P       !< Space step.
  real(R_P),    allocatable :: U(:,:)          !< Conservative (state) variables [1:Nc,1-Ng:Ni+Ng].
  real(R_P),    allocatable :: previous(:,:,:) !< Conservative (state) variables of previous time steps [1:Nc,1-Ng:Ni+Ng,1:steps].
  real(R_P),    allocatable :: cp0(:)          !< Specific heat cp of initial species [1:Ns].
  real(R_P),    allocatable :: cv0(:)          !< Specific heat cv of initial species [1:Ns].
  character(:), allocatable :: BC_L            !< Left boundary condition type.
  character(:), allocatable :: BC_R            !< Right boundary condition type.
  contains
    ! public methods
    ! auxiliary methods
    procedure, pass(self), public :: init             !< Init field.
    procedure, pass(self), public :: destroy          !< Destroy field.
    procedure, pass(self), public :: output           !< Extract Euler field.
    procedure, pass(self), public :: dt => compute_dt !< Compute the current time step, by means of CFL condition.
    ! type_integrand deferred methods
    procedure, pass(self), public :: t => dEuler_dt        !< Time derivate, residuals function.
    procedure, pass(self), public :: update_previous_steps !< Update previous time steps.
    ! operators overloading
    procedure, pass(lhs),  public :: integrand_multiply_integrand => euler_multiply_euler !< Euler * Euler operator.
    procedure, pass(lhs),  public :: integrand_multiply_real => euler_multiply_real       !< Euler * real operator.
    procedure, pass(rhs),  public :: real_multiply_integrand => real_multiply_euler       !< Real * Euler operator.
    procedure, pass(lhs),  public :: add => add_euler                                     !< Euler + Euler oprator.
    procedure, pass(lhs),  public :: assign_integrand => euler_assign_euler               !< Euler = Euler.
    procedure, pass(lhs),  public :: assign_real => euler_assign_real                     !< Euler = real.
    ! private methods
    procedure, pass(self), private :: primitive2conservative     !< Convert primitive variables to conservative ones.
    procedure, pass(self), private :: conservative2primitive     !< Convert conservative variables to primitive ones.
    procedure, pass(self), private :: impose_boundary_conditions !< Impose boundary conditions.
    procedure, pass(self), private :: riemann_solver             !< Solve the Riemann Problem at cell interfaces.
    final                          :: finalize                   !< Finalize field.
endtype euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  subroutine init(self, Ni, Ng, Ns, Dx, BC_L, BC_R, initial_state, cp0, cv0, steps, ord)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Init field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D),        intent(INOUT) :: self               !< Euler field.
  integer(I_P),           intent(IN)    :: Ni                 !< Space dimension.
  integer(I_P),           intent(IN)    :: Ng                 !< Number of ghost cells.
  integer(I_P),           intent(IN)    :: Ns                 !< Number of initial species.
  real(R_P),              intent(IN)    :: Dx                 !< Space step.
  character(*),           intent(IN)    :: BC_L               !< Left boundary condition type.
  character(*),           intent(IN)    :: BC_R               !< Right boundary condition type.
  real(R_P),              intent(IN)    :: initial_state(:,:) !< Initial state of primitive variables.
  real(R_P),              intent(IN)    :: cp0(:)             !< Initial specific heat, constant pressure.
  real(R_P),              intent(IN)    :: cv0(:)             !< Initial specific heat, constant volume.
  integer(I_P), optional, intent(IN)    :: steps              !< Time steps stored.
  integer(I_P), optional, intent(IN)    :: ord                !< Space accuracy formal order.
  integer(I_P)                          :: i                  !< Space counter.
  integer(I_P)                          :: s                  !< Steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = 0 ; if (present(steps)) self%steps = steps
  self%ord = 1 ; if (present(ord)) self%ord = ord
  self%Ni = Ni
  self%Ng = Ng
  self%Ns = Ns
  self%Nc = Ns + 2
  self%Np = Ns + 4
  self%Dx = Dx
  if (allocated(self%U)) deallocate(self%U) ; allocate(self%U  (1:self%Nc, 1:Ni))
  if (self%steps>0) then
    if (allocated(self%previous)) deallocate(self%previous) ; allocate(self%previous(1:self%Nc, 1:Ni, 1:self%steps))
  endif
  self%cp0 = cp0
  self%cv0 = cv0
  self%BC_L = BC_L
  self%BC_R = BC_R
  do i=1, Ni
    self%U(:, i) = self%primitive2conservative(initial_state(:, i))
  enddo
  if (self%steps>0) then
    do s=1, self%steps
      do i=1, Ni
        self%previous(:, i, s) = self%primitive2conservative(initial_state(:, i))
      enddo
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  pure subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(INOUT) :: self !< Euler field.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%steps = 0
  self%ord = 0
  self%Ni = 0
  self%Ng = 0
  self%Ns = 0
  self%Nc = 0
  self%Np = 0
  self%Dx = 0._R_P
  if (allocated(self%U)) deallocate(self%U)
  if (allocated(self%previous)) deallocate(self%previous)
  if (allocated(self%cp0)) deallocate(self%cp0)
  if (allocated(self%cv0)) deallocate(self%cv0)
  if (allocated(self%BC_L)) deallocate(self%BC_L)
  if (allocated(self%BC_R)) deallocate(self%BC_R)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  pure function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Euler field state (primitive variables).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN)            :: self  !< Euler field.
  real(R_P), dimension(:,:), allocatable :: state !< Euler state vector.
  integer(I_P)                           :: i     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(state(1:self%Np, 1:self%Ni))
  do i=1, self%Ni
    state(:, i) = self%conservative2primitive(self%U(:, i))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  function compute_dt(self, Nmax, Tmax, t, CFL) result(Dt)
  !--------------------------------------------------------------------------------------------------------------------------------
  !< Compute the current time step by means of CFL condition.
  !--------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN) :: self           !< Euler field.
  integer(I_P),    intent(IN) :: Nmax           !< Maximun number of iterates.
  real(R_P),       intent(IN) :: Tmax           !< Maximum time (ignored if Nmax>0).
  real(R_P),       intent(IN) :: t              !< Time.
  real(R_P),       intent(IN) :: CFL            !< CFL value.
  real(R_P)                   :: Dt             !< Time step.
  real(R_P), allocatable      :: P(:)           !< Primitive variables.
  real(R_P)                   :: vmax           !< Maximum propagation speed of signals.
  integer(I_P)                :: i              !< Counter.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  associate(Ni=>self%Ni, Ns=>self%Ns, Dx=>self%Dx)
    vmax = 0._R_P
    do i=1, Ni
      P    = self%conservative2primitive(self%U(:, i))
      vmax = max(abs(P(Ns+1)) + a(p=P(Ns+2), r=P(Ns+3), g=P(Ns+4)), vmax)
    enddo
    Dt = Dx * CFL / vmax
    if (Nmax<=0) then
      if ((t + Dt) > Tmax) Dt = Tmax - t
    endif
    return
  endassociate
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_dt

  function dEuler_dt(self, n) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Euler field, the residuals function.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D),        intent(IN) :: self      !< Euler field.
  integer(I_P), optional, intent(IN) :: n         !< Time level.
  class(integrand), allocatable      :: dState_dt !< Euler field time derivative.
  real(R_P), allocatable             :: F(:,:)    !< Fluxes of conservative variables.
  real(R_P), allocatable             :: P(:,:)    !< Primitive variables.
  integer(I_P)                       :: i         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate (Ni=>self%Ni, Ng=>self%Ng, Nc=>self%Nc, Np=>self%Np, Ns=>self%Ns, U=>self%U)
    ! allocate temporary arrays
    allocate(F(1:Nc, 0:Ni)) ; F = 0._R_P
    allocate(P(1:Np, 1-Ng:Ni+Ng)) ; P = 0._R_P
    ! compute primitive variables
    do i=1, Ni
      P(:, i) = self%conservative2primitive(U(:, i))
    enddo
    call self%impose_boundary_conditions(primitive=P)
    ! compute fluxes by solving Rimeann Problems at each interface
    do i=0, Ni
      call self%riemann_solver(p1=P(Ns+2, i  ), r1=P(Ns+3, i  ), u1=P(Ns+1, i  ), g1=P(Ns+4, i  ), &
                               p4=P(Ns+2, i+1), r4=P(Ns+3, i+1), u4=P(Ns+1, i+1), g4=P(Ns+4, i+1), &
                               F=F(:, i))
    enddo
    ! compute residuals
    allocate(euler_1D :: dState_dt)
    select type(dState_dt)
    class is(euler_1D)
      dState_dt = self
      do i=1, Ni
        dState_dt%U(:, i) = (F(:, i - 1) - F(:, i)) / self%Dx
      enddo
    endselect
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dEuler_dt

  pure subroutine update_previous_steps(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Update previous time steps.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(INOUT) :: self !< Euler field.
  integer                        :: s    !< Time steps counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (self%steps>0) then
    do s=1, self%steps - 1
      self%previous(:, :, s) = self%previous(:, :, s + 1)
    enddo
    self%previous(:, :, self%steps) = self%U
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_previous_steps

  ! operators overloading
  pure function euler_multiply_euler(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply an Euler field by another one.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D),  intent(IN)  :: lhs !< Left hand side.
  class(integrand), intent(IN)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(euler_1D :: opr)
  select type(opr)
  class is(euler_1D)
    opr = lhs
    select type(rhs)
    class is (euler_1D)
      ! opr%steps = lhs%steps
      ! opr%ord   = lhs%ord
      ! opr%Ni    = lhs%Ni
      ! opr%Ng    = lhs%Ng
      ! opr%Ns    = lhs%Ns
      ! opr%Nc    = lhs%Nc
      ! opr%Np    = lhs%Np
      ! opr%Dx    = lhs%Dx
      opr%U     = lhs%U * rhs%U
      ! opr%cp0   = lhs%cp0
      ! opr%cv0   = lhs%cv0
      ! opr%BC_L  = lhs%BC_L
      ! opr%BC_R  = lhs%BC_R
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction euler_multiply_euler

  pure function euler_multiply_real(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply an Euler field by a real scalar.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN)   :: lhs !< Left hand side.
  real(R_P),       intent(IN)   :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(euler_1D :: opr)
  select type(opr)
  class is(euler_1D)
    opr = lhs
    ! opr%steps = lhs%steps
    ! opr%ord   = lhs%ord
    ! opr%Ni    = lhs%Ni
    ! opr%Ng    = lhs%Ng
    ! opr%Ns    = lhs%Ns
    ! opr%Nc    = lhs%Nc
    ! opr%Np    = lhs%Np
    ! opr%Dx    = lhs%Dx
    opr%U     = lhs%U * rhs
    ! opr%cp0   = lhs%cp0
    ! opr%cv0   = lhs%cv0
    ! opr%BC_L  = lhs%BC_L
    ! opr%BC_R  = lhs%BC_R
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction euler_multiply_real

  pure function real_multiply_euler(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Multiply a real scalar by an Euler field.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P),       intent(IN)   :: lhs !< Left hand side.
  class(euler_1D), intent(IN)   :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(euler_1D :: opr)
  select type(opr)
  class is(euler_1D)
    opr = rhs
    ! opr%steps = rhs%steps
    ! opr%ord   = rhs%ord
    ! opr%Ni    = rhs%Ni
    ! opr%Ng    = rhs%Ng
    ! opr%Ns    = rhs%Ns
    ! opr%Nc    = rhs%Nc
    ! opr%Np    = rhs%Np
    ! opr%Dx    = rhs%Dx
    opr%U     = rhs%U * lhs
    ! opr%cp0   = rhs%cp0
    ! opr%cv0   = rhs%cv0
    ! opr%BC_L  = rhs%BC_L
    ! opr%BC_R  = rhs%BC_R
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_multiply_euler

  pure function add_euler(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Euler fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D),  intent(IN)  :: lhs !< Left hand side.
  class(integrand), intent(IN)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate (euler_1D :: opr)
  select type(opr)
  class is(euler_1D)
    opr = lhs
    select type(rhs)
    class is (euler_1D)
      ! opr%steps = lhs%steps
      ! opr%ord   = lhs%ord
      ! opr%Ni    = lhs%Ni
      ! opr%Ng    = lhs%Ng
      ! opr%Ns    = lhs%Ns
      ! opr%Nc    = lhs%Nc
      ! opr%Np    = lhs%Np
      ! opr%Dx    = lhs%Dx
      opr%U     = lhs%U + rhs%U
      ! opr%cp0   = lhs%cp0
      ! opr%cv0   = lhs%cv0
      ! opr%BC_L  = lhs%BC_L
      ! opr%BC_R  = lhs%BC_R
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction add_euler

  pure subroutine euler_assign_euler(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one Euler field to another.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D),  intent(INOUT) :: lhs !< Left hand side.
  class(integrand), intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(rhs)
  class is(euler_1D)
                                 lhs%steps    = rhs%steps
                                 lhs%ord      = rhs%ord
                                 lhs%Ni       = rhs%Ni
                                 lhs%Ng       = rhs%Ng
                                 lhs%Ns       = rhs%Ns
                                 lhs%Nc       = rhs%Nc
                                 lhs%Np       = rhs%Np
                                 lhs%Dx       = rhs%Dx
    if (allocated(rhs%U))        lhs%U        = rhs%U
    if (allocated(rhs%previous)) lhs%previous = rhs%previous
    if (allocated(rhs%cp0))      lhs%cp0      = rhs%cp0
    if (allocated(rhs%cv0))      lhs%cv0      = rhs%cv0
    if (allocated(rhs%BC_L))     lhs%BC_L     = rhs%BC_L
    if (allocated(rhs%BC_R))     lhs%BC_R     = rhs%BC_R
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine euler_assign_euler

  pure subroutine euler_assign_real(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one real to an Euler field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(INOUT) :: lhs !< Left hand side.
  real(R_P),       intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(lhs%U)) lhs%U = rhs
  if (allocated(lhs%previous)) lhs%previous = rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine euler_assign_real

  ! private methods
  pure function primitive2conservative(self, primitive) result(conservative)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert primitive variables to conservative variables.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN) :: self                    !< Euler field.
  real(R_P),       intent(IN) :: primitive(:)            !< Primitive variables.
  real(R_P)                   :: conservative(1:self%Nc) !< Conservative variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(Ns=>self%Ns)
    conservative(1:Ns) = primitive(1:Ns)
    conservative(Ns + 1) = primitive(Ns + 3) * primitive(Ns + 1)
    conservative(Ns + 2) = primitive(Ns + 2) / (primitive(Ns + 4) - 1._R_P) + &
                                0.5_R_P*primitive(Ns + 3) * primitive(Ns + 1) * primitive(Ns + 1)
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction primitive2conservative

  pure function conservative2primitive(self, conservative) result(primitive)
  !--------------------------------------------------------------------------------------------------------------------------------
  !< Convert conservative variables to primitive variables.
  !--------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN) :: self                 !< Euler field.
  real(R_P),       intent(IN) :: conservative(:)      !< Conservative variables.
  real(R_P)                   :: primitive(1:self%Np) !< Primitive variables.
  real(R_P), allocatable      :: c(:)                 !< Species concentration.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  associate(Ns=>self%Ns, cp0=>self%cp0, cv0=>self%cv0)
    primitive(1:Ns) = conservative(1:Ns)
    primitive(Ns + 3) = sum(conservative(1:Ns))
    c = primitive(1:Ns) / primitive(Ns + 3)
    primitive(Ns + 4) = dot_product(c, cp0) / dot_product(c, cv0)
    primitive(Ns + 1) = conservative(Ns + 1) / primitive(Ns + 3)
    primitive(Ns + 2) = (conservative(Ns + 2) - 0.5_R_P * primitive(Ns + 3) * primitive(Ns + 1) * primitive(Ns + 1)) * &
                        (primitive(Ns + 4) - 1._R_P)
  endassociate
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction conservative2primitive

  subroutine impose_boundary_conditions(self, primitive)
  !--------------------------------------------------------------------------------------------------------------------------------
  !< Impose boundary conditions.
  !<
  !< The boundary conditions are imposed on the primitive variables by means of the ghost cells approach.
  !--------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN)    :: self           !< Euler field.
  real(R_P),       intent(INOUT) :: primitive(1:self%Np,1-self%Ng:self%Ni+self%Ng) !< Primitive variables [1:Np,1-Ng:Ni+Ng].
  integer(I_P)                   :: i              !< Space counter.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
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
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endsubroutine impose_boundary_conditions

  subroutine riemann_solver(self, p1, r1, u1, g1, p4, r4, u4, g4, F)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Solve the Riemann problem between the state $1$ and $4$ using the (local) Lax Friedrichs (Rusanov) solver.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN) :: self         !< Euler field.
  real(R_P), intent(IN)       :: p1           !< Pressure of state 1.
  real(R_P), intent(IN)       :: r1           !< Density of state 1.
  real(R_P), intent(IN)       :: u1           !< Velocity of state 1.
  real(R_P), intent(IN)       :: g1           !< Specific heats ratio of state 1.
  real(R_P), intent(IN)       :: p4           !< Pressure of state 4.
  real(R_P), intent(IN)       :: r4           !< Density of state 4.
  real(R_P), intent(IN)       :: u4           !< Velocity of state 4.
  real(R_P), intent(IN)       :: g4           !< Specific heats ratio of state 4.
  real(R_P), intent(OUT)      :: F(1:self%Nc) !< Resulting fluxes.
  real(R_P)                   :: F1(1:3)      !< State 1 fluxes.
  real(R_P)                   :: F4(1:3)      !< State 4 fluxes.
  real(R_P)                   :: u            !< Velocity of the intermediate states.
  real(R_P)                   :: p            !< Pressure of the intermediate states.
  real(R_P)                   :: S1           !< Maximum wave speed of state 1 and 4.
  real(R_P)                   :: S4           !< Maximum wave speed of state 1 and 4.
  real(R_P)                   :: lmax         !< Maximum wave speed estimation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    function fluxes(p, r, u, g) result(Fc)
    !-------------------------------------------------------------------------------------------------------------------------------
    !< 1D Euler fluxes from primitive variables.
    !-------------------------------------------------------------------------------------------------------------------------------
    real(R_P), intent(IN) :: p       !< Pressure.
    real(R_P), intent(IN) :: r       !< Density.
    real(R_P), intent(IN) :: u       !< Velocity.
    real(R_P), intent(IN) :: g       !< Specific heats ratio.
    real(R_P)             :: Fc(1:3) !< State fluxes.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    Fc(1) = r*u
    Fc(2) = Fc(1)*u + p
    Fc(3) = Fc(1)*H(p=p, r=r, u=u, g=g)
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction fluxes
  endsubroutine riemann_solver

  subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy field.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(euler_1D), intent(INOUT) :: self !< Euler field.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  ! non type-bound procedures
  subroutine compute_inter_states(r1, p1, u1, g1, r4, p4, u4, g4, p, S, S1, S4)
  !------------------------------------------------------------------------------------------------------------------------------
  !< Compute inter states (23*-states) from state1 and state4.
  !------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN)  :: r1             !< Density of state 1.
  real(R_P), intent(IN)  :: p1             !< Pressure of state 1.
  real(R_P), intent(IN)  :: u1             !< Velocity of state 1.
  real(R_P), intent(IN)  :: g1             !< Specific heat ratio of state 1.
  real(R_P), intent(IN)  :: r4             !< Density of state 4.
  real(R_P), intent(IN)  :: p4             !< Pressure of state 4.
  real(R_P), intent(IN)  :: u4             !< Velocity of state 4.
  real(R_P), intent(IN)  :: g4             !< Specific heat ratio of state 4.
  real(R_P), intent(OUT) :: p              !< Pressure of the intermediate states.
  real(R_P), intent(OUT) :: S              !< Contact discontinuity signal velocity.
  real(R_P), intent(OUT) :: S1             !< Left fastest signal velocity.
  real(R_P), intent(OUT) :: S4             !< Right fastest signal velocity.
  real(R_P)              :: a1             !< Speed of sound of state 1.
  real(R_P)              :: a4             !< Speed of sound of state 4.
  real(R_P)              :: ram            !< Mean value of rho*a.
  real(R_P), parameter   :: toll=1e-10_R_P !< Tollerance.
  !------------------------------------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------------------------------------
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
  return
  !------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_inter_states

  elemental function p(r, a, g) result(pressure)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the pressure for an ideal calorically perfect gas.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN) :: r        !< Density.
  real(R_P), intent(IN) :: a        !< Speed of sound.
  real(R_P), intent(IN) :: g        !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P)             :: pressure !< Pressure.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pressure = r*a*a/g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction p

  elemental function r(p, a, g) result(density)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the density for an ideal calorically perfect gas.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN) :: p       !< Pressure.
  real(R_P), intent(IN) :: a       !< Speed of sound.
  real(R_P), intent(IN) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P)             :: density !< Density.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  density = g*p/(a*a)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction r

  elemental function a(p, r, g) result(ss)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the speed of sound for an ideal calorically perfect gas.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN) :: p  !< Pressure.
  real(R_P), intent(IN) :: r  !< Density.
  real(R_P), intent(IN) :: g  !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P)             :: ss !< Speed of sound.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ss = sqrt(g*p/r)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction a

  elemental function E(p, r, u, g) result(energy)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute total specific energy (per unit of mass).
  !<$$
  !<  E = \frac{p}{{\left( {\g  - 1} \right)\r }} + \frac{{u^2 }}{2}
  !<$$
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN) :: p      !< Pressure.
  real(R_P), intent(IN) :: r      !< Density.
  real(R_P), intent(IN) :: u      !< Module of velocity vector.
  real(R_P), intent(IN) :: g      !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P)             :: energy !< Total specific energy (per unit of mass).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  energy = p/((g - 1._R_P) * r) + 0.5_R_P * u * u
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction E

  elemental function H(p, r, u, g) result(entalpy)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute total specific entalpy (per unit of mass).
  !<$$
  !<  H = \frac{{\g p}}{{\left( {\g  - 1} \right)\r }} + \frac{{u^2 }}{2}
  !<$$
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R_P), intent(IN) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
  real(R_P), intent(IN) :: p       !< Pressure.
  real(R_P), intent(IN) :: r       !< Density.
  real(R_P), intent(IN) :: u       !< Module of velocity vector.
  real(R_P)             :: entalpy !< Total specific entalpy (per unit of mass).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  entalpy = g * p / ((g - 1._R_P) * r) + 0.5_R_P * u * u
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction H
endmodule type_euler_1D
