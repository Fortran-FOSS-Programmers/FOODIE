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
  private
  integer(I_P)              :: steps=0    !< Number of time steps stored.
  integer(I_P)              :: ord=0      !< Space accuracy formal order.
  integer(I_P)              :: Ni=0       !< Space dimension.
  integer(I_P)              :: Ng=0       !< Number of ghost cells for boundary conditions handling.
  integer(I_P)              :: Ns=0       !< Number of initial species.
  integer(I_P)              :: Nc=0       !< Number of conservative variables, Ns+2.
  integer(I_P)              :: Np=0       !< Number of primitive variables, Ns+4.
  real(R_P)                 :: Dx=0._R_P  !< Space step.
  real(R_P),    allocatable :: U(:,:,:)   !< Conservative (state) variables [1:Nc,1-Ng:Ni+Ng,1:steps].
  real(R_P),    allocatable :: cp0(:)     !< Specific heat cp of initial species [1:Ns].
  real(R_P),    allocatable :: cv0(:)     !< Specific heat cv of initial species [1:Ns].
  character(:), allocatable :: BC_L       !< Left boundary condition type.
  character(:), allocatable :: BC_R       !< Right boundary condition type.
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
    procedure, pass(self), private :: left_eigen_matrix          !< Compute left eigen matrix.
    procedure, pass(self), private :: right_eigen_matrix         !< Compute right eigen matrix.
    procedure, pass(self), private :: primitive2conservative     !< Convert primitive variables to conservative ones.
    procedure, pass(self), private :: conservative2primitive     !< Convert conservative variables to primitive ones.
    procedure, pass(self), private :: speed_of_sound             !< Compute speed of sound.
    procedure, pass(self), private :: impose_boundary_conditions !< Impose boundary conditions.
    procedure, nopass,     private :: compute_inter_states       !< Compute inter states (23*).
    procedure, pass(self), private :: riemann_solver             !< Solve the Riemann Problem at cell interfaces.
    final                          :: finalize                   !< Finalize field.
endtype euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  pure subroutine init(self, Ni, Ng, Ns, Dx, BC_L, BC_R, initial_state, cp0, cv0, steps, ord)
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
  self%steps = 1 ; if (present(steps)) self%steps = steps
  self%ord = 1 ; if (present(ord)) self%ord = ord
  self%Ni = Ni
  self%Ng = Ng
  self%Ns = Ns
  self%Nc = Ns + 2
  self%Np = Ns + 4
  self%Dx = Dx
  if (allocated(self%U)) deallocate(self%U) ; allocate(self%U  (1:self%Nc, 1-Ng:Ni+Ng, 1:self%steps))
  self%cp0 = cp0
  self%cv0 = cv0
  self%BC_L = BC_L
  self%BC_R = BC_R
  do i=1, Ni
    do s=1, self%steps
      self%U(:, i, s) = self%primitive2conservative(initial_state(:, i))
    enddo
  enddo
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
  class(euler_1D), intent(IN)            :: self  !< Lorenz field.
  real(R_P), dimension(:,:), allocatable :: state !< Lorenz state vector.
  integer(I_P)                           :: i     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(state(1:self%Np, 1:self%Ni))
  do i=1, self%Ni
    state(:, i) = self%conservative2primitive(self%U(:, i, self%steps))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  pure function compute_dt(self, Nmax, Tmax, t, CFL) result(Dt)
  !--------------------------------------------------------------------------------------------------------------------------------
  !< Compute the current time step, by means of CFL condition.
  !--------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN) :: self  !< Euler field.
  integer(I_P),    intent(IN) :: Nmax  !< Maximun number of iterates.
  real(R_P),       intent(IN) :: Tmax  !< Maximum time (ignored if Nmax>0).
  real(R_P),       intent(IN) :: t     !< Time.
  real(R_P),       intent(IN) :: CFL   !< CFL value.
  real(R_P)                   :: Dt    !< Time step.
  real(R_P)                   :: vmax  !< Maximum signal speed propagation.
  real(R_P), allocatable      :: Pl(:) !< Primitive variables, left cell.
  real(R_P), allocatable      :: P(:)  !< Primitive variables.
  real(R_P), allocatable      :: Pr(:) !< Primitive variables, right cell.
  real(R_P)                   :: vmean !< Cell mean speed.
  real(R_P)                   :: a     !< Cell speed of sound.
  integer(I_P)                :: i     !< Cell counter.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  associate(Ni=>self%Ni, Dx=>self%Dx, Ns=>self%Ns, U=>self%U)
    ! evaluating max speed propagation
    vmax = 0.0_R_P
    do i=1, Ni
      Pl      = self%conservative2primitive(U(:, i-1, self%steps))
      P       = self%conservative2primitive(U(:, i,   self%steps))
      Pr      = self%conservative2primitive(U(:, i+1, self%steps))
      vmean   = (Pl(Ns + 1) + P(Ns + 1) + Pr(Ns + 1)) / 3._R_P
      a       = self%speed_of_sound(P)
      vmax    = max(abs(vmean) + a, vmax)
    enddo
    ! evaluating time step
    Dt = Dx * CFL / vmax
    if (Nmax<=0) then
      if ((t + Dt) > Tmax) Dt = Tmax - t
    endif
    return
  endassociate
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_dt

  pure function dEuler_dt(self, n) result(dState_dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Time derivative of Euler field, the residuals function.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D),        intent(IN) :: self       !< Lorenz field.
  integer(I_P), optional, intent(IN) :: n          !< Time level.
  class(integrand), allocatable      :: dState_dt  !< Euler field time derivative.
  integer(I_P)                       :: dn         !< Time level, dummy variable.
  real(R_P), allocatable             :: F(:,:)     !< Fluxes of conservative variables.
  real(R_P), allocatable             :: P(:,:)     !< Primitive variables.
  real(R_P), allocatable             :: PR(:,:,:)  !< Left (1) and right (2) interface values of primitive variables.
  real(R_P), allocatable             :: Pm(:,:)    !< Mean of primitive variables.
  real(R_P), allocatable             :: LPm(:,:,:) !< Mean left eigenvectors matrix.
  real(R_P), allocatable             :: RPm(:,:,:) !< Mean right eigenvectors matrix.
  real(R_P), allocatable             :: C (:,:,:)  !< Left and right (1,2) interface characteristic variables.
  real(R_P), allocatable             :: CR(:,:)    !< Left and right (1,2) interface reconstructed characteristic variables.
  integer(I_P)                       :: i, j, k, v !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate (Ni=>self%Ni, Ng=>self%Ng, Nc=>self%Nc, Np=>self%Np, Ns=>self%Ns, U=>self%U, cp0=>self%cp0, cv0=>self%cv0, &
             ord=>self%ord)
    dn = self%steps ; if (present(n)) dn = n
    ! allocate temporary arrays
    allocate(F(1:Nc, 0:Ni)) ; F = 0._R_P
    allocate(P(1:Np, 1-Ng:Ni+Ng)) ; P = 0._R_P
    allocate(PR(1:Np, 1:2, 1-Ng:Ni+Ng)) ; PR = 0._R_P
    allocate(Pm(1:Np, 1:2)) ; Pm = 0._R_P
    allocate(LPm(1:Ns+2, 1:Ns+2, 1:2)) ; LPm = 0._R_P
    allocate(RPm(1:Ns+2, 1:Ns+2, 1:2)) ; RPm = 0._R_P
    allocate(C(1:2, 1-ord:-1+ord, 1:Ns+2)) ; C = 0._R_P
    allocate(CR(1:Ns+2, 1:2)) ; CR = 0._R_P
    ! compute primitive variables
    do i=1 - Ng, Ni + Ng
      P(:, i) = self%conservative2primitive(U(:, i, dn))
    enddo
    call self%impose_boundary_conditions(primitive=P)
    ! compute the left and right states or Riemann Problems
    do i=0, Ni + 1
      ! trasform primitive variables to local charteristic ones
      ! compute mean of primitive variables across the interfaces i+-1/2
      do k=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          Pm(:, k) = 0.5_R_P*(P(:, i + k - 2) + P(:, i + k - 1))
      enddo
      ! compute mean left and right eigenvectors matrix
      do k=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        call self%left_eigen_matrix(primitive=Pm(:, k), matrix=LPm(:, :, k))
        call self%right_eigen_matrix(primitive=Pm(:, k), matrix=RPm(:, :, k))
      enddo
      ! transform variables into local characteristic field
      do j=i + 1 - Ng, i - 1 + Ng
        do k=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          do v=1, Ns + 2
            C(k, j - i, v) = dot_product(LPm(v, 1:Ns + 2, k), P(1:Ns + 2, j))
          enddo
        enddo
      enddo
      ! compute WENO reconstruction
      do v=1, Ns + 2
        ! compute the convultion of reconstructing plynomials
        ! CR(v, 1:2) = weno(R=Ng)
      enddo
      ! trasform back local charteristic variables to primitive ones
      do k=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        do v=1, Ns + 2
          PR(v, k, i) = dot_product(RPm(v, 1:Ns + 2, k), CR(1:Ns + 2, k))
        enddo
        PR(Ns + 3, k, i) = sum(PR(1:Ns, k, i))
        PR(Ns + 4, k, i) = dot_product(PR(1:Ns, k, i) / PR(Ns + 3, k, i), cp0)/dot_product(PR(1:Ns, k, i) / PR(Ns + 3, k, i), cv0)
      enddo
    enddo
    ! solve Rimeann Problems
    do i=0, Ni
      F(:, i) = self%riemann_solver(state1 = PR(:, 2, i    ), & ! right value of cell i (letf of interface i+1/2)
                                    state4 = PR(:, 1, i + 1))   ! left value of cell i+1 (right of interface i+1/2)
    enddo
    ! compute residuals
    allocate(euler_1D :: dState_dt)
    select type(dState_dt)
    class is(euler_1D)
      dState_dt = self
      do i=1, Ni
        dState_dt%U(:, i, dn) = (F(:, i - 1) - F(:, i)) / self%Dx
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
  do s=1, self%steps - 1
    self%U(:, :, s) = self%U(:, :, s + 1)
  enddo
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
    select type(rhs)
    class is (euler_1D)
      opr%steps = lhs%steps
      opr%ord   = lhs%ord
      opr%Ni    = lhs%Ni
      opr%Ng    = lhs%Ng
      opr%Ns    = lhs%Ns
      opr%Nc    = lhs%Nc
      opr%Np    = lhs%Np
      opr%Dx    = lhs%Dx
      opr%U     = lhs%U * rhs%U
      opr%cp0   = lhs%cp0
      opr%cv0   = lhs%cv0
      opr%BC_L  = lhs%BC_L
      opr%BC_R  = lhs%BC_R
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
    opr%steps = lhs%steps
    opr%ord   = lhs%ord
    opr%Ni    = lhs%Ni
    opr%Ng    = lhs%Ng
    opr%Ns    = lhs%Ns
    opr%Nc    = lhs%Nc
    opr%Np    = lhs%Np
    opr%Dx    = lhs%Dx
    opr%U     = lhs%U * rhs
    opr%cp0   = lhs%cp0
    opr%cv0   = lhs%cv0
    opr%BC_L  = lhs%BC_L
    opr%BC_R  = lhs%BC_R
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
    opr%steps = rhs%steps
    opr%ord   = rhs%ord
    opr%Ni    = rhs%Ni
    opr%Ng    = rhs%Ng
    opr%Ns    = rhs%Ns
    opr%Nc    = rhs%Nc
    opr%Np    = rhs%Np
    opr%Dx    = rhs%Dx
    opr%U     = rhs%U * lhs
    opr%cp0   = rhs%cp0
    opr%cv0   = rhs%cv0
    opr%BC_L  = rhs%BC_L
    opr%BC_R  = rhs%BC_R
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction real_multiply_euler

  pure function add_euler(lhs, rhs) result(opr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Add two Lorenz fields.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D),  intent(IN)  :: lhs !< Left hand side.
  class(integrand), intent(IN)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate (euler_1D :: opr)
  select type(opr)
  class is(euler_1D)
    select type(rhs)
    class is (euler_1D)
      opr%steps = lhs%steps
      opr%ord   = lhs%ord
      opr%Ni    = lhs%Ni
      opr%Ng    = lhs%Ng
      opr%Ns    = lhs%Ns
      opr%Nc    = lhs%Nc
      opr%Np    = lhs%Np
      opr%Dx    = lhs%Dx
      opr%U     = lhs%U + rhs%U
      opr%cp0   = lhs%cp0
      opr%cv0   = lhs%cv0
      opr%BC_L  = lhs%BC_L
      opr%BC_R  = lhs%BC_R
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
  class is (euler_1D)
    lhs%steps = rhs%steps
    lhs%ord   = rhs%ord
    lhs%Ni    = rhs%Ni
    lhs%Ng    = rhs%Ng
    lhs%Ns    = rhs%Ns
    lhs%Nc    = rhs%Nc
    lhs%Np    = rhs%Np
    lhs%Dx    = rhs%Dx
    lhs%U     = rhs%U
    lhs%cp0   = rhs%cp0
    lhs%cv0   = rhs%cv0
    lhs%BC_L  = rhs%BC_L
    lhs%BC_R  = rhs%BC_R
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine euler_assign_real

  ! private methods
  pure subroutine left_eigen_matrix(self, primitive, matrix)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute left eigenvectors from primitive variables.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN)  :: self         !< Euler field.
  real(R_P),       intent(IN)  :: primitive(:) !< Primitive variables [1:Np].
  real(R_P),       intent(OUT) :: matrix(:,:)  !< Conservative variables [1:Ns+2,1:Ns+2].
  real(R_P)                    :: a            !< Speed of sound.
  real(R_P)                    :: gp           !< g*p.
  real(R_P)                    :: gp_a         !< g*p/a.
  integer(I_P)                 :: i            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(Ns=>self%Ns)
    a    = self%speed_of_sound(primitive)
    gp   = primitive(Ns+4) * primitive(Ns+2)
    gp_a = gp/a

    matrix = 0._R_P
    do i=1, Ns
      matrix(i, i)    =  gp
      matrix(i, Ns+2) = -primitive(i)
    enddo
    matrix(Ns+1, Ns+1) = -gp_a ; matrix(Ns+1, Ns+2) = 1._R_P
    matrix(Ns+2, Ns+1) =  gp_a ; matrix(Ns+2, Ns+2) = 1._R_P
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine left_eigen_matrix

  pure subroutine right_eigen_matrix(self, primitive, matrix)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute right eigenvectors from primitive variables.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN)  :: self         !< Euler field.
  real(R_P),       intent(IN)  :: primitive(:) !< Primitive variables [1:Np].
  real(R_P),       intent(OUT) :: matrix(:,:)  !< Conservative variables [1:Ns+2,1:Ns+2].
  real(R_P)                    :: a            !< Speed of sound.
  real(R_P)                    :: gp_inv       !< 1/(g*p).
  real(R_P)                    :: a_2gp        !< a/(2*g*p).
  integer(I_P)                 :: i            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(Ns=>self%Ns)
    a      = self%speed_of_sound(primitive)
    gp_inv = 1._R_P / (primitive(Ns+4) * primitive(Ns+2))
    a_2gp  = a * gp_inv * 0.5_R_P

    matrix = 0._R_P
    do i=1, Ns
      matrix(i, i)    = gp_inv
      matrix(i, Ns+1) = primitive(i) * gp_inv * 0.5_R_P
      matrix(i, Ns+2) = matrix(i, Ns+1)
    enddo
    matrix(Ns+1, Ns+1) = -a_2gp   ; matrix(Ns+1, Ns+2) = a_2gp
    matrix(Ns+2, Ns+1) =  0.5_R_P ; matrix(Ns+2, Ns+2) = 0.5_R_P
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine right_eigen_matrix

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

  pure function speed_of_sound(self, primitive) result(a)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the speed of sound from primitive variables.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN) :: self            !< Euler field.
  real(R_P),       intent(IN) :: primitive(:)    !< Primitive variables.
  real(R_P)                   :: a               !< Speed of sound.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  a = sqrt(primitive(self%Ns+4) * primitive(self%Ns+2) / primitive(self%Ns+3))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction speed_of_sound

  pure subroutine impose_boundary_conditions(self, primitive)
  !--------------------------------------------------------------------------------------------------------------------------------
  !< Impose boundary conditions.
  !<
  !< The boundary conditions are imposed on the primitive variables by means of the ghost cells approach.
  !--------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN)    :: self           !< Euler field.
  real(R_P),       intent(INOUT) :: primitive(:,:) !< Primitive variables [1:Np,1-Ng:Ni+Ng].
  integer(I_P)                   :: i              !< Space counter.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(trim(adjustl(self%BC_L)))
    case('TRA') ! trasmissive (non reflective) BC
      do i=1, self%Ng
        primitive(:, i) = primitive(:, 2 * self%Ng - (i - 1))
      enddo
    case('REF') ! reflective BC
      do i=1, self%Ng
        primitive(:,           i) =  primitive(:,           2 * self%Ng - (i - 1)) ! all variables
        primitive(self%Ns + 1, i) = -primitive(self%Ns + 1, 2 * self%Ng - (i - 1)) ! only velocity
      enddo
  endselect

  select case(trim(adjustl(self%BC_R)))
    case('TRA') ! trasmissive (non reflective) BC
      do i=1,self%Ng
        primitive(:, 2 * self%Ng - (i - 1)) = primitive(:, i)
      enddo
    case('REF') ! reflective BC
      do i=1,self%Ng
        primitive(:,           2 * self%Ng - (i - 1)) =  primitive(:,           i) ! all variables
        primitive(self%Ns + 1, 2 * self%Ng - (i - 1)) = -primitive(self%Ns + 1, i) ! only velocity
      enddo
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endsubroutine impose_boundary_conditions

  pure subroutine compute_inter_states(r1, p1, u1, g1, r4, p4, u4, g4, p, S, S1, S4)
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

  pure function riemann_solver(self, state1, state4) result(F)
  !--------------------------------------------------------------------------------------------------------------------------------
  !< Solve the Riemann problem between the state $1$ and $4$ using the HLLC solver.
  !--------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D), intent(IN) :: self         !< Euler field.
  real(R_P), intent(IN)       :: state1(:)    !< State 1.
  real(R_P), intent(IN)       :: state4(:)    !< State 4.
  real(R_P)                   :: F(1:self%Nc) !< Resulting fluxes.
  real(R_P)                   :: S1           !< Left fastest signal velocity.
  real(R_P)                   :: S4           !< Right fastest signal velocity.
  real(R_P)                   :: S            !< Contact discontinuity signal velocity.
  real(R_P)                   :: p            !< Pressure of the intermediate states.
  real(R_P)                   :: U1S1         !< Mass conservation of left-star-state.
  real(R_P)                   :: U1S4         !< Mass conservation of right-star-state.
  real(R_P)                   :: U2S1         !< Momentum conservation of left-star-state.
  real(R_P)                   :: U2S4         !< Momentum conservation of right-star-state.
  real(R_P)                   :: U3S1         !< Energy conservation of left-star-state.
  real(R_P)                   :: U3S4         !< Energy conservation of right-star-state.
  real(R_P)                   :: H1           !< Entalpy of left state.
  real(R_P)                   :: E1           !< Internal energy of left state.
  real(R_P)                   :: H4           !< Entalpy of right state.
  real(R_P)                   :: E4           !< Internal energy of right state.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  associate(F_r=>F(1), F_u=>F(self%Ns+1), F_E=>F(self%Ns+2), &
            r1=>state1(self%Ns+3), u1=>state1(self%Ns+1), p1=>state1(self%Ns+2),  g1=>state1(self%Ns+4), &
            r4=>state4(self%Ns+3), u4=>state4(self%Ns+1), p4=>state4(self%Ns+2),  g4=>state4(self%Ns+4))
    call self%compute_inter_states(r1=r1, p1=p1, u1=u1, g1=g1, r4=r4, p4=p4, u4=u4, g4=g4, p=p, S=S, S1=S1, S4=S4)
    select case(minloc([-S1, S1*S, S*S4, S4], dim=1))
    case(1)
      F_r = r1 * u1
      F_u = F_r * u1 + p1
      F_E = F_r * H(p=p1, r=r1, u=u1, g=g1)
    case(2)
      E1   = E(p=p1, r=r1, u=u1, g=g1)
      H1   = H(p=p1, r=r1, u=u1, g=g1)
      U1S1 = r1 * (S1 - u1) / (S1 - S)
      U2S1 = r1 * (S1 - u1) / (S1 - S) * S
      U3S1 = r1 * (S1 - u1) / (S1 - S) * (E1 + (S - u1) * (S + p1 / ( r1 * (S1 - u1))))

      F_r = r1 * u1           + S1 * (U1S1 - r1)
      F_u = p1 + r1 * u1 * u1 + S1 * (U2S1 - r1 * u1)
      F_E = r1 * u1 * H1      + S1 * (U3S1 - r1 * E1)
    case(3)
      E4   = E(p=p4, r=r4, u=u4, g=g4)
      H4   = H(p=p4, r=r4, u=u4, g=g4)
      U1S4 = r4 * (S4 - u4) / (S4 - S)
      U2S4 = r4 * (S4 - u4) / (S4 - S) * S
      U3S4 = r4 * (S4 - u4) / (S4 - S) * (E4 + (S - u4)*(S + p4 / (r4 * (S4 - u4))))

      F_r = r4 * u4           + S4 * (U1S4 - r4)
      F_u = p4 + r4 * u4 * u4 + S4 * (U2S4 - r4 * u4)
      F_E = r4 * u4 * H4      + S4 * (U3S4 - r4 * E4)
    case(4)
      F_r = r4 * u4
      F_u = F_r * u4 + p4
      F_E = F_r * H(p=p4, r=r4, u=u4, g=g4)
    endselect
  endassociate
  if ( F(1)>0._R_P) then
    F(:) = state1(:)/state1(self%Ns + 3) * F(1)
  else
    F(:) = state4(:)/state4(self%Ns + 3) * F(1)
  endif
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  contains

    elemental function E(p, r, u, g) result(energy)
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Compute total specific energy (per unit of mass).
    !-------------------------------------------------------------------------------------------------------------------------------
    real(R_P), intent(IN) :: p      !< Pressure.
    real(R_P), intent(IN) :: r      !< Density.
    real(R_P), intent(IN) :: u      !< Module of velocity vector.
    real(R_P), intent(IN) :: g      !< Specific heats ratio $\frac{{c_p}}{{c_v}}$.
    real(R_P)             :: energy !< Total specific energy (per unit of mass).
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    energy = p/((g - 1._R_P) * r) + 0.5_R_P * u * u
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction E

    elemental function H(p,r,u,g) result(entalpy)
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Compute total specific entalpy (per unit of mass).
    !-------------------------------------------------------------------------------------------------------------------------------
    real(R_P), intent(IN) :: g       ! Specific heats ratio $\frac{{c_p}}{{c_v}}$.
    real(R_P), intent(IN) :: p       ! Pressure.
    real(R_P), intent(IN) :: r       ! Density.
    real(R_P), intent(IN) :: u       ! Module of velocity vector.
    real(R_P)             :: entalpy ! Total specific entalpy (per unit of mass).
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    entalpy = g * p / ((g - 1._R_P) * r) + 0.5_R_P * u * u
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction H
  endfunction riemann_solver

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
endmodule type_euler_1D
