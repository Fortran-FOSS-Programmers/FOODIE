!< Define Euler 1D CAF (global) field that is based on Euler 1D local field.
module type_euler_1D_caf
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define Euler 1D global (CAF) field that is based on Euler 1D local field.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P
use type_euler_1D, only : euler_1D
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: euler_1D_caf
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: euler_1D_caf
  private
  integer(I_P)                        :: ord=0    !< Space accuracy formal order.
  integer(I_P)                        :: Ni=0     !< Space dimension.
  integer(I_P)                        :: Ng=0     !< Number of ghost cells for boundary conditions handling.
#ifdef CAF
  type(euler_1D), allocatable, public :: local[:] !< Local (CAF image) domain.
#else
  type(euler_1D), allocatable, public :: local    !< Local (CAF image) domain.
#endif
  contains
    ! auxiliary methods
    procedure, pass(self), public :: init        !< Init field.
    procedure, pass(self), public :: destroy     !< Destroy field.
    procedure, pass(self), public :: output      !< Extract Euler field.
    procedure, pass(self), public :: synchronize !< Synchronize CAF images.
endtype euler_1D_caf
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! auxiliary methods
  subroutine init(self, Ni, Ns, Dx, BC_L, BC_R, initial_state, cp0, cv0, ord)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Init field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D_caf),    intent(INOUT) :: self               !< Euler field.
  integer(I_P),           intent(IN)    :: Ni                 !< Space dimension (local image).
  integer(I_P),           intent(IN)    :: Ns                 !< Number of initial species.
  real(R_P),              intent(IN)    :: Dx                 !< Space step.
  character(*),           intent(IN)    :: BC_L               !< Left boundary condition type.
  character(*),           intent(IN)    :: BC_R               !< Right boundary condition type.
  real(R_P),              intent(IN)    :: initial_state(:,:) !< Initial state of primitive variables.
  real(R_P),              intent(IN)    :: cp0(:)             !< Initial specific heat, constant pressure.
  real(R_P),              intent(IN)    :: cv0(:)             !< Initial specific heat, constant volume.
  integer(I_P), optional, intent(IN)    :: ord                !< Space accuracy formal order.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%Ni = Ni
  self%ord = 1 ; if (present(ord)) self%ord = ord
  self%Ng = (self%ord + 1) / 2
#ifdef CAF
  if (allocated(self%local)) deallocate(self%local) ; allocate(self%local[*])
#else
  if (allocated(self%local)) deallocate(self%local) ; allocate(self%local)
#endif
  call self%local%init(Ni=Ni, Ns=Ns, Dx=Dx, BC_L=BC_L, BC_R=BC_R, initial_state=initial_state, cp0=cp0, cv0=cv0, ord=ord)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  pure subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy field.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D_caf), intent(INOUT) :: self !< Euler field.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%ord = 0
  self%Ni  = 0
  self%Ng  = 0
  if (allocated(self%local)) deallocate(self%local)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  pure function output(self) result(state)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Output the Euler field state (primitive variables).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D_caf), intent(IN)        :: self  !< Euler field.
  real(R_P), dimension(:,:), allocatable :: state !< Euler state vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  state = self%local%output()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction output

  subroutine synchronize(self, me, we)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Synchronize CAF images.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(euler_1D_caf), intent(INOUT) :: self          !< Euler field.
  integer(I_P),        intent(IN)    :: me            !< ID of this_image()
  integer(I_P),        intent(IN)    :: we            !< Number of CAF images used.
  real(R_P), allocatable             :: remote_U(:,:) !< Remote field.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef CAF
  ! if (we>1) then
  !   if (me==1) then
  !     sync images(me+1)
  !     remote_U = self%local[me+1]%output(iL=1, iU=self%Ng)
  !     call self%local%set_ghost_cells(U_R=remote_U)
  !   else if (me==we) then
  !     sync images(me-1)
  !     remote_U = self%local[me+1]%output(iL=self%Ni-self%Ng+1, iU=self%Ni)
  !     call self%local%set_ghost_cells(U_L=remote_U)
  !   else
  !     sync images([me-1, me+1])
  !     remote_U = self%local[me+1]%output(iL=self%Ni-self%Ng+1, iU=self%Ni)
  !     call self%local%set_ghost_cells(U_L=remote_U)
  !     remote_U = self%local[me+1]%output(iL=1, iU=self%Ng)
  !     call self%local%set_ghost_cells(U_R=remote_U)
  !   endif
  ! endif
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine synchronize
endmodule type_euler_1D_caf
