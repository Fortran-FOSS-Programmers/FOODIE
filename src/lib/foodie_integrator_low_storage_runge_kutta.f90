!< FOODiE integrator: provide an explicit class of low storage Runge-Kutta schemes, from 1st to 4th order accurate.
module foodie_integrator_low_storage_runge_kutta
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODiE integrator: provide an explicit class of low storage Runge-Kutta schemes, from 1st to 4th order accurate.
!<
!< The class of integrators provided have the low storage property allowing for an efficient use of the memory.
!< Following Williamson approach [1], the LSRK(5,4)2N (solution 3) scheme of Carpenter et al. [2] is implemented.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the class of schemes implemented are written in the form:
!<
!<$$\begin{matrix}
!< K_1 = U^n \\
!< K_2 = 0 \\
!<\left.\begin{matrix}
!< K_2 = A_s K_2 + \Delta t R(t^n + C_s \Delta t, K_1) \\
!< K_1 = K_1 + B_s K_2
!<\end{matrix}\right\} s=1,2,...N_s\\
!<U^{n+1} = K_1
!<\end{matrix}$$
!<
!< where *Ns* is the number of stages used and \(K_1, K_2\) are the 2 registers used for stages computation.
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit thus \(A_1=C_1=0\). The coefficients \(A_s, B_s, C_s\) are given in the Williamson low storage table
!< form.
!<
!< The different schemes are selected accordingly to the number of stages used. Currently the following schemes are available:
!<
!<##### 1 stage, Explicit Forward Euler, 1st order
!< This scheme is TVD and reverts to Explicit Forward Euler, it being 1st order. It is not a real low storage scheme, this being
!< meaningless for a first order scheme. However it is added for safety reason.
!< $$A = \left[ 0 \right]$$
!< $$B = \left[ 1 \right]$$
!< $$C = \left[ 0 \right]$$
!<
!<##### 5 stages, SSP, 4th order
!< This scheme is a low storage RK(5, 4), based on the *solution 3* proposed in [2].
!<$$A = \left[ {\begin{matrix}
!< 0                                   \\
!<-\frac{567301805773 }{1357537059087} \\
!<-\frac{2404267990393}{2016746695238} \\
!<-\frac{3550918686646}{2091501179385} \\
!<-\frac{1275806237668}{842570457699 }
!<\end{matrix}} \right]$$
!<$$B = \left[ {\begin{matrix}
!<\frac{1432997174477}{9575080441755 } \\
!<\frac{5161836677717}{13612068292357} \\
!<\frac{1720146321549}{2090206949498 } \\
!<\frac{3134564353537}{4481467310338 } \\
!<\frac{2277821191437}{14882151754819}
!<\end{matrix}} \right]$$
!<$$C = \left[ {\begin{matrix}
!< 0                                  \\
!<\frac{1432997174477}{9575080441755} \\
!<\frac{2526269341429}{6820363962896} \\
!<\frac{2006345519317}{3224310063776} \\
!<\frac{2802321613138}{2924317926251}
!<\end{matrix}} \right]$$
!<
!<#### Bibliography
!<
!< [1] *Low-Storage Runge-Kutta Schemes*, J. H. Williamson, Journal of Computational Physics, vol. 35, 1980, pp. 48--56.
!<
!< [2] *Fourth-Order 2N-Storage Runge-Kutta Schemes*, Mark H. Carpenter, Christopher A. Kennedy, NASA Technical Memorandum 109112,
!< June 1994.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P, I8P
use type_integrand, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: ls_runge_kutta_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: ls_runge_kutta_integrator
  !< FOODiE integrator: provide an explicit class of low storage Runge-Kutta schemes, from 1st to 4th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the RK coeficients) before used.
  integer(I_P)           :: stages=0 !< Number of stages.
  real(R_P), allocatable :: A(:)     !< Low storage *A* coefficients.
  real(R_P), allocatable :: B(:)     !< Low storage *B* coefficients.
  real(R_P), allocatable :: C(:)     !< Low storage *C* coefficients.
  contains
    procedure, pass(self), public :: destroy   !< Destroy the integrator.
    procedure, pass(self), public :: init      !< Initialize (create) the integrator.
    procedure, pass(self), public :: integrate !< Integrate integrand field.
    final                         :: finalize  !< Finalize object.
endtype ls_runge_kutta_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  elemental subroutine init(self, stages)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual RK integrator: initialize the Butcher' low storage table coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(ls_runge_kutta_integrator), intent(INOUT) :: self   !< RK integrator.
  integer(I_P),                     intent(IN)    :: stages !< Number of stages used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (stages<1) return ! error print should be added
  self%stages = stages
  if (allocated(self%A)) deallocate(self%A) ; allocate(self%A(1:stages)) ; self%A = 0._R_P
  if (allocated(self%B)) deallocate(self%B) ; allocate(self%B(1:stages)) ; self%B = 0._R_P
  if (allocated(self%C)) deallocate(self%C) ; allocate(self%C(1:stages)) ; self%C = 0._R_P
  select case(stages)
  case(1)
    ! RK(1,1) Forward-Euler
    self%B(1) = 1._R_P
  case(5)
    ! LSRK(5,4)
    self%A(1) =  0._R_P
    self%A(2) = -real(567301805773_I8P,  kind=R_P) / real(1357537059087_I8P, kind=R_P)
    self%A(3) = -real(2404267990393_I8P, kind=R_P) / real(2016746695238_I8P, kind=R_P)
    self%A(4) = -real(3550918686646_I8P, kind=R_P) / real(2091501179385_I8P, kind=R_P)
    self%A(5) = -real(1275806237668_I8P, kind=R_P) / real(842570457699_I8P,  kind=R_P)

    self%B(1) =  real(1432997174477_I8P, kind=R_P) / real(9575080441755_I8P,  kind=R_P)
    self%B(2) =  real(5161836677717_I8P, kind=R_P) / real(13612068292357_I8P, kind=R_P)
    self%B(3) =  real(1720146321549_I8P, kind=R_P) / real(2090206949498_I8P,  kind=R_P)
    self%B(4) =  real(3134564353537_I8P, kind=R_P) / real(4481467310338_I8P,  kind=R_P)
    self%B(5) =  real(2277821191437_I8P, kind=R_P) / real(14882151754819_I8P, kind=R_P)

    self%C(1) =  0._R_P
    self%C(2) =  real(1432997174477_I8P, kind=R_P) / real(9575080441755_I8P, kind=R_P)
    self%C(3) =  real(2526269341429_I8P, kind=R_P) / real(6820363962896_I8P, kind=R_P)
    self%C(4) =  real(2006345519317_I8P, kind=R_P) / real(3224310063776_I8P, kind=R_P)
    self%C(5) =  real(2802321613138_I8P, kind=R_P) / real(2924317926251_I8P, kind=R_P)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy the integrator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(ls_runge_kutta_integrator), intent(INOUT) :: self !< Integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%stages = 0
  if (allocated(self%A)) deallocate(self%A)
  if (allocated(self%B)) deallocate(self%B)
  if (allocated(self%C)) deallocate(self%C)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine integrate(self, U, stage, Dt, t)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with explicit low storage Runge-Kutta scheme.
  !<
  !< @note This method can be used **after** the integrator is created (i.e. the RK coeficients are initialized).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(ls_runge_kutta_integrator), intent(IN)    :: self       !< Actual RK integrator.
  class(integrand),                 intent(INOUT) :: U          !< Field to be integrated.
  class(integrand),                 intent(INOUT) :: stage(1:2) !< Runge-Kutta registers [1:2].
  real(R_P),                        intent(IN)    :: Dt         !< Time step.
  real(R_P),                        intent(IN)    :: t          !< Time.
  integer(I_P)                                    :: s          !< First stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(stage)
  class is(integrand)
    ! computing stages
    stage(1) = U
    stage(2) = U*0._R_P
    do s=1, self%stages
      stage(2) = stage(2) * self%A(s) + stage(1)%t(t=t + self%C(s) * Dt) * Dt
      stage(1) = stage(1) + stage(2) * self%B(s)
    enddo
    U = stage(1)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate

  ! private methods
  elemental subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Finalize object.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(ls_runge_kutta_integrator), intent(INOUT) :: self !< Integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule foodie_integrator_low_storage_runge_kutta
