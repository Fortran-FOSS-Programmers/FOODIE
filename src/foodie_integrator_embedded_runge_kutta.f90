!< FOODIE integrator: provide an explicit class of embedded Runge-Kutta schemes, from 1st to 4th order accurate.
module foodie_integrator_emd_runge_kutta
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE integrator: provide an explicit class of embedded Runge-Kutta schemes, from 1st to 4th order accurate.
!<
!< The integrators provided have the embedded pairs property allowing for automatic step size control.
!< The schemes are explicit and defined through the extended Butcher's table syntax, see[1] .
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the class of schemes implemented are written in the form:
!<
!< $$ U_p^{n+1} = U^n +\Delta t \sum_{s=1}^{Ns}\beta_p^s K^s $$
!< $$ U_{p+1}^{n+1} = U^n +\Delta t \sum_{s=1}^{Ns}\beta_{p+1}^s K^s $$
!<
!< *p* is the lower accuracy order scheme and *p+1* is the higher one; Ns is the number of stages used and \(K^s\) is
!< the \(s^{th}\) stage computed as:
!<
!< $$ K^s = R\left( t^n+\gamma^s \Delta t, U^n +\Delta t \sum_{i=1}^{s-1}\alpha^{s,i} K^i \right) $$
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit thus the above summation is up to \(s-1\). The coefficients \(\beta\), \(\alpha\) and \(\gamma\) are
!< given in the extended Butcher table form:
!<
!<```
!<  gamma^1    | alpha^{1,1}       alpha^{1,2}       ...        alpha^{1,Ns}
!<  gamma^2    | alpha^{2,1}       alpha^{2,2}       ...        alpha^{2,Ns}
!<  .          | .                 .                 .          .
!<  .          | .                 .                  .         .
!<  .          | .                 .                   .        .
!<  gamma^{Ns} | alpha^{Ns,1}      alpha^{Ns,2}      ...        alpha^{Ns,Ns}
!< ------------|-------------------------------------------------------------
!<             | beta_p^1          beta_p^2          ...        beta_p^{Ns}
!<             | beta_{p+1}^1      beta_{p+1}^2      ...        beta_{p+1}^{Ns}
!<```
!<
!< Because only explicit schemes are considered the Butcher table reduces to diagonal matrix:
!<
!<```
!<  gamma^1    | 0                 0                 ...        0
!<  gamma^2    | alpha^{2,1}       0                 ...        0
!<  .          | .                 .                 .          .
!<  .          | .                 .                  .         .
!<  .          | .                 .                   .        .
!<  gamma^{Ns} | alpha^{Ns,1}      alpha^{Ns,2}      ...        0
!< ------------|-------------------------------------------------------------
!<             | beta_p^1          beta_p^2          ...        beta_p^{Ns}
!<             | beta_{p+1}^1      beta_{p+1}^2      ...        beta_{p+1}^{Ns}
!<```
!<
!< Moreover the following relation always holds:
!< \( \gamma^s = \sum_{i=1}^{Ns}\alpha^{s,i} \)
!<
!< The different schemes are selected accordingly to the number of stages used. Currently the following schemes are available:
!<
!<##### 7 stages, 4th order
!< This scheme is Dormand and Prince scheme, see [1].
!<```
!<  0    | 0             0            0             0            0              0          0
!<  1/5  | 1/5           0            0             0            0              0          0
!<  3/10 | 3/40          9/40         0             0            0              0          0
!<  4/5  | 44/45        -56/15        32/9          0            0              0          0
!<  8/9  | 19372/6561   -25360/2187   64448/6561   -212/729      0              0          0
!<  1    | 9017/3168    -355/33       46732/5247    49/176      -5103/18656     0          0
!<  1    | 35/384        0            500/1113      125/192     -2187/6784      11/84      0
!< --------------------------------------------------------------------------------------------
!<       | 35/384        0            500/1113      125/192     -2187/6784      11/84      0
!<       | 5179/57600    0            7571/16695    393/640     -92097/339200   187/2100   1/40
!<```
!<
!<#### Bibliography
!<
!< [1] *A family of embedded Runge-Kutta formulae*, Dormand, J. R.; Prince, P. J. (1980), , Journal of Computational and
!< Applied Mathematics 6 (1): 19--26, doi:10.1016/0771-050X(80)90013-3
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use foodie_kinds, only : R_P, I_P
use foodie_adt_integrand, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: emd_runge_kutta_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: emd_runge_kutta_integrator
  !< FOODIE integrator: provide an explicit class of TVD or SSP Runge-Kutta schemes, from 1st to 4th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the RK coefficients) before used.
  real(R_P)              :: tolerance=0._R_P !< Tolerance on the local truncation error.
  real(R_P)              :: pp1_inv=0._R_P   !< 1/(p+1) where p is the accuracy order of the lower accurate scheme of the pair.
  integer(I_P)           :: stages=0         !< Number of stages.
  real(R_P), allocatable :: alph(:,:)        !< \(\alpha\) Butcher's coefficients.
  real(R_P), allocatable :: beta(:,:)        !< \(\beta\) Butcher's coefficients.
  real(R_P), allocatable :: gamm(:)          !< \(\gamma\) Butcher's coefficients.
  contains
    procedure, pass(self), public  :: destroy   !< Destroy the integrator.
    procedure, pass(self), public  :: init      !< Initialize (create) the integrator.
    procedure, pass(self), public  :: integrate !< Integrate integrand field.
    procedure, pass(self), private :: new_Dt    !< Compute new estimation of the time step Dt.
    final                          :: finalize  !< Finalize object.
endtype emd_runge_kutta_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  elemental subroutine init(self, stages, tolerance)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual RK integrator: initialize the Butcher' table coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(emd_runge_kutta_integrator), intent(INOUT) :: self      !< RK integrator.
  integer(I_P),                      intent(IN)    :: stages    !< Number of stages used.
  real(R_P), optional,               intent(IN)    :: tolerance !< Tolerance on the local truncation error (default 0.01).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (stages<1) return ! error print should be added
  if (present(tolerance)) then
    self%tolerance = tolerance
  else
    self%tolerance = 0.01_R_P
  endif
  self%stages = stages
  if (allocated(self%beta)) deallocate(self%beta) ; allocate(self%beta(1:stages, 1:2     )) ; self%beta = 0._R_P
  if (allocated(self%alph)) deallocate(self%alph) ; allocate(self%alph(1:stages, 1:stages)) ; self%alph = 0._R_P
  if (allocated(self%gamm)) deallocate(self%gamm) ; allocate(self%gamm(          1:stages)) ; self%gamm = 0._R_P
  select case(stages)
  case(7)
    ! DPRK(7,4)
    self%pp1_inv = 1._R_P/(4._R_P + 1._R_P)

    self%beta(1, 1) =   35._R_P/384._R_P    ; self%beta(1, 2) =  5179._R_P/57600._R_P
    self%beta(2, 1) =   0._R_P              ; self%beta(2, 2) =  0._R_P
    self%beta(3, 1) =   500._R_P/1113._R_P  ; self%beta(3, 2) =  7571._R_P/16695._R_P
    self%beta(4, 1) =   125._R_P/192._R_P   ; self%beta(4, 2) =  393._R_P/640._R_P
    self%beta(5, 1) =  -2187._R_P/6784._R_P ; self%beta(5, 2) = -92097._R_P/339200._R_P
    self%beta(6, 1) =   11._R_P/84._R_P     ; self%beta(6, 2) =  187._R_P/2100._R_P
    self%beta(7, 1) =   0._R_P              ; self%beta(7, 2) =  1._R_P/40._R_P

    self%alph(2, 1)=1._R_P/5._R_P
    self%alph(3, 1)=3._R_P/40._R_P      ;self%alph(3, 2)= 9._R_P/40._R_P
    self%alph(4, 1)=44._R_P/45._R_P     ;self%alph(4, 2)=-56._R_P/15._R_P     ;self%alph(4, 3)=32._R_P/9._R_P
    self%alph(5, 1)=19372._R_P/6561._R_P;self%alph(5, 2)=-25360._R_P/2187._R_P;self%alph(5, 3)=64448._R_P/6561._R_P
    self%alph(6, 1)=9017._R_P/3168._R_P ;self%alph(6, 2)=-355._R_P/33._R_P    ;self%alph(6, 3)=46732._R_P/5247._R_P
    self%alph(7, 1)=35._R_P/384._R_P    ;self%alph(7, 2)= 0._R_P              ;self%alph(7, 3)=500._R_P/1113._R_P

    self%alph(5, 4)=-212._R_P/729._R_P
    self%alph(6, 4)= 49._R_P/176._R_P ;self%alph(6, 5)=-5103._R_P/18656._R_P
    self%alph(7, 4)= 125._R_P/192._R_P;self%alph(7, 5)=-2187._R_P/6784._R_P ;self%alph(7, 6)=11._R_P/84._R_P

    self%gamm(2) = 1._R_P/5._R_P
    self%gamm(3) = 3._R_P/10._R_P
    self%gamm(4) = 4._R_P/5._R_P
    self%gamm(5) = 8._R_P/9._R_P
    self%gamm(6) = 1._R_P
    self%gamm(7) = 1._R_P
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy the integrator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(emd_runge_kutta_integrator), intent(INOUT) :: self !< Integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%tolerance = 0._R_P
  self%stages = 0
  if (allocated(self%alph)) deallocate(self%alph)
  if (allocated(self%beta)) deallocate(self%beta)
  if (allocated(self%gamm)) deallocate(self%gamm)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine integrate(self, U, stage, Dt, t)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with explicit embedded Runge-Kutta scheme.
  !<
  !< The time steps is adaptively resized using the local truncation error estimation by means of the embedded pairs of RK schemes.
  !<
  !< @note This method can be used **after** the integrator is created (i.e. the RK coefficients are initialized).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(emd_runge_kutta_integrator), intent(IN)    :: self      !< Actual RK integrator.
  class(integrand),                  intent(INOUT) :: U         !< Field to be integrated.
  class(integrand),                  intent(INOUT) :: stage(1:) !< Runge-Kutta stages [1:stages].
  real(R_P),                         intent(INOUT) :: Dt        !< Time step.
  real(R_P),                         intent(IN)    :: t         !< Time.
  class(integrand), allocatable                    :: U1        !< First U evaluation.
  class(integrand), allocatable                    :: U2        !< Second U evaluation.
  real(R_P)                                        :: error     !< Local truncation error estimation.
  integer(I_P)                                     :: s         !< First stages counter.
  integer(I_P)                                     :: ss        !< Second stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(U1, source=U)
  allocate(U2, source=U)
  error = 1e6
  do while(error>self%tolerance)
    ! compute stages
    do s=1, self%stages
      stage(s) = U
      do ss=1, s - 1
        stage(s) = stage(s) + stage(ss) * (Dt * self%alph(s, ss))
      enddo
      stage(s) = stage(s)%t(t=t + self%gamm(s) * Dt)
    enddo
    ! compute new time step
    U1 = U
    U2 = U
    do s=1, self%stages
      U1 = U1 +  stage(s) * (Dt * self%beta(s, 1))
      U2 = U2 +  stage(s) * (Dt * self%beta(s, 2))
    enddo
    error = U2.lterror.U1
    call self%new_Dt(error=error, Dt=Dt)
  enddo
  U = U1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate

  ! private methods
  elemental subroutine new_Dt(self, error, Dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute new estimation of the time step Dt.
  !<
  !< The formula employed is:
  !<
  !< $$ Dt_{new} = 0.9 Dt_{old} \left( \frac{tolerance}{error} \right)^{\frac{1}{p+1}} $$
  !<
  !< @note 0.9 is a safety factor.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(emd_runge_kutta_integrator), intent(IN)    :: self  !< Integrator.
  real(R_P),                         intent(IN)    :: error !< Local truncation error estimation.
  real(R_P),                         intent(INOUT) :: Dt    !< Time step.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (error>self%tolerance) Dt = 0.9_R_P * Dt * (self%tolerance/error)**self%pp1_inv
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine new_Dt

  elemental subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Finalize object.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(emd_runge_kutta_integrator), intent(INOUT) :: self !< Integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule foodie_integrator_emd_runge_kutta
