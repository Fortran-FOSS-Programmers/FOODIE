!< FOODiE integrator: provide an explicit class of TVD or SSP Runge-Kutta schemes, from 1st to 4th order accutate.
module foodie_integrator_tvd_runge_kutta
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODiE integrator: provide an explicit class of TVD or SSP Runge-Kutta schemes, from 1st to 4th order accutate.
!<
!< The class of integrators provided have the Total Variation Diminishing (TVD) property or the Strong Stability Preserving (SSP)
!< one. The schemes are explicit and defined through the Butcher's table syntax (see Butcher, J.C., *Coefficients for the study of
!< Runge-Kutta integration processes*, J. Austral. Math. Soc., Vol. 3, pages: 185--201, 1963).
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the class of schemes implemented are written in the form:
!<
!< $$ U^{n+1} = U^n +\Delta t \sum_{s=1}^{Ns}\beta^s K^s $$
!<
!< where Ns is the number of stages used and \(K^s\) is the \(s^{th}\) stage computed as:
!<
!< $$ K^s = R\left( t^n+\gamma^s \Delta t, U^n +\Delta t \sum_{i=1}^{s-1}\alpha^{s,i} K^i \right) $$
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit thus the above summation is up to \(s-1\). The coefficients \(\beta\), \(\alpha\) and \(\gamma\) are
!< given in the Butcher table form:
!<
!< |                 |                   |                   |            |                    |
!< |-----------------|-------------------|-------------------|------------|--------------------|
!< | \(\gamma^1\)    | \(\alpha^{1,1}\)  | \(\alpha^{1,2}\)  | \(\cdots\) | \(\alpha^{1,Ns}\)  |
!< | \(\gamma^2\)    | \(\alpha^{2,1}\)  | \(\alpha^{2,2}\)  | \(\cdots\) | \(\alpha^{2,Ns}\)  |
!< | \(\vdots\)      | \(\vdots\)        | \(\vdots\)        | \(\ddots\) | \(\vdots\)         |
!< | \(\gamma^{Ns}\) | \(\alpha^{Ns,1}\) | \(\alpha^{Ns,2}\) | \(\cdots\) | \(\alpha^{Ns,Ns}\) |
!< |                 | \(\beta^1\)       | \(\beta^2\)       | \(\cdots\) | \(\beta^{Ns}\)     |
!<
!< Because only explicit schemes are considered the Butcher table reduces to diagonal matrix:
!<
!< |                 |                   |                   |            |                    |
!< |-----------------|-------------------|-------------------|------------|--------------------|
!< | \(\gamma^1\)    | \(      0     \)  | \(      0     \)  | \(\cdots\) | \(        0    \)  |
!< | \(\gamma^2\)    | \(\alpha^{2,1}\)  | \(\alpha^{2,2}\)  | \(\cdots\) | \(\alpha^{2,Ns}\)  |
!< | \(\vdots\)      | \(\vdots\)        | \(\vdots\)        | \(\ddots\) | \(\vdots\)         |
!< | \(\gamma^{Ns}\) | \(\alpha^{Ns,1}\) | \(\alpha^{Ns,2}\) | \(\cdots\) | \(        0     \) |
!< |                 | \(\beta^1\)       | \(\beta^2\)       | \(\cdots\) | \(\beta^{Ns}\)     |
!<
!< Moreover the following relation always holds:
!< \( \gamma^s = \sum_{i=1}^{Ns}\alpha^{s,i} \)
!<
!< The different schemes are selected accordingly to the number of stages used. Currently the following schemes are available:
!<
!<##### 1 stage, Explicit Forward Euler, 1st order
!< This scheme is TVD and reverts to Explicit Forward Euler, it being 1st order.
!< $$\beta = \left[1\right]$$
!< $$\alpha = \left[ {\begin{array}{*{20}{c}} 0&0 \end{array}} \right]$$
!< $$\gamma = \left[0\right]$$

!<##### 2 stages, SSP, 2nd order
!< This scheme is an optmial SSP(2, 2) without low-storage algorithm. See *High Order Strong Stability Preserving Time
!< Discretizations*, Gottlieb, S., Ketcheson, D. I., Ketcheson, D., Journal of Scientific Computing, vol. 38, N. 3, 2009,
!< pp. 251-289.
!< $$\beta = \left[ {\begin{array}{*{20}{c}} \frac{1}{2}&\frac{1}{2} \end{array}} \right]$$
!< $$\alpha = \left[ {\begin{array}{*{20}{c}} 0&0\\ 1&0 \end{array}} \right]$$
!< $$\gamma = \left[ {\begin{array}{*{20}{c}} 0 \\ 1 \end{array}} \right]$$

!<##### 3 stages, SSP, 3rd order
!< This scheme is an optmial SSP(3, 3) without low-storage algorithm. See *High Order Strong Stability Preserving Time
!< Discretizations*, Gottlieb, S., Ketcheson, D. I., Ketcheson, D., Journal of Scientific Computing, vol. 38, N. 3, 2009,
!< pp. 251-289.
!< $$\beta = \left[ {\begin{array}{*{20}{c}} \frac{1}{6}&\frac{1}{6}&\frac{1}{3}  \end{array}} \right]$$
!< $$\alpha = \left[ {\begin{array}{*{20}{c}} 0&0&0\\ 1&0&0\\ \frac{1}{4}&\frac{1}{4}&0 \end{array}} \right]$$
!< $$\gamma = \left[ {\begin{array}{*{20}{c}} 0 \\ 1 \\ \frac{1}{2} \end{array}} \right]$$

!<##### 5 stages, SSP, 4th order
!< This scheme is an optmial SSP(5, 4) without low-storage algorithm. See *High Order Strong Stability Preserving Time
!< Discretizations*, Gottlieb, S., Ketcheson, D. I., Ketcheson, D., Journal of Scientific Computing, vol. 38, N. 3, 2009,
!< pp. 251-289.
!< $$\beta = \left[ {\begin{array}{*{20}{c}} 0.14681187618661&0.24848290924556&0.10425883036650&0.27443890091960&
!< 0.22600748319395 \end{array}} \right]$$
!< $$\alpha = \left[ {\begin{array}{*{20}{c}}
!< 0&0&0&0&0 \\ 0.39175222700392&0&0&0&0 \\ 0.21766909633821&0.36841059262959&0&0&0 \\ 0.08269208670950&
!< 0.13995850206999&0.25189177424738&0&0 \\ 0.06796628370320&0.11503469844438&0.20703489864929&0.54497475021237&0
!< \end{array}} \right]$$
!< $$\gamma = \left[ {\begin{array}{*{20}{c}} 0\\0.39175222700392\\0.58607968896780\\0.47454236302687\\0.93501063100924
!<           \end{array}} \right]$$
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : R_P, I_P
use type_integrand, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: tvd_runge_kutta_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: tvd_runge_kutta_integrator
  !< FOODiE integrator: provide an explicit class of TVD or SSP Runge-Kutta schemes, from 1st to 4th order accutate.
  !<
  !< @note The integrator must be created or initialized (initialize the RK coeficients) before used.
  real(R_P), allocatable:: alph(:,:) !< \(\alpha\) Butcher's coefficients.
  real(R_P), allocatable:: beta(:)   !< \(\beta\) Butcher's coefficients.
  real(R_P), allocatable:: gamm(:)   !< \(\gamma\) Butcher's coefficients.
  contains
    procedure, pass(self), public :: destroy   !< Destroy the integrator.
    procedure, pass(self), public :: init      !< Initialize (create) the integrator.
    procedure, pass(self), public :: integrate !< Integrate integrand field.
    final                         :: finalize  !< Finalize object.
endtype tvd_runge_kutta_integrator
interface tvd_runge_kutta_integrator
  !< Overload tvd_runge_kutta_integrator name with the constructor function *create*.
  module procedure create
endinterface tvd_runge_kutta_integrator
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function create(stages) result(rk)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual RK integrators: initialize the Butcher' table coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)         :: stages !< Number of stages used.
  type(tvd_runge_kutta_integrator) :: rk     !< Actual RK integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call rk%init(stages=stages)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction create

  ! public methods
  elemental subroutine init(self, stages)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the actual RK integrators: initialize the Butcher' table coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(tvd_runge_kutta_integrator), intent(INOUT) :: self   !< RK integrator.
  integer(I_P),                      intent(IN)    :: stages !< Number of stages used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%beta)) deallocate(self%beta) ; allocate(self%beta(1:stages          )) ; self%beta = 0._R_P
  if (allocated(self%alph)) deallocate(self%alph) ; allocate(self%alph(1:stages, 1:stages)) ; self%alph = 0._R_P
  if (allocated(self%gamm)) deallocate(self%gamm) ; allocate(self%gamm(          1:stages)) ; self%gamm = 0._R_P
  select case(stages)
  case(1)
    ! RK(1,1) Forward-Euler
    self%beta(1) = 1._R_P
  case(2)
    ! SSPRK(2,2)
    self%beta(1) = 0.5_R_P
    self%beta(2) = 0.5_R_P

    self%alph(2, 1) = 1._R_P

    self%gamm(2) = 1._R_P
  case(3)
    ! SSPRK(3,3)
    self%beta(1) = 1._R_P/6._R_P
    self%beta(2) = 1._R_P/6._R_P
    self%beta(3) = 2._R_P/3._R_P

    self%alph(2, 1) = 1._R_P
    self%alph(3, 1) = 0.25_R_P ; self%alph(3, 2) = 0.25_R_P

    self%gamm(2) = 1._R_P
    self%gamm(3) = 0.5_R_P
  case(5)
    ! SSPRK(5,4)
    self%beta(1) = 0.14681187618661_R_P
    self%beta(2) = 0.24848290924556_R_P
    self%beta(3) = 0.10425883036650_R_P
    self%beta(4) = 0.27443890091960_R_P
    self%beta(5) = 0.22600748319395_R_P

    self%alph(2, 1)=0.39175222700392_R_P
    self%alph(3, 1)=0.21766909633821_R_P;self%alph(3, 2)=0.36841059262959_R_P
    self%alph(4, 1)=0.08269208670950_R_P;self%alph(4, 2)=0.13995850206999_R_P;self%alph(4, 3)=0.25189177424738_R_P
    self%alph(5, 1)=0.06796628370320_R_P;self%alph(5, 2)=0.11503469844438_R_P;self%alph(5, 3)=0.20703489864929_R_P
      self%alph(5, 4)=0.54497475021237_R_P

    self%gamm(2) = 0.39175222700392_R_P
    self%gamm(3) = 0.58607968896780_R_P
    self%gamm(4) = 0.47454236302687_R_P
    self%gamm(5) = 0.93501063100924_R_P
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy the integrator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(tvd_runge_kutta_integrator), intent(INOUT) :: self !< Integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%alph)) deallocate(self%alph)
  if (allocated(self%beta)) deallocate(self%beta)
  if (allocated(self%gamm)) deallocate(self%gamm)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine integrate(self, field, stage, Dt)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Integrate field with explicit TVD (or SSP) Runge-Kutta scheme.
  !<
  !< @note This method can be used **after** the integrator is created (i.e. the RK coeficients are initialized).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(tvd_runge_kutta_integrator), intent(IN)    :: self      !< Actual RK integrator.
  class(integrand),                  intent(INOUT) :: field     !< Field to be integrated.
  class(integrand),                  intent(INOUT) :: stage(1:) !< Runge-Kutta stages [1:stages].
  real(R_P),                         intent(IN)    :: Dt        !< Time step.
  integer(I_P)                                     :: s         !< First stages counter.
  integer(I_P)                                     :: ss        !< Second stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(stage)
  class is(integrand)
    ! computing stages
    do s=1, size(stage)
      stage(s) = field
      do ss=1, s - 1
        stage(s) = stage(s) + stage(ss) * (Dt * self%alph(s, ss))
      enddo
      stage(s) = stage(s)%t()
    enddo
    ! computing new time step
    do s=1, size(stage)
      field = field +  stage(s) * (Dt * self%beta(s))
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine integrate

  ! private methods
  elemental subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Finalize object.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(tvd_runge_kutta_integrator), intent(INOUT) :: self !< Integrator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule foodie_integrator_tvd_runge_kutta
