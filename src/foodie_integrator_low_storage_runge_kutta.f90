!< FOODIE integrator: provide an explicit class of low storage Runge-Kutta schemes, from 1st to 4th order accurate.
module foodie_integrator_low_storage_runge_kutta
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE integrator: provide an explicit class of low storage Runge-Kutta schemes, from 1st to 4th order accurate.
!<
!< The integrators provided have the low storage property allowing for an efficient use of the memory.
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
!<#### 1 stage, Explicit Forward Euler, 1st order
!< This scheme is TVD and reverts to Explicit Forward Euler, it being 1st order. It is not a real low storage scheme, this being
!< meaningless for a first order scheme. However it is added for safety reason.
!<
!< | ` Stage ` | ` A ` | ` B ` | ` C ` |
!< |-----------|-------|-------|-------|
!< | ` 1 `     | ` 0 ` | ` 1 ` | ` 0 ` |
!<
!<#### 5 stages, SSP, 4th order
!< This scheme is a low storage RK(5, 4), based on the *solution 3* proposed in [2].
!<
!< | ` Stage ` | ` A `                            | ` B `                            | ` C `                           |
!< |-----------|----------------------------------|----------------------------------|---------------------------------|
!< | ` 1 `     | `  0                           ` | ` 1432997174477/9575080441755  ` | ` 0                           ` |
!< | ` 2 `     | ` -567301805773 /1357537059087 ` | ` 5161836677717/13612068292357 ` | ` 1432997174477/9575080441755 ` |
!< | ` 3 `     | ` -2404267990393/2016746695238 ` | ` 1720146321549/2090206949498  ` | ` 2526269341429/6820363962896 ` |
!< | ` 4 `     | ` -3550918686646/2091501179385 ` | ` 3134564353537/4481467310338  ` | ` 2006345519317/3224310063776 ` |
!< | ` 5 `     | ` -1275806237668/842570457699  ` | ` 2277821191437/14882151754819 ` | ` 2802321613138/2924317926251 ` |
!<
!<
!<#### 6 stages, 4th order
!< This scheme is a low storage RK(6, 4), by [3].
!<
!< | ` Stage ` | ` A `               | ` B `              | ` C `              |
!< |-----------|---------------------|--------------------|--------------------|
!< | ` 1 `     | `  0              ` | ` 0.122000000000 ` | ` 0              ` |
!< | ` 2 `     | ` -0.691750960670 ` | ` 0.477263056358 ` | ` 0.122000000000 ` |
!< | ` 3 `     | ` -1.727127405211 ` | ` 0.381941220320 ` | ` 0.269115878630 ` |
!< | ` 4 `     | ` -0.694890150986 ` | ` 0.447757195744 ` | ` 0.447717183551 ` |
!< | ` 5 `     | ` -1.039942756197 ` | ` 0.498614246822 ` | ` 0.749979795490 ` |
!< | ` 6 `     | ` -1.531977447611 ` | ` 0.186648570846 ` | ` 0.898555413085 ` |
!<
!<#### 7 stages, 4th order
!< This scheme is a low storage RK(7, 4), by [3].
!<
!< | ` Stage ` | ` A `               | ` B `              | ` C `              |
!< |-----------|---------------------|--------------------|--------------------|
!< | ` 1 `     | `  0              ` | ` 0.117322146869 ` | ` 0              ` |
!< | ` 2 `     | ` -0.647900745934 ` | ` 0.503270262127 ` | ` 0.117322146869 ` |
!< | ` 3 `     | ` -2.704760863204 ` | ` 0.233663281658 ` | ` 0.294523230758 ` |
!< | ` 4 `     | ` -0.460080550118 ` | ` 0.283419634625 ` | ` 0.305658622131 ` |
!< | ` 5 `     | ` -0.500581787785 ` | ` 0.540367414023 ` | ` 0.582864148403 ` |
!< | ` 6 `     | ` -1.906532255913 ` | ` 0.371499414620 ` | ` 0.858664273599 ` |
!< | ` 7 `     | ` -1.450000000000 ` | ` 0.136670099385 ` | ` 0.868664273599 ` |
!<
!<#### 12 stages, 4th order
!< This scheme is a low storage RK(12, 4), by [4].
!<
!< | ` Stage ` | ` A `                   | ` B `                  | ` C `                  |
!< |-----------|-------------------------|------------------------|------------------------|
!< | ` 1  `    | `  0                  ` | ` 0.0650008435125904 ` | ` 0                  ` |
!< | ` 2  `    | ` -0.0923311242368072 ` | ` 0.0161459902249842 ` | ` 0.0650008435125904 ` |
!< | ` 3  `    | ` -0.9441056581158819 ` | ` 0.5758627178358159 ` | ` 0.0796560563081853 ` |
!< | ` 4  `    | ` -4.3271273247576394 ` | ` 0.1649758848361671 ` | ` 0.1620416710085376 ` |
!< | ` 5  `    | ` -2.1557771329026072 ` | ` 0.3934619494248182 ` | ` 0.2248877362907778 ` |
!< | ` 6  `    | ` -0.9770727190189062 ` | ` 0.0443509641602719 ` | ` 0.2952293985641261 ` |
!< | ` 7  `    | ` -0.7581835342571139 ` | ` 0.2074504268408778 ` | ` 0.3318332506149405 ` |
!< | ` 8  `    | ` -1.7977525470825499 ` | ` 0.6914247433015102 ` | ` 0.4094724050198658 ` |
!< | ` 9  `    | ` -2.6915667972700770 ` | ` 0.3766646883450449 ` | ` 0.6356954475753369 ` |
!< | ` 10 `    | ` -4.6466798960268143 ` | ` 0.0757190350155483 ` | ` 0.6806551557645497 ` |
!< | ` 11 `    | ` -0.1539613783825189 ` | ` 0.2027862031054088 ` | ` 0.7143773712418350 ` |
!< | ` 12 `    | ` -0.5943293901830616 ` | ` 0.2167029365631842 ` | ` 0.9032588871651854 ` |
!<
!<#### 13 stages, 4th order
!< This scheme is a low storage RK(13, 4), by [4].
!<
!< | ` Stage ` | ` A `                   | ` B `                  | ` C `                  |
!< |-----------|-------------------------|------------------------|------------------------|
!< | ` 1  `    | `  0                  ` | ` 0.0271990297818803 ` | ` 0                  ` |
!< | ` 2  `    | ` -0.6160178650170565 ` | ` 0.1772488819905108 ` | ` 0.0271990297818803 ` |
!< | ` 3  `    | ` -0.4449487060774118 ` | ` 0.0378528418949694 ` | ` 0.0952594339119365 ` |
!< | ` 4  `    | ` -1.0952033345276178 ` | ` 0.6086431830142991 ` | ` 0.1266450286591127 ` |
!< | ` 5  `    | ` -1.2256030785959187 ` | ` 0.2154313974316100 ` | ` 0.1825883045699772 ` |
!< | ` 6  `    | ` -0.2740182222332805 ` | ` 0.2066152563885843 ` | ` 0.3737511439063931 ` |
!< | ` 7  `    | ` -0.0411952089052647 ` | ` 0.0415864076069797 ` | ` 0.5301279418422206 ` |
!< | ` 8  `    | ` -0.1797084899153560 ` | ` 0.0219891884310925 ` | ` 0.5704177433952291 ` |
!< | ` 9  `    | ` -1.1771530652064288 ` | ` 0.9893081222650993 ` | ` 0.5885784947099155 ` |
!< | ` 10 `    | ` -0.4078831463120878 ` | ` 0.0063199019859826 ` | ` 0.6160769826246714 ` |
!< | ` 11 `    | ` -0.8295636426191777 ` | ` 0.3749640721105318 ` | ` 0.6223252334314046 ` |
!< | ` 12 `    | ` -4.7895970584252288 ` | ` 1.6080235151003195 ` | ` 0.6897593128753419 ` |
!< | ` 13 `    | ` -0.6606671432964504 ` | ` 0.0961209123818189 ` | ` 0.9126827615920843 ` |
!<
!<#### 14 stages, 4th order
!< This scheme is a low storage RK(14, 4), by [4].
!<
!< | ` Stage ` | ` A `                   | ` B `                  | ` C `                  |
!< |-----------|-------------------------|------------------------|------------------------|
!< | ` 1  `    | `  0                  ` | ` 0.0367762454319673 ` | ` 0                  ` |
!< | ` 2  `    | ` -0.7188012108672410 ` | ` 0.3136296607553959 ` | ` 0.0367762454319673 ` |
!< | ` 3  `    | ` -0.7785331173421570 ` | ` 0.1531848691869027 ` | ` 0.1249685262725025 ` |
!< | ` 4  `    | ` -0.0053282796654044 ` | ` 0.0030097086818182 ` | ` 0.2446177702277698 ` |
!< | ` 5  `    | ` -0.8552979934029281 ` | ` 0.3326293790646110 ` | ` 0.2476149531070420 ` |
!< | ` 6  `    | ` -3.9564138245774565 ` | ` 0.2440251405350864 ` | ` 0.2969311120382472 ` |
!< | ` 7  `    | ` -1.5780575380587385 ` | ` 0.3718879239592277 ` | ` 0.3978149645802642 ` |
!< | ` 8  `    | ` -2.0837094552574054 ` | ` 0.6204126221582444 ` | ` 0.5270854589440328 ` |
!< | ` 9  `    | ` -0.7483334182761610 ` | ` 0.1524043173028741 ` | ` 0.6981269994175695 ` |
!< | ` 10 `    | ` -0.7032861106563359 ` | ` 0.0760894927419266 ` | ` 0.8190890835352128 ` |
!< | ` 11 `    | `  0.0013917096117681 ` | ` 0.0077604214040978 ` | ` 0.8527059887098624 ` |
!< | ` 12 `    | ` -0.0932075369637460 ` | ` 0.0024647284755382 ` | ` 0.8604711817462826 ` |
!< | ` 13 `    | ` -0.9514200470875948 ` | ` 0.0780348340049386 ` | ` 0.8627060376969976 ` |
!< | ` 14 `    | ` -7.1151571693922548 ` | ` 5.5059777270269628 ` | ` 0.8734213127600976 ` |
!<
!<#### Bibliography
!<
!< [1] *Low-Storage Runge-Kutta Schemes*, J. H. Williamson, Journal of Computational Physics, vol. 35, 1980, pp. 48--56.
!<
!< [2] *Fourth-Order 2N-Storage Runge-Kutta Schemes*, Mark H. Carpenter, Christopher A. Kennedy, NASA Technical Memorandum 109112,
!< June 1994.
!<
!< [3] *High-accuracy large-step explicit Runge–Kutta (HALE-RK) schemes for computational aeroacoustics*, Vasanth Allampalli and
!< Ray Hixon and M. Nallasamy and Scott D. Sawyer, Journal of Computational Physics, vol. 228, 2009, pp. 3837--3850.
!<
!< [4] *Efficient low-storage Runge–Kutta schemes with optimized stability regions*, Jens Niegemann and Richard Diehl and Kurt
!< Busch, Journal of Computational Physics, vol. 231, 2012, pp. 364--372.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use, intrinsic:: ISO_FORTRAN_ENV, only: stderr => ERROR_UNIT
use foodie_kinds, only : R_P, I_P, I8P
use foodie_adt_integrand, only : integrand
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: ls_runge_kutta_integrator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: ls_runge_kutta_integrator
  !< FOODIE integrator: provide an explicit class of low storage Runge-Kutta schemes, from 1st to 4th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the RK coefficients) before used.
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
  subroutine init(self, stages)
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
  case(1) ! RK(1,1) Forward-Euler
    self%B(1) = 1._R_P
  case(5) ! LSRK(5,4)
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
  case(6) ! LSRK(6,4)
    self%A(1) =  0._R_P             ; self%B(1) = 0.122000000000_R_P ; self%C(1) = 0._R_P
    self%A(2) = -0.691750960670_R_P ; self%B(2) = 0.477263056358_R_P ; self%C(2) = 0.122000000000_R_P
    self%A(3) = -1.727127405211_R_P ; self%B(3) = 0.381941220320_R_P ; self%C(3) = 0.269115878630_R_P
    self%A(4) = -0.694890150986_R_P ; self%B(4) = 0.447757195744_R_P ; self%C(4) = 0.447717183551_R_P
    self%A(5) = -1.039942756197_R_P ; self%B(5) = 0.498614246822_R_P ; self%C(5) = 0.749979795490_R_P
    self%A(6) = -1.531977447611_R_P ; self%B(6) = 0.186648570846_R_P ; self%C(6) = 0.898555413085_R_P
  case(7) ! LSRK(7,4)
    self%A(1) =  0._R_P             ; self%B(1) = 0.117322146869_R_P ; self%C(1) = 0._R_P
    self%A(2) = -0.647900745934_R_P ; self%B(2) = 0.503270262127_R_P ; self%C(2) = 0.117322146869_R_P
    self%A(3) = -2.704760863204_R_P ; self%B(3) = 0.233663281658_R_P ; self%C(3) = 0.294523230758_R_P
    self%A(4) = -0.460080550118_R_P ; self%B(4) = 0.283419634625_R_P ; self%C(4) = 0.305658622131_R_P
    self%A(5) = -0.500581787785_R_P ; self%B(5) = 0.540367414023_R_P ; self%C(5) = 0.582864148403_R_P
    self%A(6) = -1.906532255913_R_P ; self%B(6) = 0.371499414620_R_P ; self%C(6) = 0.858664273599_R_P
    self%A(7) = -1.450000000000_R_P ; self%B(7) = 0.136670099385_R_P ; self%C(7) = 0.868664273599_R_P
  case(12) ! LSRK(12,4)
    self%A(1 ) =  0._R_P                 ; self%B(1 ) = 0.0650008435125904_R_P ; self%C(1 ) = 0._R_P
    self%A(2 ) = -0.0923311242368072_R_P ; self%B(2 ) = 0.0161459902249842_R_P ; self%C(2 ) = 0.0650008435125904_R_P
    self%A(3 ) = -0.9441056581158819_R_P ; self%B(3 ) = 0.5758627178358159_R_P ; self%C(3 ) = 0.0796560563081853_R_P
    self%A(4 ) = -4.3271273247576394_R_P ; self%B(4 ) = 0.1649758848361671_R_P ; self%C(4 ) = 0.1620416710085376_R_P
    self%A(5 ) = -2.1557771329026072_R_P ; self%B(5 ) = 0.3934619494248182_R_P ; self%C(5 ) = 0.2248877362907778_R_P
    self%A(6 ) = -0.9770727190189062_R_P ; self%B(6 ) = 0.0443509641602719_R_P ; self%C(6 ) = 0.2952293985641261_R_P
    self%A(7 ) = -0.7581835342571139_R_P ; self%B(7 ) = 0.2074504268408778_R_P ; self%C(7 ) = 0.3318332506149405_R_P
    self%A(8 ) = -1.7977525470825499_R_P ; self%B(8 ) = 0.6914247433015102_R_P ; self%C(8 ) = 0.4094724050198658_R_P
    self%A(9 ) = -2.6915667972700770_R_P ; self%B(9 ) = 0.3766646883450449_R_P ; self%C(9 ) = 0.6356954475753369_R_P
    self%A(10) = -4.6466798960268143_R_P ; self%B(10) = 0.0757190350155483_R_P ; self%C(10) = 0.6806551557645497_R_P
    self%A(11) = -0.1539613783825189_R_P ; self%B(11) = 0.2027862031054088_R_P ; self%C(11) = 0.7143773712418350_R_P
    self%A(12) = -0.5943293901830616_R_P ; self%B(12) = 0.2167029365631842_R_P ; self%C(12) = 0.9032588871651854_R_P
  case(13) ! LSRK(13,4)
    self%A(1 ) =  0._R_P                 ; self%B(1 ) = 0.0271990297818803_R_P ; self%C(1 ) = 0._R_P
    self%A(2 ) = -0.6160178650170565_R_P ; self%B(2 ) = 0.1772488819905108_R_P ; self%C(2 ) = 0.0271990297818803_R_P
    self%A(3 ) = -0.4449487060774118_R_P ; self%B(3 ) = 0.0378528418949694_R_P ; self%C(3 ) = 0.0952594339119365_R_P
    self%A(4 ) = -1.0952033345276178_R_P ; self%B(4 ) = 0.6086431830142991_R_P ; self%C(4 ) = 0.1266450286591127_R_P
    self%A(5 ) = -1.2256030785959187_R_P ; self%B(5 ) = 0.2154313974316100_R_P ; self%C(5 ) = 0.1825883045699772_R_P
    self%A(6 ) = -0.2740182222332805_R_P ; self%B(6 ) = 0.2066152563885843_R_P ; self%C(6 ) = 0.3737511439063931_R_P
    self%A(7 ) = -0.0411952089052647_R_P ; self%B(7 ) = 0.0415864076069797_R_P ; self%C(7 ) = 0.5301279418422206_R_P
    self%A(8 ) = -0.1797084899153560_R_P ; self%B(8 ) = 0.0219891884310925_R_P ; self%C(8 ) = 0.5704177433952291_R_P
    self%A(9 ) = -1.1771530652064288_R_P ; self%B(9 ) = 0.9893081222650993_R_P ; self%C(9 ) = 0.5885784947099155_R_P
    self%A(10) = -0.4078831463120878_R_P ; self%B(10) = 0.0063199019859826_R_P ; self%C(10) = 0.6160769826246714_R_P
    self%A(11) = -0.8295636426191777_R_P ; self%B(11) = 0.3749640721105318_R_P ; self%C(11) = 0.6223252334314046_R_P
    self%A(12) = -4.7895970584252288_R_P ; self%B(12) = 1.6080235151003195_R_P ; self%C(12) = 0.6897593128753419_R_P
    self%A(13) = -0.6606671432964504_R_P ; self%B(13) = 0.0961209123818189_R_P ; self%C(13) = 0.9126827615920843_R_P
  case(14) ! LSRK(14,4)
    self%A(1 ) =  0._R_P                 ; self%B(1 ) = 0.0367762454319673_R_P ; self%C(1 ) = 0._R_P
    self%A(2 ) = -0.7188012108672410_R_P ; self%B(2 ) = 0.3136296607553959_R_P ; self%C(2 ) = 0.0367762454319673_R_P
    self%A(3 ) = -0.7785331173421570_R_P ; self%B(3 ) = 0.1531848691869027_R_P ; self%C(3 ) = 0.1249685262725025_R_P
    self%A(4 ) = -0.0053282796654044_R_P ; self%B(4 ) = 0.0030097086818182_R_P ; self%C(4 ) = 0.2446177702277698_R_P
    self%A(5 ) = -0.8552979934029281_R_P ; self%B(5 ) = 0.3326293790646110_R_P ; self%C(5 ) = 0.2476149531070420_R_P
    self%A(6 ) = -3.9564138245774565_R_P ; self%B(6 ) = 0.2440251405350864_R_P ; self%C(6 ) = 0.2969311120382472_R_P
    self%A(7 ) = -1.5780575380587385_R_P ; self%B(7 ) = 0.3718879239592277_R_P ; self%C(7 ) = 0.3978149645802642_R_P
    self%A(8 ) = -2.0837094552574054_R_P ; self%B(8 ) = 0.6204126221582444_R_P ; self%C(8 ) = 0.5270854589440328_R_P
    self%A(9 ) = -0.7483334182761610_R_P ; self%B(9 ) = 0.1524043173028741_R_P ; self%C(9 ) = 0.6981269994175695_R_P
    self%A(10) = -0.7032861106563359_R_P ; self%B(10) = 0.0760894927419266_R_P ; self%C(10) = 0.8190890835352128_R_P
    self%A(11) =  0.0013917096117681_R_P ; self%B(11) = 0.0077604214040978_R_P ; self%C(11) = 0.8527059887098624_R_P
    self%A(12) = -0.0932075369637460_R_P ; self%B(12) = 0.0024647284755382_R_P ; self%C(12) = 0.8604711817462826_R_P
    self%A(13) = -0.9514200470875948_R_P ; self%B(13) = 0.0780348340049386_R_P ; self%C(13) = 0.8627060376969976_R_P
    self%A(14) = -7.1151571693922548_R_P ; self%B(14) = 5.5059777270269628_R_P ; self%C(14) = 0.8734213127600976_R_P
  case default
    write(stderr, '(A,I3,A)')' Error: ls_runge_kutta_integrator%init: ', stages, ' is an invalid number-of-stages!'
    write(stderr, '(A)')' Valid number of stages are:'
    write(stderr, '(A)')'   1  => LSRK(1, 1)'
    write(stderr, '(A)')'   5  => LSRK(5, 4)'
    write(stderr, '(A)')'   6  => LSRK(6, 4)'
    write(stderr, '(A)')'   7  => LSRK(7, 4)'
    write(stderr, '(A)')'   12 => LSRK(12, 4)'
    write(stderr, '(A)')'   13 => LSRK(13, 4)'
    write(stderr, '(A)')'   14 => LSRK(14, 4)'
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy the integrator.
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
  !< @note This method can be used **after** the integrator is created (i.e. the RK coefficients are initialized).
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
