!< FOODIE integrator: provide an explicit class of embedded Runge-Kutta schemes, from 2nd to 6th order accurate.
module foodie_integrator_emd_runge_kutta
!-----------------------------------------------------------------------------------------------------------------------------------
!< FOODIE integrator: provide an explicit class of embedded Runge-Kutta schemes, from 2nd to 6th order accurate.
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
!<             | beta_{p+1}^1      beta_{p+1}^2      ...        beta_{p+1}^{Ns}
!<             | beta_p^1          beta_p^2          ...        beta_p^{Ns}
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
!<             | beta_{p+1}^1      beta_{p+1}^2      ...        beta_{p+1}^{Ns}
!<             | beta_p^1          beta_p^2          ...        beta_p^{Ns}
!<```
!<
!< Moreover the following relation always holds:
!< \( \gamma^s = \sum_{i=1}^{Ns}\alpha^{s,i} \)
!<
!< The different schemes are selected accordingly to the number of stages used. Currently the following schemes are available:
!<
!<##### 2 stages, 2th order
!< This scheme is due to Heun-Euler.
!<```
!<  0  | 0
!<  1  | 1     0
!< ----------------
!<     | 1/2   1/2
!<     | 1      0
!<```
!<
!<##### 6 stages, 5th order
!< This scheme is due to Cash and Karp, see [3].
!<```
!<  0    | 0
!<  1/5	 | 1/5
!<  3/10 | 3/40	         9/40
!<  3/5	 | 3/10	         -9/10	      6/5
!<  1	   | -11/54	       5/2	        -70/27	    35/27
!<  7/8	 | 1631/55296    175/512      575/13824   44275/110592     253/4096     0
!< ----------------------------------------------------------------------------------------
!<       | 37/378        0           250/621      125/594          0            512/1771
!<       | 2825/27648    0           18575/48384  13525/55296      277/14336    1/4
!<```
!<
!<##### 7 stages, 4th order
!< This scheme is due to Dormand and Prince, see [1].
!<```
!<  0    | 0
!<  1/5  | 1/5
!<  3/10 | 3/40          9/40
!<  4/5  | 44/45        -56/15        32/9
!<  8/9  | 19372/6561   -25360/2187   64448/6561   -212/729
!<  1    | 9017/3168    -355/33       46732/5247    49/176      -5103/18656
!<  1    | 35/384        0            500/1113      125/192     -2187/6784      11/84      0
!< --------------------------------------------------------------------------------------------
!<       | 5179/57600    0            7571/16695    393/640     -92097/339200   187/2100   1/40
!<       | 35/384        0            500/1113      125/192     -2187/6784      11/84      0
!<```
!<
!<##### 9 stages, 6th order
!< This scheme is due to Calvo et al., see [2].
!<```
!<  0                 | 0
!<  2/15              | 2/15
!<  1/5               | 1/20                  3/20
!<  3/10              | 3/40                  0                      9/40
!<  14/25             | 86727015/196851553    -60129073/52624712     957436434/1378352377    83886832/147842441
!<  19/25             | -86860849/45628967    111022885/25716487     108046682/101167669     -141756746/36005461
!<  35226607/35688279 | 77759591/16096467     -49252809/6452555      -381680111/51572984     879269579/66788831
!<  1                 | 237564263/39280295    -100523239/10677940    -265574846/27330247     317978411/18988713
!<  1                 | 17572349/289262523    0                      57513011/201864250      15587306/354501571
!< --------------------------------------------------------------------------------------------------------------
!<                    | 17572349/289262523    0                      57513011/201864250      15587306/354501571
!<                    | 15231665/510830334    0                      59452991/116050448      -28398517/122437738
!< ...continued...
!<  0                 |
!<  2/15              |
!<  1/5               |
!<  3/10              |
!<  14/25             |
!<  19/25             | 73139862/60170633
!<  35226607/35688279 | -90453121/33722162     111179552/157155827
!<  1                 | -124494385/35453627    86822444/100138635     -12873523/724232625
!<  1                 | 71783021/234982865     29672000/180480167     65567621/127060952     -79074570/210557597    0
!< -----------------------------------------------------------------------------------------------------------------------
!<                    | 71783021/234982865     29672000/180480167     65567621/127060952     -79074570/210557597    0
!<                    | 56673824/137010559     68003849/426673583     7097631/37564021       -71226429/583093742    1/20
!<```
!<
!<##### 17 stages, 10th order
!< This scheme is due to Feagin, see [4].
!<```
!<  0                        |  0
!<  0.1                      |  0.1
!<  0.539357840802981787532  | -0.915176561375291440520  1.454534402178273228052
!<  0.809036761204472681298  |  0.202259190301118170324  0                        0.606777570903354510974
!<  0.309036761204472681298  |  0.184024714708643575149  0                        0.197966831227192369068
!<  0.981074190219795268254  |  0.087900734020668133731  0                        0
!<  0.833333333333333333333  |  0.085970050490246030218  0                        0
!<  0.354017365856802376329  |  0.120930449125333720660  0                        0
!<  0.882527661964732346425  |  0.110854379580391483508  0                        0
!<  0.642615758240322548157  |  0.112054414752879004829  0                        0
!<  0.357384241759677451842  |  0.113976783964185986138  0                        0
!<  0.117472338035267653574  |  0.079831452828019604635  0                        0
!<  0.833333333333333333333  |  0.985115610164857280120  0                        0
!<  0.309036761204472681298  |  0.895080295771632891049  0                        0.197966831227192369068
!<  0.539357840802981787532  | -0.915176561375291440520  1.454534402178273228052  0
!<  0.1                      |  0.1                      0                       -0.157178665799771163367
!<  1                        |  0.181781300700095283888  0.675                    0.342758159847189839942
!< ------------------------------------------------------------------------------------------------------
!<                           |  0.033333333333333333333  0.025                    0.033333333333333333333
!< ...continued...
!<  0                        |
!<  0.1                      |
!<  0.539357840802981787532  |
!<  0.809036761204472681298  |
!<  0.309036761204472681298  | -0.072954784731363262918
!<  0.981074190219795268254  |  0.410459702520260645318  0.482713753678866489204
!<  0.833333333333333333333  |  0.330885963040722183948  0.489662957309450192844 -0.073185637507085073678
!<  0.354017365856802376329  |  0                        0.260124675758295622809  0.032540262154909133015
!<  0.882527661964732346425  |  0                        0                       -0.060576148825500558762
!<  0.642615758240322548157  |  0                        0                       -0.144942775902865915672
!<  0.357384241759677451842  |  0                        0                       -0.076881336420335693858
!<  0.117472338035267653574  |  0                        0                       -0.052032968680060307651
!<  0.833333333333333333333  |  0.330885963040722183948  0.489662957309450192844 -1.378964865748435675821
!<  0.309036761204472681298  | -0.072954784731363262918  0                       -0.851236239662007619739
!<  0.539357840802981787532  |  0                       -0.777333643644968233538  0
!<  0.1                      |  0                        0                        0
!<  1                        |  0                        0.259111214548322744512 -0.358278966717952089048
!< ------------------------------------------------------------------------------------------------------
!<                           |  0                        0.05                     0
!< ...continued...
!<  0                        |
!<  0.1                      |
!<  0.539357840802981787532  |
!<  0.809036761204472681298  |
!<  0.309036761204472681298  |
!<  0.981074190219795268254  |
!<  0.833333333333333333333  |
!<  0.354017365856802376329  | -0.059578021181736100156
!<  0.882527661964732346425  |  0.321763705601778390100  0.510485725608063031577
!<  0.642615758240322548157  | -0.333269719096256706589  0.499269229556880061353  0.509504608929686104236
!<  0.357384241759677451842  |  0.239527360324390649107  0.397774662368094639047  0.010755895687360745555
!<  0.117472338035267653574  | -0.057695414616854888173  0.194781915712104164976  0.145384923188325069727
!<  0.833333333333333333333  | -0.861164195027635666673  5.784288136375372200229  3.288077619851035668904
!<  0.309036761204472681298  |  0.398320112318533301719  3.639372631810356060294  1.548228770398303223653
!<  0.539357840802981787532  | -0.091089566215517606959  0                        0
!<  0.1                      |  0                        0                        0
!<  1                        | -1.045948959408833060950  0.930327845415626983292  1.779509594317081024461
!< ------------------------------------------------------------------------------------------------------
!<                           |  0.04                     0                        0.189237478148923490158
!< ...continued...
!<  0                        |
!<  0.1                      |
!<  0.539357840802981787532  |
!<  0.809036761204472681298  |
!<  0.309036761204472681298  |
!<  0.981074190219795268254  |
!<  0.833333333333333333333  |
!<  0.354017365856802376329  |
!<  0.882527661964732346425  |
!<  0.642615758240322548157  |
!<  0.357384241759677451842  | -0.327769124164018874147
!<  0.117472338035267653574  | -0.078294271035167077755 -0.114503299361098912184
!<  0.833333333333333333333  | -2.386339050931363840134 -3.254793424836439186545 -2.163435416864229823539
!<  0.309036761204472681298  | -2.122217147040537160260 -1.583503985453261727133 -1.715616082859362649220
!<  0.539357840802981787532  |  0                        0                        0
!<  0.1                      |  0                        0                        0
!<  1                        |  0.1                     -0.282547569539044081612 -0.159327350119972549169
!< ------------------------------------------------------------------------------------------------------
!<                           |  0.277429188517743176508  0.277429188517743176508  0.189237478148923490158
!< ...continued...
!<  0                        |
!<  0.1                      |
!<  0.539357840802981787532  |
!<  0.809036761204472681298  |
!<  0.309036761204472681298  |
!<  0.981074190219795268254  |
!<  0.833333333333333333333  |
!<  0.354017365856802376329  |
!<  0.882527661964732346425  |
!<  0.642615758240322548157  |
!<  0.357384241759677451842  |
!<  0.117472338035267653574  |
!<  0.833333333333333333333  |
!<  0.309036761204472681298  | -0.024403640575012745213
!<  0.539357840802981787532  |  0.091089566215517606959  0.777333643644968233538
!<  0.1                      |  0                        0                        0.157178665799771163367
!<  1                        | -0.145515894647001510860 -0.259111214548322744512 -0.342758159847189839942 -0.675
!< -------------------------------------------------------------------------------------------------------------
!<                           | -0.04                    -0.05                    -0.033333333333333333333 -0.025
!< ...continued...
!<  0                        |
!<  0.1                      |
!<  0.539357840802981787532  |
!<  0.809036761204472681298  |
!<  0.309036761204472681298  |
!<  0.981074190219795268254  |
!<  0.833333333333333333333  |
!<  0.354017365856802376329  |
!<  0.882527661964732346425  |
!<  0.642615758240322548157  |
!<  0.357384241759677451842  |
!<  0.117472338035267653574  |
!<  0.833333333333333333333  |
!<  0.309036761204472681298  |
!<  0.539357840802981787532  |
!<  0.1                      |
!<  1                        |
!< ---------------------------------------------------
!<                           | 0.033333333333333333333
!<```
!<
!<#### Bibliography
!<
!< [1] *A family of embedded Runge-Kutta formulae*, Dormand, J. R., Prince, P. J. (1980), Journal of Computational and
!< Applied Mathematics 6 (1): 19--26, doi:10.1016/0771-050X(80)90013-3.
!<
!< [2] *A New Embedded Pair of Runge-Kutta Formulas of orders 5 and 6*, M. Calvo, J.I. Montijano, L. Randez, Computers & Mathematics
!< with Applications, Volume 20, Issue 1, 1990, Pages 15--24, ISSN 0898-1221, http://dx.doi.org/10.1016/0898-1221(90)90064-Q.
!<
!< [3] *A variable order Runge-Kutta method for initial value problems with rapidly varying right-hand sides*, J. R. Cash,
!< A. H. Karp, ACM Transactions on Mathematical Software, vol. 16,  pp. 201--222, 1990, doi:10.1145/79505.79507.
!<
!< [4] *A tenth-order Runge-Kutta method with error estimate*, Feagin, T., Proceedings of the IAENG Conf. on Scientific
!< Computing. 2007.
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
  !<
  !< ### List of errors status
  !<+ error=0 => no error;
  !<+ error=1 => bad (unsupported) number of required time steps;
  real(R_P)              :: tolerance=0._R_P !< Tolerance on the local truncation error.
  real(R_P)              :: pp1_inv=0._R_P   !< 1/(p+1) where p is the accuracy order of the lower accurate scheme of the pair.
  integer(I_P)           :: stages=0         !< Number of stages.
  real(R_P), allocatable :: alph(:,:)        !< \(\alpha\) Butcher's coefficients.
  real(R_P), allocatable :: beta(:,:)        !< \(\beta\) Butcher's coefficients.
  real(R_P), allocatable :: gamm(:)          !< \(\gamma\) Butcher's coefficients.
  integer(I_P)           :: error=0          !< Error status flag: trap occurrences of errors.
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
  case(2)
    ! HERK(2,2)
    self%pp1_inv = 1._R_P/(2._R_P + 1._R_P)

    self%beta(1, 1) =  0.5_R_P ; self%beta(1, 2) =  5._R_P
    self%beta(2, 1) =  1._R_P  ; self%beta(2, 2) =  0._R_P

    self%alph(2, 1) = 1._R_P

    self%gamm(2) = 1._R_P
  case(6)
    ! CKRK(6,5)
    self%pp1_inv = 1._R_P/(5._R_P + 1._R_P)

    self%beta(1, 1) = 37._R_P/378._R_P   ; self%beta(1, 2) = 2825._R_P/27648._R_P
    self%beta(2, 1) = 0._R_P             ; self%beta(2, 2) = 0._R_P
    self%beta(3, 1) = 250._R_P/621._R_P  ; self%beta(3, 2) = 18575._R_P/48384._R_P
    self%beta(4, 1) = 125._R_P/594._R_P  ; self%beta(4, 2) = 13525._R_P/55296._R_P
    self%beta(5, 1) = 0._R_P             ; self%beta(5, 2) = 277._R_P/14336._R_P
    self%beta(6, 1) = 512._R_P/1771._R_P ; self%beta(6, 2) = 1._R_P/4._R_P

    self%alph(2, 1)=1._R_P/5._R_P
    self%alph(3, 1)=3._R_P/40._R_P      ;self%alph(3, 2)= 9._R_P/40._R_P
    self%alph(4, 1)=3._R_P/10._R_P      ;self%alph(4, 2)= -9._R_P/10._R_P ;self%alph(4, 3)=6._R_P/5._R_P
    self%alph(5, 1)=-11._R_P/54._R_P    ;self%alph(5, 2)= 5._R_P/2._R_P   ;self%alph(5, 3)=-70._R_P/27._R_P
    self%alph(6, 1)=1631._R_P/55296._R_P;self%alph(6, 2)=175._R_P/512._R_P;self%alph(6, 3)=575._R_P/13824._R_P

    self%alph(5, 4)=35._R_P/27._R_P
    self%alph(6, 4)=44275._R_P/110592._R_P;self%alph(6, 5)=253._R_P/4096._R_P

    self%gamm(2) = 1._R_P/5._R_P
    self%gamm(3) = 3._R_P/10._R_P
    self%gamm(4) = 3._R_P/5._R_P
    self%gamm(5) = 1._R_P
    self%gamm(6) = 7._R_P/8._R_P
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
  case(9)
    ! CMRK(9,6)
    self%pp1_inv = 1._R_P/(6._R_P + 1._R_P)

    self%beta(1, 1) = 17572349._R_P/289262523._R_P  ; self%beta(1, 2) = 15231665._R_P/510830334._R_P
    self%beta(2, 1) = 0._R_P                        ; self%beta(2, 2) = 0._R_P
    self%beta(3, 1) = 57513011._R_P/201864250._R_P  ; self%beta(3, 2) = 59452991._R_P/116050448._R_P
    self%beta(4, 1) = 15587306._R_P/354501571._R_P  ; self%beta(4, 2) = -28398517._R_P/122437738._R_P
    self%beta(5, 1) = 71783021._R_P/234982865._R_P  ; self%beta(5, 2) = 56673824._R_P/137010559._R_P
    self%beta(6, 1) = 29672000._R_P/180480167._R_P  ; self%beta(6, 2) = 68003849._R_P/426673583._R_P
    self%beta(7, 1) = 65567621._R_P/127060952._R_P  ; self%beta(7, 2) = 7097631._R_P/37564021._R_P
    self%beta(8, 1) = -79074570._R_P/210557597._R_P ; self%beta(8, 2) = -71226429._R_P/583093742._R_P
    self%beta(9, 1) = 0._R_P                        ; self%beta(9, 2) = 1._R_P/20._R_P

    self%alph(2, 1)=2._R_P/15._R_P
    self%alph(3, 1)=1._R_P/20._R_P               ; self%alph(3, 2)=3._R_P/20._R_P
    self%alph(4, 1)=3._R_P/40._R_P               ; self%alph(4, 2)=0._R_P
    self%alph(5, 1)=86727015._R_P/196851553._R_P ; self%alph(5, 2)=-60129073._R_P/52624712._R_P
    self%alph(6, 1)=-86860849._R_P/45628967._R_P ; self%alph(6, 2)=111022885._R_P/25716487._R_P
    self%alph(7, 1)=77759591._R_P/16096467._R_P  ; self%alph(7, 2)=-49252809._R_P/6452555._R_P
    self%alph(8, 1)=237564263._R_P/39280295._R_P ; self%alph(8, 2)=-100523239._R_P/10677940._R_P
    self%alph(9, 1)=17572349._R_P/289262523._R_P ; self%alph(9, 2)=0._R_P

    self%alph(4, 3)=9._R_P/40._R_P
    self%alph(5, 3)=957436434._R_P/1378352377._R_P ; self%alph(5, 4)=83886832._R_P/147842441._R_P
    self%alph(6, 3)=108046682._R_P/101167669._R_P  ; self%alph(6, 4)=-141756746._R_P/36005461._R_P
    self%alph(7, 3)=-381680111._R_P/51572984._R_P  ; self%alph(7, 4)=879269579._R_P/66788831._R_P
    self%alph(8, 3)=-265574846._R_P/27330247._R_P  ; self%alph(8, 4)=317978411._R_P/18988713._R_P
    self%alph(9, 3)=57513011._R_P/201864250._R_P   ; self%alph(9, 4)=15587306._R_P/354501571._R_P

    self%alph(6, 5)=73139862._R_P/60170633._R_P
    self%alph(7, 5)=-90453121._R_P/33722162._R_P  ; self%alph(7, 6)=111179552._R_P/157155827._R_P
    self%alph(8, 5)=-124494385._R_P/35453627._R_P ; self%alph(8, 6)=86822444._R_P/100138635._R_P
    self%alph(9, 5)=71783021._R_P/234982865._R_P  ; self%alph(9, 6)=29672000._R_P/180480167._R_P

    self%alph(8, 7)=-12873523._R_P/724232625._R_P
    self%alph(9, 7)=65567621._R_P/127060952._R_P ;self%alph(9, 8)=-79074570._R_P/210557597._R_P

    self%gamm(2) = 2._R_P/15._R_P
    self%gamm(3) = 1._R_P/5._R_P
    self%gamm(4) = 3._R_P/10._R_P
    self%gamm(5) = 14._R_P/25._R_P
    self%gamm(6) = 19._R_P/25._R_P
    self%gamm(7) = 35226607._R_P/35688279._R_P
    self%gamm(8) = 1._R_P
    self%gamm(9) = 1._R_P
  case(17)
    ! FRK(17,10)
    self%pp1_inv = 1._R_P/(10._R_P + 1._R_P)

    self%beta(1 , 1) = 0.033333333333333333333_R_P; self%beta(1 , 2) = 0.033333333333333333333_R_P
    self%beta(2 , 1) = 0.025_R_P                  ; self%beta(2 , 2) = 1._R_P/36._R_P
    self%beta(3 , 1) = 0.033333333333333333333_R_P; self%beta(3 , 2) = 0.033333333333333333333_R_P
    self%beta(4 , 1) = 0._R_P                     ; self%beta(4 , 2) = 0._R_P
    self%beta(5 , 1) = 0.05_R_P                   ; self%beta(5 , 2) = 0.05_R_P
    self%beta(6 , 1) = 0._R_P                     ; self%beta(6 , 2) = 0._R_P
    self%beta(7 , 1) = 0.04_R_P                   ; self%beta(7 , 2) = 0.04_R_P
    self%beta(8 , 1) = 0._R_P                     ; self%beta(8 , 2) = 0._R_P
    self%beta(9 , 1) = 0.189237478148923490158_R_P; self%beta(9 , 2) = 0.189237478148923490158_R_P
    self%beta(10, 1) = 0.277429188517743176508_R_P; self%beta(10, 2) = 0.277429188517743176508_R_P
    self%beta(11, 1) = 0.277429188517743176508_R_P; self%beta(11, 2) = 0.277429188517743176508_R_P
    self%beta(12, 1) = 0.189237478148923490158_R_P; self%beta(12, 2) = 0.189237478148923490158_R_P
    self%beta(13, 1) =-0.04_R_P                   ; self%beta(13, 2) =-0.04_R_P
    self%beta(14, 1) =-0.05_R_P                   ; self%beta(14, 2) =-0.05_R_P
    self%beta(15, 1) =-0.033333333333333333333_R_P; self%beta(15, 2) =-0.033333333333333333333_R_P
    self%beta(16, 1) =-0.025_R_P                  ; self%beta(16, 2) =-1._R_P/36._R_P
    self%beta(17, 1) = 0.033333333333333333333_R_P; self%beta(17, 2) = 0.033333333333333333333_R_P

    self%alph(2 , 1)= 0.1_R_P
    self%alph(3 , 1)=-0.915176561375291440520_R_P; self%alph(3 , 2)=1.454534402178273228052_R_P
    self%alph(4 , 1)= 0.202259190301118170324_R_P; self%alph(4 , 2)=0._R_P
    self%alph(5 , 1)= 0.184024714708643575149_R_P; self%alph(5 , 2)=0._R_P
    self%alph(6 , 1)= 0.087900734020668133731_R_P; self%alph(6 , 2)=0._R_P
    self%alph(7 , 1)= 0.085970050490246030218_R_P; self%alph(7 , 2)=0._R_P
    self%alph(8 , 1)= 0.120930449125333720660_R_P; self%alph(8 , 2)=0._R_P
    self%alph(9 , 1)= 0.110854379580391483508_R_P; self%alph(9 , 2)=0._R_P
    self%alph(10, 1)= 0.112054414752879004829_R_P; self%alph(10, 2)=0._R_P
    self%alph(11, 1)= 0.113976783964185986138_R_P; self%alph(11, 2)=0._R_P
    self%alph(12, 1)= 0.079831452828019604635_R_P; self%alph(12, 2)=0._R_P
    self%alph(13, 1)= 0.985115610164857280120_R_P; self%alph(13, 2)=0._R_P
    self%alph(14, 1)= 0.895080295771632891049_R_P; self%alph(14, 2)=0._R_P
    self%alph(15, 1)=-0.915176561375291440520_R_P; self%alph(15, 2)=1.454534402178273228052_R_P
    self%alph(16, 1)= 0.1_R_P                    ; self%alph(16, 2)=0._R_P
    self%alph(17, 1)= 0.181781300700095283888_R_P; self%alph(17, 2)=0.675_R_P

    self%alph(4 , 3)= 0.606777570903354510974_R_P
    self%alph(5 , 3)= 0.197966831227192369068_R_P; self%alph(5 , 4)=-0.072954784731363262918_R_P
    self%alph(6 , 3)= 0._R_P                     ; self%alph(6 , 4)= 0.410459702520260645318_R_P
    self%alph(7 , 3)= 0._R_P                     ; self%alph(7 , 4)= 0.330885963040722183948_R_P
    self%alph(8 , 3)= 0._R_P                     ; self%alph(8 , 4)= 0._R_P
    self%alph(9 , 3)= 0._R_P                     ; self%alph(9 , 4)= 0._R_P
    self%alph(10, 3)= 0._R_P                     ; self%alph(10, 4)= 0._R_P
    self%alph(11, 3)= 0._R_P                     ; self%alph(11, 4)= 0._R_P
    self%alph(12, 3)= 0._R_P                     ; self%alph(12, 4)= 0._R_P
    self%alph(13, 3)= 0._R_P                     ; self%alph(13, 4)= 0.330885963040722183948_R_P
    self%alph(14, 3)= 0.197966831227192369068_R_P; self%alph(14, 4)=-0.072954784731363262918_R_P
    self%alph(15, 3)= 0._R_P                     ; self%alph(15, 4)= 0._R_P
    self%alph(16, 3)=-0.157178665799771163367_R_P; self%alph(16, 4)= 0._R_P
    self%alph(17, 3)= 0.342758159847189839942_R_P; self%alph(17, 4)= 0._R_P

    self%alph(6 , 5)= 0.482713753678866489204_R_P
    self%alph(7 , 5)= 0.489662957309450192844_R_P; self%alph(7 , 6)=-0.073185637507085073678_R_P
    self%alph(8 , 5)= 0.260124675758295622809_R_P; self%alph(8 , 6)= 0.032540262154909133015_R_P
    self%alph(9 , 5)= 0._R_P                     ; self%alph(9 , 6)=-0.060576148825500558762_R_P
    self%alph(10, 5)= 0._R_P                     ; self%alph(10, 6)=-0.144942775902865915672_R_P
    self%alph(11, 5)= 0._R_P                     ; self%alph(11, 6)=-0.076881336420335693858_R_P
    self%alph(12, 5)= 0._R_P                     ; self%alph(12, 6)=-0.052032968680060307651_R_P
    self%alph(13, 5)= 0.489662957309450192844_R_P; self%alph(13, 6)=-1.378964865748435675821_R_P
    self%alph(14, 5)= 0._R_P                     ; self%alph(14, 6)=-0.851236239662007619739_R_P
    self%alph(15, 5)=-0.777333643644968233538_R_P; self%alph(15, 6)= 0._R_P
    self%alph(16, 5)= 0._R_P                     ; self%alph(16, 6)= 0._R_P
    self%alph(17, 5)= 0.259111214548322744512_R_P; self%alph(17, 6)=-0.358278966717952089048_R_P

    self%alph(8 , 7)=-0.059578021181736100156_R_P
    self%alph(9 , 7)= 0.321763705601778390100_R_P; self%alph(9 , 8)=0.510485725608063031577_R_P
    self%alph(10, 7)=-0.333269719096256706589_R_P; self%alph(10, 8)=0.499269229556880061353_R_P
    self%alph(11, 7)= 0.239527360324390649107_R_P; self%alph(11, 8)=0.397774662368094639047_R_P
    self%alph(12, 7)=-0.057695414616854888173_R_P; self%alph(12, 8)=0.194781915712104164976_R_P
    self%alph(13, 7)=-0.861164195027635666673_R_P; self%alph(13, 8)=5.784288136375372200229_R_P
    self%alph(14, 7)= 0.398320112318533301719_R_P; self%alph(14, 8)=3.639372631810356060294_R_P
    self%alph(15, 7)=-0.091089566215517606959_R_P; self%alph(15, 8)=0._R_P
    self%alph(16, 7)= 0._R_P                     ; self%alph(16, 8)=0._R_P
    self%alph(17, 7)=-1.045948959408833060950_R_P; self%alph(17, 8)=0.930327845415626983292_R_P

    self%alph(10, 9)=0.509504608929686104236_R_P
    self%alph(11, 9)=0.010755895687360745555_R_P; self%alph(11, 10)=-0.327769124164018874147_R_P
    self%alph(12, 9)=0.145384923188325069727_R_P; self%alph(12, 10)=-0.078294271035167077755_R_P
    self%alph(13, 9)=3.288077619851035668904_R_P; self%alph(13, 10)=-2.386339050931363840134_R_P
    self%alph(14, 9)=1.548228770398303223653_R_P; self%alph(14, 10)=-2.122217147040537160260_R_P
    self%alph(15, 9)=0._R_P                     ; self%alph(15, 10)= 0._R_P
    self%alph(16, 9)=0._R_P                     ; self%alph(16, 10)= 0._R_P
    self%alph(17, 9)=1.779509594317081024461_R_P; self%alph(17, 10)= 0.1_R_P

    self%alph(12, 11)=-0.114503299361098912184_R_P
    self%alph(13, 11)=-3.254793424836439186545_R_P; self%alph(13, 12)=-2.163435416864229823539_R_P
    self%alph(14, 11)=-1.583503985453261727133_R_P; self%alph(14, 12)=-1.715616082859362649220_R_P
    self%alph(15, 11)= 0._R_P                     ; self%alph(15, 12)= 0._R_P
    self%alph(16, 11)= 0._R_P                     ; self%alph(16, 12)= 0._R_P
    self%alph(17, 11)=-0.282547569539044081612_R_P; self%alph(17, 12)=-0.159327350119972549169_R_P

    self%alph(14, 13)=-0.024403640575012745213_R_P
    self%alph(15, 13)= 0.091089566215517606959_R_P; self%alph(15, 14)= 0.777333643644968233538_R_P
    self%alph(16, 13)= 0._R_P                     ; self%alph(16, 14)= 0._R_P
    self%alph(17, 13)=-0.145515894647001510860_R_P; self%alph(17, 14)=-0.259111214548322744512_R_P

    self%alph(16, 15)= 0.157178665799771163367_R_P
    self%alph(17, 15)=-0.342758159847189839942_R_P; self%alph(17, 16)=-0.675_R_P

    self%gamm(2 ) = 0.1_R_P
    self%gamm(3 ) = 0.539357840802981787532_R_P
    self%gamm(4 ) = 0.809036761204472681298_R_P
    self%gamm(5 ) = 0.309036761204472681298_R_P
    self%gamm(6 ) = 0.981074190219795268254_R_P
    self%gamm(7 ) = 0.833333333333333333333_R_P
    self%gamm(8 ) = 0.354017365856802376329_R_P
    self%gamm(9 ) = 0.882527661964732346425_R_P
    self%gamm(10) = 0.642615758240322548157_R_P
    self%gamm(11) = 0.357384241759677451842_R_P
    self%gamm(12) = 0.117472338035267653574_R_P
    self%gamm(13) = 0.833333333333333333333_R_P
    self%gamm(14) = 0.309036761204472681298_R_P
    self%gamm(15) = 0.539357840802981787532_R_P
    self%gamm(16) = 0.1_R_P
    self%gamm(17) = 1._R_P
  case default
    ! bad (unsupported) number of required time steps
    self%error = 1
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
  self%error = 0
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
