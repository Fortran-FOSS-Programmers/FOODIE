!< FOODIE integrator: provide an explicit class of SSP Runge-Kutta schemes, from 1st to 4th order accurate.

module foodie_integrator_runge_kutta_ssp
!< FOODIE integrator: provide an explicit class of SSP Runge-Kutta schemes, from 1st to 4th order accurate.
!<
!< The integrators provided have the Total Variation Diminishing (TVD) property or the Strong Stability Preserving (SSP)
!< one. The schemes are explicit and defined through the Butcher's table syntax, see[1] .
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
!<```
!<  gamma^1    | alpha^{1,1}       alpha^{1,2}       ...        alpha^{1,Ns}
!<  gamma^2    | alpha^{2,1}       alpha^{2,2}       ...        alpha^{2,Ns}
!<  .          | .                 .                 .          .
!<  .          | .                 .                  .         .
!<  .          | .                 .                   .        .
!<  gamma^{Ns} | alpha^{Ns,1}      alpha^{Ns,2}      ...        alpha^{Ns,Ns}
!< ------------|-------------------------------------------------------------
!<             | beta^1            beta^2            ...        beta^{Ns}
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
!<             | beta^1            beta^2            ...        beta^{Ns}
!<```
!<
!< Moreover the following relation always holds:
!< \( \gamma^s = \sum_{i=1}^{Ns}\alpha^{s,i} \)
!<
!< The different schemes are selected accordingly to the number of stages used. Currently the following schemes are available:
!<
!<##### 1 stage, Explicit Forward Euler, 1st order
!< This scheme is TVD and reverts to Explicit Forward Euler, it being 1st order.
!<```
!<  0 | 0
!< ---|---
!<    | 1
!<```
!<
!<##### 2 stages, SSP, 2nd order
!< This scheme is an optmial SSP(2, 2) without low-storage algorithm, see [2].
!<```
!<  0 | 0     0
!<  1 | 1     0
!< ---|-----------
!<    | 1/2   1/2
!<```
!<
!<##### 3 stages, SSP, 3rd order
!< This scheme is an optmial SSP(3, 3) without low-storage algorithm, see [2].
!<```
!<  0   | 0     0     0
!<  1   | 1     0     0
!<  1/2 | 1/4   1/4   0
!< -----|-----------------
!<      | 1/6   1/6   1/3
!<```
!<
!<##### 5 stages, SSP, 4th order
!< This scheme is an optmial SSP(5, 4) without low-storage algorithm, see [2].
!<```
!<  0                | 0                  0                  0                  0                  0
!<  0.39175222700392 | 0.39175222700392   0                  0                  0                  0
!<  0.58607968896780 | 0.21766909633821   0.36841059262959   0                  0                  0
!<  0.47454236302687 | 0.08269208670950   0.13995850206999   0.25189177424738   0                  0
!<  0.93501063100924 | 0.06796628370320   0.11503469844438   0.20703489864929   0.54497475021237   0
!< ------------------|---------------------------------------------------------------------------------------------
!<                   | 0.14681187618661   0.24848290924556   0.10425883036650   0.27443890091960   0.22600748319395
!<```
!<
!<#### Bibliography
!<
!< [1] *Coefficients for the study of Runge-Kutta integration processes*, Butcher, J.C., J. Austral. Math. Soc., Vol. 3,
!< pages: 185--201, 1963.
!<
!< [2] *High Order Strong Stability Preserving Time Discretizations*, Gottlieb, S., Ketcheson, D. I., Shu, C.W., Journal of
!< Scientific Computing, vol. 38, N. 3, 2009, pp. 251-289.

use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrand_object, only : integrand_object
use foodie_integrator_multistage_object, only : integrator_multistage_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_runge_kutta_ssp

character(len=99), parameter :: class_name_='runge_kutta_ssp'                                    !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:4)=[trim(class_name_)//'_stages_1_order_1', &
                                                         trim(class_name_)//'_stages_2_order_2', &
                                                         trim(class_name_)//'_stages_3_order_3', &
                                                         trim(class_name_)//'_stages_5_order_4'] !< List of supported schemes.

logical, parameter :: has_fast_mode_=.true. !< Flag to check if integrator provides *fast mode* integrate.

type, extends(integrator_multistage_object) :: integrator_runge_kutta_ssp
  !< FOODIE integrator: provide an explicit class of SSP Runge-Kutta schemes, from 1st to 4th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the RK coefficients) before used.
  private
  real(R_P), allocatable :: alph(:,:) !< \(\alpha\) Butcher's coefficients.
  real(R_P), allocatable :: beta(:)   !< \(\beta\) Butcher's coefficients.
  real(R_P), allocatable :: gamm(:)   !< \(\gamma\) Butcher's coefficients.
  contains
    ! deferred methods
    procedure, pass(self) :: class_name           !< Return the class name of schemes.
    procedure, pass(self) :: has_fast_mode        !< Return .true. if the integrator class has *fast mode* integrate.
    procedure, pass(lhs)  :: integr_assign_integr !< Operator `=`.
    procedure, pass(self) :: integrate            !< Integrate integrand field.
    procedure, pass(self) :: integrate_fast       !< Integrate integrand field, fast mode.
    procedure, pass(self) :: is_supported         !< Return .true. if the integrator class support the given scheme.
    procedure, pass(self) :: supported_schemes    !< Return the list of supported schemes.
    ! public methods
    procedure, pass(self) :: destroy    !< Destroy the integrator.
    procedure, pass(self) :: initialize !< Initialize (create) the integrator.
endtype integrator_runge_kutta_ssp

contains
   ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_runge_kutta_ssp), intent(in) :: self       !< Integrator.
  character(len=99)                             :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  elemental function has_fast_mode(self)
  !< Return .true. if the integrator class has *fast mode* integrate.
  class(integrator_runge_kutta_ssp), intent(in) :: self          !< Integrator.
  logical                                       :: has_fast_mode !< Inquire result.

  has_fast_mode = has_fast_mode_
  endfunction has_fast_mode

  pure subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_runge_kutta_ssp), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),          intent(in)    :: rhs !< Right hand side.

  call lhs%assign_abstract(rhs=rhs)
  select type(rhs)
  class is (integrator_runge_kutta_ssp)
                             lhs%stages = rhs%stages
    if (allocated(rhs%alph)) lhs%alph   = rhs%alph
    if (allocated(rhs%beta)) lhs%beta   = rhs%beta
    if (allocated(rhs%gamm)) lhs%gamm   = rhs%gamm
  endselect
  endsubroutine integr_assign_integr

  subroutine integrate(self, U, Dt, t, new_Dt)
  !< Integrate field with explicit SSP Runge-Kutta scheme.
  class(integrator_runge_kutta_ssp), intent(inout) :: self      !< Integrator.
  class(integrand_object),           intent(inout) :: U         !< Field to be integrated.
  real(R_P),                         intent(in)    :: Dt        !< Time step.
  real(R_P),                         intent(in)    :: t         !< Time.
  real(R_P), optional,               intent(out)   :: new_Dt    !< New adapted time step.
  integer(I_P)                                     :: s         !< First stages counter.
  integer(I_P)                                     :: ss        !< Second stages counter.

  ! computing stages
  do s=1, self%stages
     self%stage(s) = U
     do ss=1, s - 1
        self%stage(s) = self%stage(s) + (self%stage(ss) * (Dt * self%alph(s, ss)))
     enddo
     self%stage(s) = self%stage(s)%t(t=t + self%gamm(s) * Dt)
  enddo
  ! computing new time step
  do s=1, self%stages
     U = U + (self%stage(s) * (Dt * self%beta(s)))
  enddo
  if (present(new_Dt)) new_Dt = Dt
  endsubroutine integrate

  subroutine integrate_fast(self, U, Dt, t, new_Dt)
  !< Integrate field with explicit SSP Runge-Kutta scheme.
  class(integrator_runge_kutta_ssp), intent(inout) :: self   !< Integrator.
  class(integrand_object),           intent(inout) :: U      !< Field to be integrated.
  real(R_P),                         intent(in)    :: Dt     !< Time step.
  real(R_P),                         intent(in)    :: t      !< Time.
  real(R_P), optional,               intent(out)   :: new_Dt !< New adapted time step.
  integer(I_P)                                     :: s      !< First stages counter.
  integer(I_P)                                     :: ss     !< Second stages counter.

  ! computing stages
  self%buffer = U
  do s=1, self%stages
     self%stage(s) = U
     do ss=1, s - 1
        call self%buffer%multiply_fast(lhs=self%stage(ss), rhs=Dt * self%alph(s, ss))
        call self%stage(s)%add_fast(lhs=self%stage(s), rhs=self%buffer)
     enddo
     call self%stage(s)%t_fast(t=t + self%gamm(s) * Dt)
  enddo
  ! computing new time step
  do s=1, self%stages
     call self%buffer%multiply_fast(lhs=self%stage(s), rhs=Dt * self%beta(s))
     call U%add_fast(lhs=U, rhs=self%buffer)
  enddo
  if (present(new_Dt)) new_Dt = Dt
  endsubroutine integrate_fast

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_runge_kutta_ssp), intent(in) :: self         !< Integrator.
  character(*),                      intent(in) :: scheme       !< Selected scheme.
  logical                                       :: is_supported !< Inquire result.
  integer(I_P)                                  :: s            !< Counter.

  is_supported = .false.
  do s=lbound(supported_schemes_, dim=1), ubound(supported_schemes_, dim=1)
    if (trim(adjustl(scheme)) == trim(adjustl(supported_schemes_(s)))) then
      is_supported = .true.
      return
    endif
  enddo
  endfunction is_supported

  pure function supported_schemes(self) result(schemes)
  !< Return the list of supported schemes.
  class(integrator_runge_kutta_ssp), intent(in) :: self       !< Integrator.
  character(len=99), allocatable                :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_runge_kutta_ssp), intent(inout) :: self !< Integrator.

  call self%destroy_multistage
  if (allocated(self%alph)) deallocate(self%alph)
  if (allocated(self%beta)) deallocate(self%beta)
  if (allocated(self%gamm)) deallocate(self%gamm)
  endsubroutine destroy

  subroutine initialize(self, scheme, U, stop_on_fail)
  !< Create the actual RK integrator: initialize the Butcher' table coefficients.
  class(integrator_runge_kutta_ssp), intent(inout)        :: self         !< Integrator.
  character(*),                      intent(in)           :: scheme       !< Selected scheme.
  class(integrand_object),           intent(in), optional :: U            !< Integrand molding prototype.
  logical,                           intent(in), optional :: stop_on_fail !< Stop execution if initialization fail.

  if (self%is_supported(scheme=scheme)) then
    call self%destroy
    self%description_ = trim(adjustl(scheme))
    select case(trim(adjustl(scheme)))
    case('runge_kutta_ssp_stages_1_order_1')
      self%stages = 1
      allocate(self%beta(1:self%stages               )) ; self%beta = 0._R_P
      allocate(self%alph(1:self%stages, 1:self%stages)) ; self%alph = 0._R_P
      allocate(self%gamm(               1:self%stages)) ; self%gamm = 0._R_P

      self%beta(1) = 1._R_P
    case('runge_kutta_ssp_stages_2_order_2')
      self%stages = 2
      allocate(self%beta(1:self%stages               )) ; self%beta = 0._R_P
      allocate(self%alph(1:self%stages, 1:self%stages)) ; self%alph = 0._R_P
      allocate(self%gamm(               1:self%stages)) ; self%gamm = 0._R_P

      self%beta(1) = 0.5_R_P
      self%beta(2) = 0.5_R_P

      self%alph(2, 1) = 1._R_P

      self%gamm(2) = 1._R_P
    case('runge_kutta_ssp_stages_3_order_3')
      self%stages = 3
      allocate(self%beta(1:self%stages               )) ; self%beta = 0._R_P
      allocate(self%alph(1:self%stages, 1:self%stages)) ; self%alph = 0._R_P
      allocate(self%gamm(               1:self%stages)) ; self%gamm = 0._R_P

      self%beta(1) = 1._R_P/6._R_P
      self%beta(2) = 1._R_P/6._R_P
      self%beta(3) = 2._R_P/3._R_P

      self%alph(2, 1) = 1._R_P
      self%alph(3, 1) = 0.25_R_P ; self%alph(3, 2) = 0.25_R_P

      self%gamm(2) = 1._R_P
      self%gamm(3) = 0.5_R_P
    case('runge_kutta_ssp_stages_5_order_4')
      self%stages = 5
      allocate(self%beta(1:self%stages               )) ; self%beta = 0._R_P
      allocate(self%alph(1:self%stages, 1:self%stages)) ; self%alph = 0._R_P
      allocate(self%gamm(               1:self%stages)) ; self%gamm = 0._R_P

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
    self%registers = self%stages
    if (present(U)) call self%allocate_integrand_members(U=U)
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=stop_on_fail)
  endif
  endsubroutine initialize
endmodule foodie_integrator_runge_kutta_ssp
