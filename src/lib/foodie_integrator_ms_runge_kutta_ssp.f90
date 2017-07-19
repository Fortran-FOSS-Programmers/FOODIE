!< FOODIE integrator: provide an explicit class of Multi-step Runge-Kutta Methods with Strong Stability Preserving property, from
!< 2nd to 3rd order accurate.

module foodie_integrator_ms_runge_kutta_ssp
!< FOODIE integrator: provide an explicit class of Multi-step Runge-Kutta Methods with Strong Stability Preserving property, from
!< 2nd to 3rd order accurate.
!<
!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the LMM-SSP class scheme implemented is:
!<
!< $$
!< \begin{align}
!< y_1^n & = u^n \\
!< y_i^n & = \sum_{l=1}^{k} d_{il} u^{n-k+l} + \Delta{t}\sum_{l=1}^{k-1} \hat{a}_{il} F(u^{n-k+l}) +
!<                                             \Delta{t}\sum_{j=1}^{i-1} a_{ij} F(y_j^n) \; \; \; \;  2 \leq i \leq s \\
!< u^{n+1} & = \sum_{l=1}^{k} \theta_l u^{n-k+l} + \Delta{t}\sum_{l=1}^{k-1} \hat{b}_{l} F(u^{n-k+l}) +
!<                                                  \Delta{t}\sum_{j=1}^s b_j F(y_j^n).
!< \end{align}
!< $$
!<
!<where \(N_s\) is the number of previous steps considered.
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit. The coefficients *a,b* define the actual scheme, that is selected accordingly to the number of
!< **steps** used. Currently, the schemes provided have steps number in *[3, 5]*. The formal order of accuracy varies
!< consistently in *[2nd, 3rd]* order.
!<
!<#### Bibliography
!< [1] *Explicit Strong Stability Preserving Multistep Runge-Kutta Methods*, C. Bresten, S. Gottlieb, Z. Grant, D. Higgs,
!< D. Ketcheson, A. Németh, 2016, Mathematics of Computations,

use foodie_error_codes, only : ERROR_UNSUPPORTED_SCHEME
use foodie_integrand_object, only : integrand_object
use foodie_integrator_multistage_multistep_object, only : integrator_multistage_multistep_object
use foodie_integrator_object, only : integrator_object
use penf, only : I_P, R_P

implicit none
private
public :: integrator_ms_runge_kutta_ssp

character(len=99), parameter :: class_name_='ms_runge_kutta_ssp' !< Name of the class of schemes.
character(len=99), parameter :: supported_schemes_(1:3)=[trim(class_name_)//'_steps_2_stages_2_order_3', &
                                                         trim(class_name_)//'_steps_3_stages_2_order_3', &
                                                         trim(class_name_)//'_steps_4_stages_5_order_8'] !< List of supported
                                                                                                         !< schemes.

logical, parameter :: has_fast_mode_=.true. !< Flag to check if integrator provides *fast mode* integrate.

type, extends(integrator_multistage_multistep_object) :: integrator_ms_runge_kutta_ssp
  !< FOODIE integrator: provide an explicit class of Multi-step Runge-Kutta Methods with Strong Stability Preserving property,
  !< from 3rd to 8th order accurate.
  !<
  !< @note The integrator must be created or initialized (initialize the *A,Ahat,B,Bhat,D,T* coefficients) before used.
  private
  real(R_P), allocatable :: A(:,:)    !< *A* coefficients.
  real(R_P), allocatable :: Ahat(:,:) !< *Ahat* coefficients.
  real(R_P), allocatable :: B(:)      !< *B* coefficients.
  real(R_P), allocatable :: Bhat(:)   !< *Bhat* coefficients.
  real(R_P), allocatable :: D(:,:)    !< *D* coefficients.
  real(R_P), allocatable :: Q(:)      !< *T* coefficients.
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
endtype integrator_ms_runge_kutta_ssp

contains
  ! deferred methods
  pure function class_name(self)
  !< Return the class name of schemes.
  class(integrator_ms_runge_kutta_ssp), intent(in) :: self       !< Integrator.
  character(len=99)                                :: class_name !< Class name.

  class_name = trim(adjustl(class_name_))
  endfunction class_name

  elemental function has_fast_mode(self)
  !< Return .true. if the integrator class has *fast mode* integrate.
  class(integrator_ms_runge_kutta_ssp), intent(in) :: self          !< Integrator.
  logical                                          :: has_fast_mode !< Inquire result.

  has_fast_mode = has_fast_mode_
  endfunction has_fast_mode

  subroutine integr_assign_integr(lhs, rhs)
  !< Operator `=`.
  class(integrator_ms_runge_kutta_ssp), intent(inout) :: lhs !< Left hand side.
  class(integrator_object),             intent(in)    :: rhs !< Right hand side.

  call lhs%assign_multistage_multistep(rhs=rhs)
  select type(rhs)
  class is (integrator_ms_runge_kutta_ssp)
    if (allocated(rhs%A   )) lhs%A    = rhs%A
    if (allocated(rhs%Ahat)) lhs%Ahat = rhs%Ahat
    if (allocated(rhs%B   )) lhs%B    = rhs%B
    if (allocated(rhs%Bhat)) lhs%Bhat = rhs%Bhat
    if (allocated(rhs%D   )) lhs%D    = rhs%D
    if (allocated(rhs%Q   )) lhs%Q    = rhs%Q
  endselect
  endsubroutine integr_assign_integr

  subroutine integrate(self, U, Dt, t)
  !< Integrate field with LMM-SSP class scheme.
  class(integrator_ms_runge_kutta_ssp), intent(inout) :: self  !< Integrator.
  class(integrand_object),              intent(inout) :: U     !< Field to be integrated.
  real(R_P),                            intent(in)    :: Dt    !< Time step.
  real(R_P),                            intent(in)    :: t     !< Time.
  integer(I_P)                                        :: k, kk !< Stages counters.
  integer(I_P)                                        :: s     !< Steps counter.

  ! computing stages
  self%stage(1) = U
  do k=2, self%stages
    self%stage(k) = U * 0._R_P
    do s=1, self%steps
      if (self%D(k, s) /= 0._R_P) self%stage(k) = self%stage(k) + (self%previous(s) * self%D(k, s))
    enddo
    do s=1, self%steps - 1
      if (self%Ahat(k, s) /= 0._R_P) self%stage(k) = self%stage(k) + (self%previous(s)%t(t=self%t(s)) * (Dt * self%Ahat(k, s)))
    enddo
    do kk=1, k - 1
      if (self%A(k, kk) /= 0._R_P) self%stage(k) = self%stage(k) + (self%stage(kk)%t(t=t) * (Dt * self%A(k, kk)))
    enddo
  enddo
  ! computing new time step
  U = U * 0._R_P
  do s=1, self%steps
    if (self%Q(s) /= 0._R_P) U = U + (self%previous(s) * self%Q(s))
  enddo
  do s=1, self%steps - 1
    if (self%Bhat(s) /= 0._R_P) U = U + (self%previous(s)%t(t=self%t(s)) * (Dt * self%Bhat(s)))
  enddo
  do k=1, self%stages
    if (self%B(k) /= 0._R_P) U = U + (self%stage(k)%t(t=t) * (Dt * self%B(k)))
  enddo
  if (self%autoupdate) call self%update_previous(U=U, previous=self%previous, Dt=Dt, t=t, previous_t=self%t)
  endsubroutine integrate

  subroutine integrate_fast(self, U, Dt, t)
  !< Integrate field with LMM-SSP class scheme, fast mode.
  class(integrator_ms_runge_kutta_ssp), intent(inout) :: self  !< Integrator.
  class(integrand_object),              intent(inout) :: U     !< Field to be integrated.
  real(R_P),                            intent(in)    :: Dt    !< Time step.
  real(R_P),                            intent(in)    :: t     !< Time.
  integer(I_P)                                        :: k, kk !< Stages counters.
  integer(I_P)                                        :: s     !< Steps counter.

  ! computing stages
  self%stage(1) = U
  do k=2, self%stages
    call self%stage(k)%multiply_fast(lhs=U, rhs=0._R_P)
    do s=1, self%steps
      if (self%D(k, s) /= 0._R_P) then
         call self%buffer%multiply_fast(lhs=self%previous(s), rhs=self%D(k, s))
         call self%stage(k)%add_fast(lhs=self%stage(k), rhs=self%buffer)
      endif
    enddo
    do s=1, self%steps - 1
      if (self%Ahat(k, s) /= 0._R_P) then
         self%buffer = self%previous(s)
         call self%buffer%t_fast(t=self%t(s))
         call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%Ahat(k, s))
         call self%stage(k)%add_fast(lhs=self%stage(k), rhs=self%buffer)
      endif
    enddo
    do kk=1, k - 1
      if (self%A(k, kk) /= 0._R_P) then
         self%buffer = self%stage(kk)
         call self%buffer%t_fast(t=t)
         call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%A(k, kk))
         call self%stage(k)%add_fast(lhs=self%stage(k), rhs=self%buffer)
      endif
    enddo
  enddo
  ! computing new time step
  U = U * 0._R_P
  do s=1, self%steps
    if (self%Q(s) /= 0._R_P) then
      call self%buffer%multiply_fast(lhs=self%previous(s), rhs=self%Q(s))
      call U%add_fast(lhs=U, rhs=self%buffer)
    endif
  enddo
  do s=1, self%steps - 1
    if (self%Bhat(s) /= 0._R_P) then
       self%buffer = self%previous(s)
       call self%buffer%t_fast(t=self%t(s))
       call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%Bhat(s))
       call U%add_fast(lhs=U, rhs=self%buffer)
    endif
  enddo
  do k=1, self%stages
    if (self%B(k) /= 0._R_P) U = U + (self%stage(k)%t(t=t) * (Dt * self%B(k)))
    if (self%B(k) /= 0._R_P) then
       self%buffer = self%stage(k)
       call self%buffer%t_fast(t=t)
       call self%buffer%multiply_fast(lhs=self%buffer, rhs=Dt * self%B(k))
       call U%add_fast(lhs=U, rhs=self%buffer)
    endif
  enddo
  if (self%autoupdate) call self%update_previous(U=U, previous=self%previous, Dt=Dt, t=t, previous_t=self%t)
  endsubroutine integrate_fast

  elemental function is_supported(self, scheme)
  !< Return .true. if the integrator class support the given scheme.
  class(integrator_ms_runge_kutta_ssp), intent(in) :: self         !< Integrator.
  character(*),                         intent(in) :: scheme       !< Selected scheme.
  logical                                          :: is_supported !< Inquire result.
  integer(I_P)                                     :: s            !< Counter.

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
  class(integrator_ms_runge_kutta_ssp), intent(in) :: self       !< Integrator.
  character(len=99), allocatable                   :: schemes(:) !< Queried scheme.

  allocate(schemes(lbound(supported_schemes_, dim=1):ubound(supported_schemes_, dim=1)))
  schemes = supported_schemes_
  endfunction supported_schemes

  ! public methods
  elemental subroutine destroy(self)
  !< Destroy the integrator.
  class(integrator_ms_runge_kutta_ssp), intent(inout) :: self !< Integrator.

  call self%destroy_multistage_multistep
  if (allocated(self%A   )) deallocate(self%A   )
  if (allocated(self%Ahat)) deallocate(self%Ahat)
  if (allocated(self%B   )) deallocate(self%B   )
  if (allocated(self%Bhat)) deallocate(self%Bhat)
  if (allocated(self%D   )) deallocate(self%D   )
  if (allocated(self%Q   )) deallocate(self%Q   )
  endsubroutine destroy

  subroutine initialize(self, scheme, iterations, autoupdate, U, stop_on_fail)
  !< Create the actual MS-RK-SSP integrator: initialize the *A,Ahat,B,Bhat,D,T* coefficients.
  class(integrator_ms_runge_kutta_ssp), intent(inout)        :: self         !< Integrator.
  character(*),                         intent(in)           :: scheme       !< Selected scheme.
  integer(I_P),                         intent(in), optional :: iterations   !< Implicit iterations.
  logical,                              intent(in), optional :: autoupdate   !< Enable cyclic autoupdate of previous time steps.
  class(integrand_object),              intent(in), optional :: U            !< Integrand molding prototype.
  logical,                              intent(in), optional :: stop_on_fail !< Stop execution if initialization fail.

  if (self%is_supported(scheme=scheme)) then
    call self%destroy
    self%description_ = trim(adjustl(scheme))
    select case(trim(adjustl(scheme)))
    case('ms_runge_kutta_ssp_steps_2_stages_2_order_3')
      self%steps = 2
      self%stages = 2

      allocate(self%A(1:self%stages, 1:self%stages)) ; self%A = 0._R_P
      self%A(2, 1) = 0.910683602522959_R_P

      allocate(self%Ahat(1:self%stages, 1:self%steps)) ; self%Ahat = 0._R_P
      self%Ahat(2, 1) = -1.11985107194706e-19_R_P

      allocate(self%B(1:self%stages)) ; self%B = 0._R_P
      self%B(1) = 0.535898384862245_R_P
      self%B(2) = 0.803847577293368_R_P

      allocate(self%Bhat(1:self%steps)) ; self%Bhat = 0._R_P
      self%Bhat(1) = 0.267949192431123_R_P

      allocate(self%D(1:self%stages, 1:self%steps)) ; self%D = 0._R_P
      self%D(2, 1) = 1._R_P / 3._R_P
      self%D(2, 2) = 2._R_P / 3._R_P

      allocate(self%Q(1:self%steps)) ; self%Q = 0._R_P
      self%Q(1) = 0.607695154586736_R_P
      self%Q(2) = 0.392304845413264_R_P

    case('ms_runge_kutta_ssp_steps_3_stages_2_order_3')
      self%steps = 3
      self%stages = 2

      allocate(self%A(1:self%stages, 1:self%stages)) ; self%A = 0._R_P
      self%A(2, 1) = 0.731058363135786_R_P

      allocate(self%Ahat(1:self%stages, 1:self%steps)) ; self%Ahat = 0._R_P
      self%Ahat(2, 1) = 0.127467809251820_R_P

      allocate(self%B(1:self%stages)) ; self%B = 0._R_P
      self%B(1) = 0.618048297723782_R_P
      self%B(2) = 0.759677988437936_R_P

      allocate(self%Bhat(1:self%steps)) ; self%Bhat = 0._R_P
      self%Bhat(1) = 0.246670340394148_R_P

      allocate(self%D(1:self%stages, 1:self%steps)) ; self%D = 0._R_P
      self%D(2, 1) = 0.186433848116852_R_P
      self%D(2, 2) = 1.80945758995975e-24_R_P
      self%D(2, 3) = 0.813566151883148_R_P

      allocate(self%Q(1:self%steps)) ; self%Q = 0._R_P
      self%Q(1) = 0.312198313277933_R_P
      self%Q(2) = 2.58493941422821e-24_R_P
      self%Q(3) = 0.687801686722067_R_P

    case('ms_runge_kutta_ssp_steps_4_stages_5_order_8')
      self%steps = 4
      self%stages = 5

      allocate(self%A(1:self%stages, 1:self%stages)) ; self%A = 0._R_P
      self%A(2, 1) = 0.0112355687952080_R_P
      self%A(3, 1) = 0.782384118905967_R_P
      self%A(4, 1) = 0.0997788285846345_R_P
      self%A(5, 1) = 0.0173871875042219_R_P

      self%A(4, 2) = 0.775239818309315_R_P
      self%A(5, 2) = 0.781995253645396_R_P

      self%A(4, 3) = 0.522304633131092_R_P
      self%A(5, 3) = 0.0817254029032851_R_P

      self%A(5, 4) = 0.654483113500859_R_P

      allocate(self%Ahat(1:self%stages, 1:self%steps)) ; self%Ahat = 0._R_P
      self%Ahat(2, 1) = 0.000422703250336666_R_P
      self%Ahat(3, 1) = 0.0959783036454617_R_P
      self%Ahat(4, 1) = 0.0140562464699573_R_P
      self%Ahat(5, 1) = 0.0519851819388547_R_P

      self%Ahat(2, 2) = 0.259546324808661_R_P
      self%Ahat(3, 2) = 0.382496291927802_R_P
      self%Ahat(4, 2) = 0.195323228972419_R_P
      self%Ahat(5, 2) = 0.435648262830826_R_P

      self%Ahat(2, 3) = 0.752684925098657_R_P
      self%Ahat(3, 3) = 0.563081036068107_R_P
      self%Ahat(4, 3) = 0.209815168854422_R_P
      self%Ahat(5, 3) = 0.151720560507208_R_P

      allocate(self%B(1:self%stages)) ; self%B = 0._R_P
      self%B(1) = 0.711472565648602_R_P
      self%B(2) = 0.0953138922500395_R_P
      self%B(3) = 0.0808915576045876_R_P
      self%B(4) = 0.214580109044146_R_P
      self%B(5) = 0.351640244526174_R_P

      allocate(self%Bhat(1:self%steps)) ; self%Bhat = 0._R_P
      self%Bhat(1) = 0.00614782373612238_R_P
      self%Bhat(2) = 0.138226341918060_R_P
      self%Bhat(3) = 0.541410372692778_R_P

      allocate(self%D(1:self%stages, 1:self%steps)) ; self%D = 0._R_P
      self%D(2, 1) = 0.0273928990974108_R_P
      self%D(3, 1) = 0.283607987548794_R_P
      self%D(4, 1) = 0.0642241421937960_R_P
      self%D(5, 1) = 0.194681462814288_R_P

      self%D(2, 2) = 0.554039201229659_R_P
      self%D(3, 2) = 0.0914454235177934_R_P
      self%D(4, 2) = 0.198601843371796_R_P
      self%D(5, 2) = 0.293086617882372_R_P

      self%D(2, 3) = 0.348927848402249_R_P
      self%D(3, 3) = 0.437897855084625_R_P
      self%D(4, 3) = 0.662266498804591_R_P
      self%D(5, 3) = 0.158740367819382_R_P

      self%D(2, 4) = 0.0696400512706807_R_P
      self%D(3, 4) = 0.187048733848788_R_P
      self%D(4, 4) = 0.0749075156298171_R_P
      self%D(5, 4) = 0.353491551483958_R_P

      allocate(self%Q(1:self%steps)) ; self%Q = 0._R_P
      self%Q(1) = 0.0273988434707855_R_P
      self%Q(2) = 0.286296288278021_R_P
      self%Q(3) = 0.484893800452111_R_P
      self%Q(4) = 0.201411067799082_R_P
    endselect
    self%autoupdate = .true. ; if (present(autoupdate)) self%autoupdate = autoupdate
    self%iterations = 1 ; if (present(iterations)) self%iterations = iterations
    self%registers_stages = self%stages
    self%registers_steps = self%steps
    if (present(U)) call self%allocate_integrand_members(U=U)
  else
    call self%trigger_error(error=ERROR_UNSUPPORTED_SCHEME,                                   &
                            error_message='"'//trim(adjustl(scheme))//'" unsupported scheme', &
                            is_severe=stop_on_fail)
  endif
  endsubroutine initialize
endmodule foodie_integrator_ms_runge_kutta_ssp
