module type_weno_interpolator_upwind
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing upwind biased WENO interpolator object and constructor,
!<
!< @note The provided WENO interpolator implements the *Efficient Implementation of Weighted ENO Schemes*,
!< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : I_P, R_P, str
use type_weno_interpolator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_interpolator_upwind, weno_constructor_upwind
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(weno_constructor) :: weno_constructor_upwind
  !< Upwind biased WENO interpolator constructor,
  !<
  !< @note The constructed WENO interpolator implements the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  integer(I_P) :: S = 0               !< Stencils dimension.
  real(R_P)    :: eps = 10._R_P**(-6) !< Parameter for avoiding divided by zero when computing smoothness indicators.
endtype weno_constructor_upwind
interface weno_constructor_upwind
  procedure weno_constructor_upwind_init
endinterface

type, extends(weno_interpolator) :: weno_interpolator_upwind
  !< Upwind biased WENO interpolator object,
  !<
  !< @note The WENO interpolator implemented is the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  !<
  !< @note The supported accuracy formal order are: 3rd, 5th, 7th corresponding to use 2, 3, 4 stencils composed of 2, 3, 4 values,
  !< respectively.
  private
  integer(I_P)           :: S = 0_I_P          !< Stencil dimension.
  real(R_P)              :: eps = 0._R_P       !< Parameter for avoiding divided by zero when computing smoothness indicators.
  real(R_P), allocatable :: weights_opt(:,:)   !< Optimal weights                    [1:2,0:S-1].
  real(R_P), allocatable :: poly_coef(:,:,:)   !< Polynomials coefficients           [1:2,0:S-1,0:S-1].
  real(R_P), allocatable :: smooth_coef(:,:,:) !< Smoothness indicators coefficients [0:S-1,0:S-1,0:S-1].
  contains
    ! public methods
    procedure, pass(self), public :: destroy
    procedure, pass(self), public :: create
    procedure, pass(self), public :: description
    procedure, pass(self), public :: interpolate
    generic,               public :: assignment(=) => assign_interpolator !< Overloading = assignament.
    ! private methods
    procedure, pass(lhs), private :: assign_interpolator !< Assignament operator.
    final :: finalize
endtype weno_interpolator_upwind
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! weno_constructor_upwind
  elemental function  weno_constructor_upwind_init(S, eps) result(constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create (initialize) the WENO interpolator.
  !<
  !< @note For this class of interpolators it is sufficient to provide the maximum number of stencils used (that is also the
  !< dimension, i.e. number of values, of each stencil). During the actual interpolation phase the client code can specify, for each
  !< intepolation a different number of stencil bounded by this maximum value. This is useful for coupling the interpolator with
  !< algorithm like the Recursive Order Reduction (ROR) strategy.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)           :: S           !< Maximum stencils dimension.
  real(R_P),    intent(IN), optional :: eps         !< Parameter for avoiding divided by zero when computing smoothness indicators.
  type(weno_constructor_upwind)      :: constructor !<WENO constructor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  constructor%S = S
  if (present(eps)) constructor%eps = eps
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction  weno_constructor_upwind_init

  ! weno_interpolator_upwind
  ! public methods
  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy the WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(INOUT) :: self !< WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%S = 0_I_P
  self%eps = 0._R_P
  if (allocated(self%weights_opt)) deallocate(self%weights_opt)
  if (allocated(self%poly_coef  )) deallocate(self%poly_coef  )
  if (allocated(self%smooth_coef)) deallocate(self%smooth_coef)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine create(self, constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(INOUT) :: self        !< WENO interpolator.
  class(weno_constructor),         intent(IN)    :: constructor !< WENO constructor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(constructor)
  type is(weno_constructor_upwind)
    call self%destroy
    self%S = constructor%S
    self%eps = constructor%eps
    allocate(self%weights_opt(1:2, 0:self%S - 1))
    allocate(self%poly_coef(1:2, 0:self%S - 1, 0:self%S - 1))
    allocate(self%smooth_coef(0:self%S - 1, 0:self%S - 1, 0:self%S - 1))
    call set_weights_optimal
    call set_polynomial_coefficients
    call set_smoothness_indicators_coefficients
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    subroutine set_weights_optimal()
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Set the values of optimial weights.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    select case(self%S)
    case(2) ! 3rd order
      ! 1 => left interface (i-1/2)
      self%weights_opt(1, 0) = 2._R_P/3._R_P ! stencil 0
      self%weights_opt(1, 1) = 1._R_P/3._R_P ! stencil 1
      ! 2 => right interface (i+1/2)
      self%weights_opt(2, 0) = 1._R_P/3._R_P ! stencil 0
      self%weights_opt(2, 1) = 2._R_P/3._R_P ! stencil 1
    case(3) ! 5th order
      ! 1 => left interface (i-1/2)
      self%weights_opt(1, 0) = 0.3_R_P ! stencil 0
      self%weights_opt(1, 1) = 0.6_R_P ! stencil 1
      self%weights_opt(1, 2) = 0.1_R_P ! stencil 2
      ! 2 => right interface (i+1/2)
      self%weights_opt(2, 0) = 0.1_R_P ! stencil 0
      self%weights_opt(2, 1) = 0.6_R_P ! stencil 1
      self%weights_opt(2, 2) = 0.3_R_P ! stencil 2
    case(4) ! 7th order
      ! 1 => left interface (i-1/2)
      self%weights_opt(1, 0) =  4._R_P/35._R_P ! stencil 0
      self%weights_opt(1, 1) = 18._R_P/35._R_P ! stencil 1
      self%weights_opt(1, 2) = 12._R_P/35._R_P ! stencil 2
      self%weights_opt(1, 3) =  1._R_P/35._R_P ! stencil 3
      ! 2 => right interface (i+1/2)
      self%weights_opt(2, 0) =  1._R_P/35._R_P ! stencil 0
      self%weights_opt(2, 1) = 12._R_P/35._R_P ! stencil 1
      self%weights_opt(2, 2) = 18._R_P/35._R_P ! stencil 2
      self%weights_opt(2, 3) =  4._R_P/35._R_P ! stencil 3
    endselect
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_weights_optimal

    subroutine set_polynomial_coefficients()
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Set the values of polynomial_coefficient.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    associate(coef => self%poly_coef)
      select case(self%S)
      case(2) ! 3rd order
        ! 1 => left interface (i-1/2)
        !  cell  0               ;    cell  1
        coef(1, 0, 0) =  0.5_R_P ; coef(1, 1, 0) =  0.5_R_P ! stencil 0
        coef(1, 0, 1) = -0.5_R_P ; coef(1, 1, 1) =  1.5_R_P ! stencil 1
        ! 2 => right interface (i+1/2)
        !  cell  0               ;    cell  1
        coef(2, 0, 0) =  1.5_R_P ; coef(2, 1, 0) = -0.5_R_P ! stencil 0
        coef(2, 0, 1) =  0.5_R_P ; coef(2, 1, 1) =  0.5_R_P ! stencil 1
      case(3) ! 5th order
        ! 1 => left interface (i-1/2)
        !  cell  0                     ;    cell  1                     ;    cell  2
        coef(1, 0, 0) =  1._R_P/3._R_P ; coef(1, 1, 0) =  5._R_P/6._R_P ; coef(1, 2, 0) = -1._R_P/6._R_P ! stencil 0
        coef(1, 0, 1) = -1._R_P/6._R_P ; coef(1, 1, 1) =  5._R_P/6._R_P ; coef(1, 2, 1) =  1._R_P/3._R_P ! stencil 1
        coef(1, 0, 2) =  1._R_P/3._R_P ; coef(1, 1, 2) = -7._R_P/6._R_P ; coef(1, 2, 2) = 11._R_P/6._R_P ! stencil 2
        ! 2 => right interface (i+1/2)
        !  cell  0                     ;    cell  1                     ;    cell  2
        coef(2, 0, 0) = 11._R_P/6._R_P ; coef(2, 1, 0) = -7._R_P/6._R_P ; coef(2, 2, 0) =  1._R_P/3._R_P ! stencil 0
        coef(2, 0, 1) =  1._R_P/3._R_P ; coef(2, 1, 1) =  5._R_P/6._R_P ; coef(2, 2, 1) = -1._R_P/6._R_P ! stencil 1
        coef(2, 0, 2) = -1._R_P/6._R_P ; coef(2, 1, 2) =  5._R_P/6._R_P ; coef(2, 2, 2) =  1._R_P/3._R_P ! stencil 2
      case(4) ! 7th order
        ! 1 => left interface (i-1/2)
        !  cell  0                 ;   cell  1                 ;   cell  2                  ;   cell  3
        coef(1,0,0)= 1._R_P/4._R_P ;coef(1,1,0)=13._R_P/12._R_P;coef(1,2,0)= -5._R_P/12._R_P;coef(1,3,0)= 1._R_P/12._R_P ! sten 0
        coef(1,0,1)=-1._R_P/12._R_P;coef(1,1,1)= 7._R_P/12._R_P;coef(1,2,1)=  7._R_P/12._R_P;coef(1,3,1)=-1._R_P/12._R_P ! sten 1
        coef(1,0,2)= 1._R_P/12._R_P;coef(1,1,2)=-5._R_P/12._R_P;coef(1,2,2)= 13._R_P/12._R_P;coef(1,3,2)= 1._R_P/4._R_P  ! sten 2
        coef(1,0,3)=-1._R_P/4._R_P ;coef(1,1,3)=13._R_P/12._R_P;coef(1,2,3)=-23._R_P/12._R_P;coef(1,3,3)=25._R_P/12._R_P ! sten 3
        ! 2 => right interface (i+1/2)
        !  cell  0                 ;   cell  1                  ;   cell  2                 ;   cell  3
        coef(2,0,0)=25._R_P/12._R_P;coef(2,1,0)=-23._R_P/12._R_P;coef(2,2,0)=13._R_P/12._R_P;coef(2,3,0)=-1._R_P/4._R_P  ! sten 0
        coef(2,0,1)= 1._R_P/4._R_P ;coef(2,1,1)= 13._R_P/12._R_P;coef(2,2,1)=-5._R_P/12._R_P;coef(2,3,1)= 1._R_P/12._R_P ! sten 1
        coef(2,0,2)=-1._R_P/12._R_P;coef(2,1,2)=  7._R_P/12._R_P;coef(2,2,2)= 7._R_P/12._R_P;coef(2,3,2)=-1._R_P/12._R_P ! sten 2
        coef(2,0,3)= 1._R_P/12._R_P;coef(2,1,3)= -5._R_P/12._R_P;coef(2,2,3)=13._R_P/12._R_P;coef(2,3,3)= 1._R_P/4._R_P  ! sten 3
      endselect
    endassociate
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_polynomial_coefficients

    subroutine set_smoothness_indicators_coefficients()
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Set the values of smoothness indicators coefficients.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    associate(coef => self%smooth_coef)
      select case(self%S)
      case(2) ! 3rd order
        ! stencil 0
        !      i*i             ;       (i-1)*i
        coef(0, 0, 0) = 1._R_P ; coef(1, 0, 0) = -2._R_P
        !      /               ;       (i-1)*(i-1)
        coef(0, 1, 0) = 0._R_P ; coef(1, 1, 0) = 1._R_P
        ! stencil 1
        !     (i+1)*(i+1)      ;       (i+1)*i
        coef(0, 0, 1) = 1._R_P ; coef(1, 0, 1) = -2._R_P
        !      /               ;        i*i
        coef(0, 1, 1) = 0._R_P ; coef(1, 1, 1) = 1._R_P
      case(3) ! 5th order
        ! stencil 0
        !      i*i                      ;       (i-1)*i                   ;       (i-2)*i
        coef(0, 0, 0) =  10._R_P/3._R_P ; coef(1, 0, 0) = -31._R_P/3._R_P ; coef(2, 0, 0) =  11._R_P/3._R_P
        !      /                        ;       (i-1)*(i-1)               ;       (i-2)*(i-1)
        coef(0, 1, 0) =   0._R_P        ; coef(1, 1, 0) =  25._R_P/3._R_P ; coef(2, 1, 0) = -19._R_P/3._R_P
        !      /                        ;        /                        ;       (i-2)*(i-2)
        coef(0, 2, 0) =   0._R_P        ; coef(1, 2, 0) =   0._R_P        ; coef(2, 2, 0) =   4._R_P/3._R_P
        ! stencil 1
        !     (i+1)*(i+1)               ;        i*(i+1)                  ;       (i-1)*(i+1)
        coef(0, 0, 1) =   4._R_P/3._R_P ; coef(1, 0, 1) = -13._R_P/3._R_P ; coef(2, 0, 1) =   5._R_P/3._R_P
        !      /                        ;        i*i                      ;       (i-1)*i
        coef(0, 1, 1) =   0._R_P        ; coef(1, 1, 1) =  13._R_P/3._R_P ; coef(2, 1, 1) = -13._R_P/3._R_P
        !      /                        ;        /                        ;       (i-1)*(i-1)
        coef(0, 2, 1) =   0._R_P        ; coef(1, 2, 1) =   0._R_P        ; coef(2, 2, 1) =   4._R_P/3._R_P
        ! stencil 2
        !     (i+2)*(i+2)               ;       (i+1)*(i+2)               ;        i*(i+2)
        coef(0, 0, 2) =   4._R_P/3._R_P ; coef(1, 0, 2) = -19._R_P/3._R_P ; coef(2, 0, 2) =  11._R_P/3._R_P
        !      /                        ;       (i+1)*(i+1)               ;        i*(i+1)
        coef(0, 1, 2) =   0._R_P        ; coef(1, 1, 2) =  25._R_P/3._R_P ; coef(2, 1, 2) = -31._R_P/3._R_P
        !      /                        ;        /                        ;        i*i
        coef(0, 2, 2) =   0._R_P        ; coef(1, 2, 2) =   0._R_P        ; coef(2, 2, 2) =  10._R_P/3._R_P
      case(4) ! 7th order
        ! stencil 0
        !      i*i                ;       (i-1)*i             ;       (i-2)*i              ;       (i-3)*i
        coef(0, 0, 0) = 2107._R_P ; coef(1, 0, 0) =-9402._R_P ; coef(2, 0, 0) =  7042._R_P ; coef(3, 0, 0) = -1854._R_P
        !      /                  ;       (i-1)*(i-1)         ;       (i-2)*(i-1)          ;       (i-3)*(i-1)
        coef(0, 1, 0) =    0._R_P ; coef(1, 1, 0) =11003._R_P ; coef(2, 1, 0) =-17246._R_P ; coef(3, 1, 0) =  4642._R_P
        !      /                  ;        /                  ;       (i-2)*(i-2)          ;       (i-3)*(i-2)
        coef(0, 2, 0) =    0._R_P ; coef(1, 2, 0) =    0._R_P ; coef(2, 2, 0) =  7043._R_P ; coef(3, 2, 0) = -3882._R_P
        !      /                  ;        /                  ;        /                   ;       (i-3)*(i-3)
        coef(0, 3, 0) =    0._R_P ; coef(1, 3, 0) =    0._R_P ; coef(2, 3, 0) =     0._R_P ; coef(3, 3, 0) =   547._R_P
        ! stencil 1
        !     (i+1)*(i+1)         ;        i*(i+1)            ;       (i-1)*(i+1)          ;       (i-2)*(i+1)
        coef(0, 0, 1) =  547._R_P ; coef(1, 0, 1) =-2522._R_P ; coef(2, 0, 1) =  1922._R_P ; coef(3, 0, 1) =  -494._R_P
        !       /                 ;          i*i              ;       (i-1)*i              ;       (i-2)*i
        coef(0, 1, 1) =    0._R_P ; coef(1, 1, 1) = 3443._R_P ; coef(2, 1, 1) = -5966._R_P ; coef(3, 1, 1) =  1602._R_P
        !       /                 ;          /                ;       (i-1)*(i-1)          ;       (i-2)*(i-1)
        coef(0, 2, 1) =    0._R_P ; coef(1, 2, 1) =    0._R_P ; coef(2, 2, 1) =  2843._R_P ; coef(3, 2, 1) = -1642._R_P
        !       /                 ;          /                ;        /                   ;       (i-2)*(i-2)
        coef(0, 3, 1) =    0._R_P ; coef(1, 3, 1) =    0._R_P ; coef(2, 3, 1) =     0._R_P ; coef(3, 3, 1) =   267._R_P
        ! stencil 2
        !     (i+2)*(i+2)         ;       (i+1)*(i+2)         ;            i*(i+2)         ;       (i-1)*(i+2)
        coef(0, 0, 2) =  267._R_P ; coef(1, 0, 2) =-1642._R_P ; coef(2, 0, 2) =  1602._R_P ; coef(3, 0, 2) =  -494._R_P
        !      /                  ;       (i+1)*(i+1)         ;        i*(i+1)             ;       (i-1)*(i+1)
        coef(0, 1, 2) =    0._R_P ; coef(1, 1, 2) = 2843._R_P ; coef(2, 1, 2) = -5966._R_P ; coef(3, 1, 2) =  1922._R_P
        !      /                  ;        /                  ;        i*i                 ;       (i-1)*i
        coef(0, 2, 2) =    0._R_P ; coef(1, 2, 2) =    0._R_P ; coef(2, 2, 2) =  3443._R_P ; coef(3, 2, 2) = -2522._R_P
        !      /                  ;        /                  ;        /                   ;       (i-1)*(i-1)
        coef(0, 3, 2) =    0._R_P ; coef(1, 3, 2) =    0._R_P ; coef(2, 3, 2) =     0._R_P ; coef(3, 3, 2) =   547._R_P
        ! stencil 3
        !     (i+3)*(i+3)         ;       (i+2)*(i+3)         ;           (i+1)*(i+3)      ;        i*(i+3)
        coef(0, 0, 3) =  547._R_P ; coef(1, 0, 3) =-3882._R_P ; coef(2, 0, 3) =  4642._R_P ; coef(3, 0, 3) = -1854._R_P
        !      /                  ;       (i+2)*(i+2)         ;       (i+1)*(i+2)          ;        i*(i+2)
        coef(0, 1, 3) =    0._R_P ; coef(1, 1, 3) = 7043._R_P ; coef(2, 1, 3) =-17246._R_P ; coef(3, 1, 3) =  7042._R_P
        !      /                  ;        /                  ;       (i+1)*(i+1)          ;        i*(i+1)
        coef(0, 2, 3) =    0._R_P ; coef(1, 2, 3) =    0._R_P ; coef(2, 2, 3) = 11003._R_P ; coef(3, 2, 3) = -9402._R_P
        !      /                  ;        /                  ;        /                   ;        i*i
        coef(0, 3, 3) =    0._R_P ; coef(1, 3, 3) =    0._R_P ; coef(2, 3, 3) =     0._R_P ; coef(3, 3, 3) =  2107._R_P
      endselect
    endassociate
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_smoothness_indicators_coefficients
  endsubroutine create

  pure subroutine description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing the WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(IN)  :: self             !< WENO interpolator.
  character(len=:), allocatable,   intent(OUT) :: string           !< String returned.
  character(len=1)                             :: dummy            !< Dummy string.
  character(len=1), parameter                  :: nl=new_line('a') !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO upwind-biased interpolator'//nl
  string = string//'  Based on the scheme proposed by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130'//nl
  string = string//'  Provide a formal order of accuracy equals to: '//trim(str(.true.,2*self%S - 1))//nl
  string = string//'  Use '//trim(str(.true.,self%S))//' stencils composed by '//trim(str(.true.,self%S))//' values'//nl
  string = string//'  The eps value used for avoiding division by zero is '//trim(str(.true.,self%eps))//nl
  string = string//'  The "interpolate" method has the following public API'//nl
  string = string//'    interpolate(S, stencil, location, interpolation)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(IN), the number of stencils actually used'//nl
  string = string//'    stencil(1:, 1-S:-1+S): real(R_P), intent(IN), the stencils used'//nl
  string = string//'    location: character(*), intent(IN), the location of interpolation {left, right, both}'//nl
  string = string//'    interpolation(1:, 1-S:-1+S): realR_P, intent(OUT), the interpolated values'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description

  pure subroutine interpolate(self, S, stencil, location, interpolation)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Interpolate the stecil input values computing the actual interpolation.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(IN)  :: self                      !< WENO interpolator.
  integer,                         intent(IN)  :: S                         !< Number of stencils actually used.
  real(R_P),                       intent(IN)  :: stencil(1:, 1 - S:)       !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  character(*),                    intent(IN)  :: location                  !< Location of interpolated value(s): left, right, both.
  real(R_P),                       intent(OUT) :: interpolation(1:)         !< Result of the interpolation, [1:2].
  real(R_P)                                    :: polynomials(1:2, 0:S - 1) !< Polynomial reconstructions.
  real(R_P)                                    :: weights(1:2, 0:S - 1)     !< Weights of the stencils.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(location)
  case('both', 'b')
    call compute_polynomials(f1=1_I_P, f2=2_I_P, ff=0_I_P, polynomials=polynomials)
    call compute_weights(f1=1_I_P, f2=2_I_P, ff=0_I_P, weights=weights)
    call compute_convolution(f1=1_I_P, f2=2_I_P, ff=0_I_P, interpolation=interpolation)
  case('left', 'l')
    call compute_polynomials(f1=1_I_P, f2=1_I_P, ff=0_I_P, polynomials=polynomials)
    call compute_weights(f1=1_I_P, f2=1_I_P, ff=0_I_P, weights=weights)
    call compute_convolution(f1=1_I_P, f2=1_I_P, ff=0_I_P, interpolation=interpolation)
  case('right', 'r')
    call compute_polynomials(f1=2_I_P, f2=2_I_P, ff=-1_I_P, polynomials=polynomials)
    call compute_weights(f1=2_I_P, f2=2_I_P, ff=-1_I_P, weights=weights)
    call compute_convolution(f1=2_I_P, f2=2_I_P, ff=-1_I_P, interpolation=interpolation)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure subroutine compute_polynomials(f1, f2, ff, polynomials)
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Compute the polynomials reconstructions.
    !-------------------------------------------------------------------------------------------------------------------------------
    integer(I_P), intent(IN)  :: f1, f2, ff          !< Faces to be computed.
    real(R_P),    intent(OUT) :: polynomials(1:, 0:) !< Polynomial reconstructions.
    integer(I_P)              :: s1, s2, f           !< Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    polynomials = 0.
    do s1 = 0, S - 1 ! stencils loop
      do s2 = 0, S - 1 ! values loop
        do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          polynomials(f, s1) = polynomials(f, s1) + self%poly_coef(f, s2, s1) * stencil(f + ff, -s2 + s1)
        enddo
      enddo
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine compute_polynomials

    pure subroutine compute_weights(f1, f2, ff, weights)
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Compute the stencils weights.
    !-------------------------------------------------------------------------------------------------------------------------------
    integer,   intent(IN)  :: f1, f2, ff       !< Faces to be computed.
    real(R_P), intent(OUT) :: weights(1:, 0:)  !< Weights of the stencils, [1:2, 0:S - 1 ].
    real(R_P)              :: IS(1:2, 0:S - 1) !< Smoothness indicators of the stencils.
    real(R_P)              :: a(1:2, 0:S - 1)  !< Alpha coefficients for the weights.
    real(R_P)              :: a_tot(1:2)       !< Sum of the alpha coefficients.
    integer(I_P)           :: s1, s2, s3, f    !< Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    ! computing smoothness indicators
    do s1 = 0, S - 1 ! stencils loop
      do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        IS(f, s1) = 0.
        do s2 = 0, S - 1
          do s3 = 0, S - 1
            IS(f, s1) = IS(f, s1) + self%smooth_coef(s3, s2, s1) * stencil(f + ff, s1 - s3) * stencil(f + ff, s1 - s2)
          enddo
        enddo
      enddo
    enddo
    ! computing alfa coefficients
    a_tot = 0.
    do s1 = 0, S - 1 ! stencil loops
      do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        a(f, s1) = self%weights_opt(f, s1) * (1./(self%eps + IS(f, s1))**S) ; a_tot(f) = a_tot(f) + a(f, s1)
      enddo
    enddo
    ! computing the weights
    do s1 = 0, S - 1 ! stencils loop
      do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        weights(f, s1) = a(f, s1) / a_tot(f)
      enddo
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine compute_weights

    pure subroutine compute_convolution(f1, f2, ff, interpolation)
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Compute the polynomials convolution.
    !-------------------------------------------------------------------------------------------------------------------------------
    integer(I_P), intent(IN)  :: f1, f2, ff        !< Faces to be computed.
    real(R_P),    intent(OUT) :: interpolation(1:) !< Left and right (1,2) interface value of reconstructed.
    integer(I_P)              :: k, f              !< Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    ! computing the convultion
    interpolation = 0.
    do k = 0, S - 1 ! stencils loop
      do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        interpolation(f + ff) = interpolation(f + ff) + weights(f, k) * polynomials(f, k)
      enddo
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine compute_convolution
  endsubroutine interpolate

  ! private methods
  pure subroutine assign_interpolator(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assign one interpolator to another.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(INOUT) :: lhs !< Left hand side.
  class(weno_interpolator),        intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(rhs)
  class is(weno_interpolator_upwind)
    lhs%S = rhs%S
    lhs%eps = rhs%eps
    if (allocated(rhs%weights_opt)) then
      if (allocated(lhs%weights_opt)) deallocate(lhs%weights_opt) ; allocate(lhs%weights_opt(1:2, 0:lhs%S - 1))
      lhs%weights_opt = rhs%weights_opt
    endif
    if (allocated(rhs%poly_coef)) then
      if (allocated(lhs%poly_coef)) deallocate(lhs%poly_coef) ; allocate(lhs%poly_coef(1:2, 0:lhs%S - 1, 0:lhs%S - 1))
      lhs%poly_coef = rhs%poly_coef
    endif
    if (allocated(rhs%smooth_coef)) then
      if (allocated(lhs%smooth_coef)) deallocate(lhs%smooth_coef) ; allocate(lhs%smooth_coef(0:lhs%S - 1, 0:lhs%S - 1, 0:lhs%S - 1))
      lhs%smooth_coef = rhs%smooth_coef
    endif
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_interpolator

  elemental subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Finalize object.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(weno_interpolator_upwind), intent(INOUT) :: self !< WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule type_weno_interpolator_upwind
