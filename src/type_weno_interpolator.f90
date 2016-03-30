module type_weno_interpolator
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO interpolator object,
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_interpolator, weno_constructor
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_constructor
  !< Abstract type used for create new concrete WENO interpolators.
  !<
  !< @note Every concrete WENO interpolator implementations must define their own constructor type.
  private
endtype weno_constructor

type, abstract :: weno_interpolator
  !< WENO interpolator object.
  !<
  !< @note Do not implement any real interpolator: provide the interface for the different interpolators implemented.
  private
  contains
    procedure(abstract_destructor),  pass(self), deferred, public :: destroy
    procedure(abstract_constructor), pass(self), deferred, public :: create
    procedure(abstract_description), pass(self), deferred, public :: description
    procedure(abstract_interpolate), pass(self), deferred, public :: interpolate
endtype weno_interpolator

abstract interface
  elemental subroutine abstract_destructor(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator
  class(weno_interpolator), intent(INOUT) :: self !< WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_destructor

  subroutine abstract_constructor(self, constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create a WENO interpolator.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator
  import :: weno_constructor
  class(weno_interpolator), intent(INOUT) :: self        !< WENO interpolator.
  class(weno_constructor),  intent(IN)    :: constructor !< WENO constructor.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_constructor

  pure subroutine abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator
  class(weno_interpolator),      intent(IN)  :: self   !< WENO interpolator.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_description

  pure subroutine abstract_interpolate(self, S, stencil, location, interpolation)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Interpolate the stecil input values computing the actual interpolation.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator, I_P, R_P
  class(weno_interpolator), intent(IN)  :: self                !< WENO interpolator.
  integer(I_P),             intent(IN)  :: S                   !< Number of stencils used.
  real(R_P),                intent(IN)  :: stencil(1:, 1 - S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  character(*),             intent(IN)  :: location            !< Location of interpolated value(s): central, left, right, both.
  real(R_P),                intent(OUT) :: interpolation(1:)   !< Result of the interpolation, [1:2].
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_interpolate
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_interpolator
