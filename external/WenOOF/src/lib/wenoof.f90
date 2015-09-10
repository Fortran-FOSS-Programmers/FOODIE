module wenoof
!-----------------------------------------------------------------------------------------------------------------------------------
!< WenOOF, WENO interpolation Object Oriented Fortran library
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use type_weno_interpolator, only : weno_constructor, weno_interpolator
use type_weno_interpolator_upwind, only : weno_constructor_upwind, weno_interpolator_upwind
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_factory, weno_constructor, weno_interpolator
public :: weno_constructor_upwind
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: weno_factory
  !< WENO factory aimed to create and return a concrete WENO interpolator to the client code without exposing the concrete
  !< interpolators classes.
  contains
    procedure, nopass :: create
endtype
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine create(constructor, interpolator)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create and return a concrete WENO interpolator object being an extension of the abstract *weno_interpolator* type.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_constructor),               intent(IN)  :: constructor  !< The concrete WENO constructor selected by client code.
  class(weno_interpolator), allocatable, intent(OUT) :: interpolator !< The concrete WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(constructor)
  type is(weno_constructor_upwind)
    allocate(weno_interpolator_upwind :: interpolator)
    call interpolator%create(constructor=constructor)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create
endmodule wenoof
