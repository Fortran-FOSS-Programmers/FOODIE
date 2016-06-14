!< Test WenOOF with the reconstruction of sin function.
program sin_reconstruction
!-----------------------------------------------------------------------------------------------------------------------------------
!< Test WenOOF with the reconstruction of sin function.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P, str
use wenoof, only : weno_factory, weno_constructor_upwind, weno_interpolator
use pyplot_module, only :  pyplot
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(weno_factory)                    :: factory                    !< WENO factory.
class(weno_interpolator), allocatable :: interpolator               !< WENO interpolator.
integer(I_P), parameter               :: S = 3_I_P                  !< Stencils used.
integer(I_P), parameter               :: Nv = 30_I_P                !< Number of discretized values to be interpolated.
real(R_P),    parameter               :: pi = 4._R_P * atan(1._R_P) !< Extent of domain.
real(R_P)                             :: x(1-S:Nv+S)                !< Whole domain.
real(R_P)                             :: fx(1-S:Nv+S)               !< Discretized values to be interpolated.
real(R_P)                             :: xi(1:Nv)                   !< Domain of the interpolation.
real(R_P)                             :: fx_ref(1:Nv)               !< Reference values.
real(R_P)                             :: interpolation(1:1, 1:Nv)   !< Interpolated values.
type(pyplot)                          :: plt                        !< Plotter handler.
integer                               :: i, j, f                    !< Counters.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! build the values used for the reconstruction of sin function: nodal values
x = 0.
do i = 1 - S, Nv + S
  x(i) = i * 2 * pi / Nv
  fx(i) = sin(x(i))
enddo
! face values to which the reconstruction should tend
do i = 1, Nv
  xi(i) = x(i) + pi / Nv
  fx_ref(i) = sin(xi(i))
enddo

! prepare the weno interpolator
call factory%create(constructor=weno_constructor_upwind(S=S, eps=10._R_P**(-40)), interpolator=interpolator)

! interpolate values
interpolation = 0.
do i = 1, Nv ! interpolated values loop
  call interpolator%interpolate(S=S,                                                      &
                                stencil=reshape(source=fx(i+1-S:i-1+S), shape=[1,2*S-1]), &
                                location='right',                                         &
                                interpolation=interpolation(1:1, i))
enddo

! print results
print "(A)", '# x, sin(x), weno_interpolation(x)'
do i = 1, Nv
  print "(A)", str(n=xi(i))//' '//str(n=fx_ref(i))//' '//str(n=interpolation(1, i))
enddo

! plotting graph to image file
call plt%initialize(grid=.true., xlabel='angle (rad)', title='WENO interpolation of $\sin(x)$', legend=.true.)
call plt%add_plot(x=xi(1:Nv),        &
                  y=fx_ref(1:Nv),    &
                  label='$\sin(x)$', &
                  linestyle='k-',    &
                  linewidth=2)
call plt%add_plot(x=xi(1:Nv),                 &
                  y=interpolation(1, 1:Nv),   &
                  label='WENO interpolation', &
                  linestyle='ro',             &
                  markersize=6)
call plt%savefig('sin_reconstruction.png')
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram sin_reconstruction
