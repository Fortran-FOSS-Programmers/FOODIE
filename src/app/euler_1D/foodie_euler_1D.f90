!< Euler 1D PDE solver based on FOODIE library.
program foodie_euler_1D
!< Euler 1D PDE solver based on FOODIE library.

use foodie_integrand_euler_1D, only : integrand_euler_1D
implicit none

type(integrand_euler_1D) :: euler_1D !< Euler 1D integrand field.

call euler_1D%integrate
endprogram foodie_euler_1D
