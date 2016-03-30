!< PENF's testing program.
program test_all
!-----------------------------------------------------------------------------------------------------------------------------------
!< PENF's testing program.
!<
!<### Usage
!<```bash
!< ./test_all
!<```
!-----------------------------------------------------------------------------------------------------------------------------------
use penf
use, intrinsic :: ISO_FORTRAN_ENV, only : stdout=>OUTPUT_UNIT
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
call penf_Init
call penf_Print(unit=stdout)
print "(A)", ''
print "(A)", 'Testing IR_Precision tools'
print "(A)", 'Casting real-to-string: '//str(n=1._R8P)
print "(A)", 'Casting integer-to-string: '//str(n=11_I8P)
print "(A,"//FR8P//")", 'Casting string-to-real: ',cton(str='2.2d0', knd=1._R8P)
print "(A,"//FI4P//")", 'Casting integer-to-string: ',cton(str='43', knd=1_I4P)
print "(A)", 'Casting integer-to-string with zero padding: '//trim(strz(nz_pad=3, n=34_I8P))
#ifndef __GFORTRAN__
print "(A)", 'Casting real-to-bit_string: '//bstr(n=1._R4P)
#endif
print "(A)", 'Casting integer-to-bit_string: '//bstr(n=1_I4P)
#ifndef __GFORTRAN__
print "(A,"//FR4P//")", 'Casting bit_string-to-real: ', bcton(bstr='00111111100000000000000000000000', knd=1._R4P)
#endif
print "(A,"//FI4P//")", 'Casting bit_string-to-integer: ',bcton(bstr='00000000000000000000000000000001', knd=1_I4P)
print "(A)", 'Number of digit of 1023: '//str(n=digit(1023_I4P))
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram test_all
