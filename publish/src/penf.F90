!< Portability Environment for Fortran poor people.

module penf_global_parameters_variables
!-----------------------------------------------------------------------------------------------------------------------------------
!< PENF global (exposed) parameters and variables.
!<
!< @note All module defined entities are public.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
public
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
logical            :: is_initialized = .false. !< Check the initialization of some variables that must be initialized.
integer, parameter :: endianL        = 1       !< Little endian parameter.
integer, parameter :: endianB        = 0       !< Big endian parameter.
integer            :: endian         = endianL !< Bit ordering: Little endian (endianL), or Big endian (endianB).
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Portable kind parameters
#ifdef r16p
integer, parameter :: R16P = selected_real_kind(33,4931) !< 33 digits, range \([10^{-4931}, 10^{+4931} - 1]\); 128 bits.
#else
integer, parameter :: R16P = selected_real_kind(15,307)  !< Defined as R8P; 64 bits.
#endif
integer, parameter :: R8P  = selected_real_kind(15,307)  !< 15 digits, range \([10^{-307} , 10^{+307}  - 1]\); 64 bits.
integer, parameter :: R4P  = selected_real_kind(6,37)    !< 6  digits, range \([10^{-37}  , 10^{+37}   - 1]\); 32 bits.
integer, parameter :: R_P  = R8P                         !< Default real precision.
integer, parameter :: I8P  = selected_int_kind(18)       !< Range \([-2^{63},+2^{63} - 1]\), 19 digits plus sign; 64 bits.
integer, parameter :: I4P  = selected_int_kind(9)        !< Range \([-2^{31},+2^{31} - 1]\), 10 digits plus sign; 32 bits.
integer, parameter :: I2P  = selected_int_kind(4)        !< Range \([-2^{15},+2^{15} - 1]\), 5  digits plus sign; 16 bits.
integer, parameter :: I1P  = selected_int_kind(2)        !< Range \([-2^{7} ,+2^{7}  - 1]\), 3  digits plus sign; 8  bits.
integer, parameter :: I_P  = I4P                         !< Default integer precision.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Format parameters
#ifdef r16p
character(*), parameter :: FR16P  = '(E42.33E4)' !< Output format for kind=R16P real.
#else
character(*), parameter :: FR16P  = '(E23.15E3)' !< Output format for kind=R16P real.
#endif
character(*), parameter :: FR8P   = '(E23.15E3)' !< Output format for kind=R8P real.
character(*), parameter :: FR4P   = '(E13.6E2)'  !< Output format for kind=R4P real.
character(*), parameter :: FR_P   = FR8P         !< Output format for kind=R_P real.
character(*), parameter :: FI8P   = '(I20)'      !< Output format for kind=I8P integer.
character(*), parameter :: FI8PZP = '(I20.19)'   !< Output format for kind=I8P integer with zero prefixing.
character(*), parameter :: FI4P   = '(I11)'      !< Output format for kind=I4P integer.
character(*), parameter :: FI4PZP = '(I11.10)'   !< Output format for kind=I4P integer with zero prefixing.
character(*), parameter :: FI2P   = '(I6)'       !< Output format for kind=I2P integer.
character(*), parameter :: FI2PZP = '(I6.5)'     !< Output format for kind=I2P integer with zero prefixing.
character(*), parameter :: FI1P   = '(I4)'       !< Output format for kind=I1P integer.
character(*), parameter :: FI1PZP = '(I4.3)'     !< Output format for kind=I1P integer with zero prefixing.
character(*), parameter :: FI_P   = FI4P         !< Output format for kind=I_P integer.
character(*), parameter :: FI_PZP = FI4PZP       !< Output format for kind=I_P integer with zero prefixing.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Length (number of digits) of formatted numbers
#ifdef r16p
integer, parameter :: DR16P = 42   !< Number of digits of output format FR16P.
#else
integer, parameter :: DR16P = 23   !< Number of digits of output format FR16P.
#endif
integer, parameter :: DR8P  = 23   !< Number of digits of output format FR8P.
integer, parameter :: DR4P  = 13   !< Number of digits of output format FR4P.
integer, parameter :: DR_P  = DR8P !< Number of digits of output format FR_P.
integer, parameter :: DI8P  = 20   !< Number of digits of output format I8P.
integer, parameter :: DI4P  = 11   !< Number of digits of output format I4P.
integer, parameter :: DI2P  = 6    !< Number of digits of output format I2P.
integer, parameter :: DI1P  = 4    !< Number of digits of output format I1P.
integer, parameter :: DI_P  = DI4P !< Number of digits of output format I_P.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! List of kinds
integer,      parameter :: REAL_KINDS_LIST(1:4)      = [R16P, R8P, R4P, R_P]                    !< List of real kinds.
character(*), parameter :: REAL_FORMATS_LIST(1:4)    = [FR16P, FR8P, FR4P//' ', FR_P]           !< List of real formats.
integer,      parameter :: INTEGER_KINDS_LIST(1:5)   = [I8P, I4P, I2P, I1P,I_P]                 !< List of integer kinds.
character(*), parameter :: INTEGER_FORMATS_LIST(1:5) = [FI8P, FI4P, FI2P//' ', FI1P//' ', FI_P] !< List of integer formats.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Minimum and maximum (representable) values
real(R16P),   parameter :: MinR16P = -huge(1._R16P) !< Minimum value of kind=R16P real.
real(R16P),   parameter :: MaxR16P =  huge(1._R16P) !< Maximum value of kind=R16P real.
real(R8P),    parameter :: MinR8P  = -huge(1._R8P ) !< Minimum value of kind=R8P real.
real(R8P),    parameter :: MaxR8P  =  huge(1._R8P ) !< Maximum value of kind=R8P real.
real(R4P),    parameter :: MinR4P  = -huge(1._R4P ) !< Minimum value of kind=R4P real.
real(R4P),    parameter :: MaxR4P  =  huge(1._R4P ) !< Maximum value of kind=R4P real.
real(R_P),    parameter :: MinR_P  = MinR8P         !< Minimum value of kind=R_P real.
real(R_P),    parameter :: MaxR_P  = MaxR8P         !< Maximum value of kind=R_P real.
integer(I8P), parameter :: MinI8P  = -huge(1_I8P)   !< Minimum value of kind=I8P integer.
integer(I4P), parameter :: MinI4P  = -huge(1_I4P)   !< Minimum value of kind=I4P integer.
integer(I2P), parameter :: MinI2P  = -huge(1_I2P)   !< Minimum value of kind=I2P integer.
integer(I1P), parameter :: MinI1P  = -huge(1_I1P)   !< Minimum value of kind=I1P integer.
integer(I_P), parameter :: MinI_P  = MinI4P         !< Minimum value of kind=I_P integer.
integer(I8P), parameter :: MaxI8P  =  huge(1_I8P)   !< Maximum value of kind=I8P integer.
integer(I4P), parameter :: MaxI4P  =  huge(1_I4P)   !< Maximum value of kind=I4P integer.
integer(I2P), parameter :: MaxI2P  =  huge(1_I2P)   !< Maximum value of kind=I2P integer.
integer(I1P), parameter :: MaxI1P  =  huge(1_I1P)   !< Maximum value of kind=I1P integer.
integer(I_P), parameter :: MaxI_P  =  MaxI4P        !< Maximum value of kind=I_P integer.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Real smallest (representable) values
real(R16P), parameter :: smallR16P = tiny(1._R16P) !< Smallest representable value of kind=R16P real.
real(R8P),  parameter :: smallR8P  = tiny(1._R8P ) !< Smallest representable value of kind=R8P real.
real(R4P),  parameter :: smallR4P  = tiny(1._R4P ) !< Smallest representable value of kind=R4P real.
real(R_P),  parameter :: smallR_P  = smallR8P      !< Smallest representable value of kind=R_P real.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Smallest real representable difference by the running calculator
real(R16P), parameter :: ZeroR16 = nearest(1._R16P, 1._R16P) - &
                                   nearest(1._R16P,-1._R16P) !< Smallest representable difference of kind=R16P real.
real(R8P),  parameter :: ZeroR8  = nearest(1._R8P, 1._R8P) - &
                                   nearest(1._R8P,-1._R8P)   !< Smallest representable difference of kind=R8P real.
real(R4P),  parameter :: ZeroR4  = nearest(1._R4P, 1._R4P) - &
                                   nearest(1._R4P,-1._R4P)   !< Smallest representable difference of kind=R4P real.
real(R_P),  parameter :: Zero    = ZeroR8                    !< Smallest representable difference of kind=R_P real.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Bits/bytes memory requirements (real variables must be computed at runtime)
integer(I2P)            :: BIR16P                         !< Number of bits of kind=R16P real.
integer(I1P)            :: BIR8P                          !< Number of bits of kind=R8P real.
integer(I1P)            :: BIR4P                          !< Number of bits of kind=R4P real.
integer(I1P)            :: BIR_P                          !< Number of bits of kind=R_P real.
integer(I2P)            :: BYR16P                         !< Number of bytes of kind=R16P real.
integer(I1P)            :: BYR8P                          !< Number of bytes of kind=R8P real.
integer(I1P)            :: BYR4P                          !< Number of bytes of kind=R4P real.
integer(I1P)            :: BYR_P                          !< Number of bytes of kind=R_P real.
integer(I8P), parameter :: BII8P = bit_size(MaxI8P)       !< Number of bits of kind=I8P integer.
integer(I4P), parameter :: BII4P = bit_size(MaxI4P)       !< Number of bits of kind=I4P integer.
integer(I2P), parameter :: BII2P = bit_size(MaxI2P)       !< Number of bits of kind=I2P integer.
integer(I1P), parameter :: BII1P = bit_size(MaxI1P)       !< Number of bits of kind=I1P integer.
integer(I_P), parameter :: BII_P = bit_size(MaxI_P)       !< Number of bits of kind=I_P integer.
integer(I8P), parameter :: BYI8P = bit_size(MaxI8P)/8_I8P !< Number of bytes of kind=I8P integer.
integer(I4P), parameter :: BYI4P = bit_size(MaxI4P)/8_I4P !< Number of bytes of kind=I4P integer.
integer(I2P), parameter :: BYI2P = bit_size(MaxI2P)/8_I2P !< Number of bytes of kind=I2P integer.
integer(I1P), parameter :: BYI1P = bit_size(MaxI1P)/8_I1P !< Number of bytes of kind=I1P integer.
integer(I_P), parameter :: BYI_P = bit_size(MaxI_P)/8_I_P !< Number of bytes of kind=I_P integer.
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule penf_global_parameters_variables

module penf_b_size
!-----------------------------------------------------------------------------------------------------------------------------------
!< PENF bit/byte size functions.
!-----------------------------------------------------------------------------------------------------------------------------------
use penf_global_parameters_variables
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: bit_size, byte_size
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
interface bit_size
  !< Overloading of the intrinsic *bit_size* function for computing the number of bits of (also) real and character variables.
  module procedure                &
#ifdef r16p
                   bit_size_R16P, &
#endif
                   bit_size_R8P,  &
                   bit_size_R4P,  &
                   bit_size_chr
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
interface byte_size
  !< Compute the number of bytes of a variable.
  module procedure                 &
                   byte_size_I8P,  &
                   byte_size_I4P,  &
                   byte_size_I2P,  &
                   byte_size_I1P,  &
#ifdef r16p
                   byte_size_R16p, &
#endif
                   byte_size_R8P,  &
                   byte_size_R4P,  &
                   byte_size_chr
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function bit_size_R16P(i) result(bits)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bits of a real variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R16P), intent(in) :: i       !< Real variable whose number of bits must be computed.
  integer(I2P)           :: bits    !< Number of bits of r.
  integer(I1P)           :: mold(1) !< "Molding" dummy variable for bits counting.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bits = size(transfer(i, mold), dim=1, kind=I2P) * 8_I2P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bit_size_R16P

  elemental function bit_size_R8P(i) result(bits)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bits of a real variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R8P), intent(in) :: i       !< Real variable whose number of bits must be computed.
  integer(I1P)          :: bits    !< Number of bits of r.
  integer(I1P)          :: mold(1) !< "Molding" dummy variable for bits counting.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bits = size(transfer(i, mold), dim=1, kind=I1P) * 8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bit_size_R8P

  elemental function bit_size_R4P(i) result(bits)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bits of a real variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R4P), intent(in) :: i       !< Real variable whose number of bits must be computed.
  integer(I1P)          :: bits    !< Number of bits of r.
  integer(I1P)          :: mold(1) !< "Molding" dummy variable for bits counting.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bits = size(transfer(i, mold), dim=1, kind=I1P) * 8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bit_size_R4P

  elemental function bit_size_chr(i) result(bits)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bits of a character variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN) :: i       !< Character variable whose number of bits must be computed.
  integer(I4P)             :: bits    !< Number of bits of c.
  integer(I1P)             :: mold(1) !< "Molding" dummy variable for bits counting.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bits = size(transfer(i, mold), dim=1, kind=I4P) * 8_I4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bit_size_chr

  elemental function byte_size_I8P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bytes of an integer variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(in) :: i     !< Integer variable whose number of bytes must be computed.
  integer(I1P)             :: bytes !< Number of bytes of i.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_I8P

  elemental function byte_size_I4P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bytes of an integer variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in) :: i     !< Integer variable whose number of bytes must be computed.
  integer(I1P)             :: bytes !< Number of bytes of i.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_I4P

  elemental function byte_size_I2P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bytes of an integer variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I2P), intent(in) :: i     !< Integer variable whose number of bytes must be computed.
  integer(I1P)             :: bytes !< Number of bytes of i.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_I2P

  elemental function byte_size_I1P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bytes of an integer variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I1P), intent(in) :: i     !< Integer variable whose number of bytes must be computed.
  integer(I1P)             :: bytes !< Number of bytes of i.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_I1P

  elemental function byte_size_R16P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bytes of a real variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R16P), intent(in) :: i     !< Real variable whose number of bytes must be computed.
  integer(I1P)           :: bytes !< Number of bytes of r.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_R16P

  elemental function byte_size_R8P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bytes of a real variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R8P), intent(in) :: i     !< Real variable whose number of bytes must be computed.
  integer(I1P)          :: bytes !< Number of bytes of r.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_R8P

  elemental function byte_size_R4P(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bytes of a real variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R4P), intent(in) :: i     !< Real variable whose number of bytes must be computed.
  integer(I1P)          :: bytes !< Number of bytes of r.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I1P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_R4P

  elemental function byte_size_chr(i) result(bytes)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of bytes of a character variable.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: i     !< Character variable whose number of bytes must be computed.
  integer(I4P)             :: bytes !< Number of bytes of c.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bytes = bit_size(i)/8_I4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction byte_size_chr
endmodule penf_b_size

module penf_stringify
!-----------------------------------------------------------------------------------------------------------------------------------
!< PENF string-to-number (and viceversa) facility.
!-----------------------------------------------------------------------------------------------------------------------------------
use, intrinsic :: ISO_FORTRAN_ENV, only : stderr => ERROR_UNIT
use penf_b_size
use penf_global_parameters_variables
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: str, strz, cton
public :: bstr, bcton
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
interface str
  !< Convert number (real and integer) to string (number to string type casting).
  module procedure                       &
#ifdef r16p
                   strf_R16P,str_R16P,   &
#endif
                   strf_R8P ,str_R8P,    &
                   strf_R4P ,str_R4P,    &
                   strf_I8P ,str_I8P,    &
                   strf_I4P ,str_I4P,    &
                   strf_I2P ,str_I2P,    &
                   strf_I1P ,str_I1P,    &
                             str_bol,    &
#ifdef r16p
                             str_a_R16P, &
#endif
                             str_a_R8P,  &
                             str_a_R4P,  &
                             str_a_I8P,  &
                             str_a_I4P,  &
                             str_a_I2P,  &
                             str_a_I1P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
interface strz
  !< Convert integer, to string, prefixing with the right number of zeros (integer to string type casting with zero padding).
  module procedure strz_I8P, strz_I4P, strz_I2P, strz_I1P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
interface cton
  !< Convert string to number (real and integer, string to number type casting).
  module procedure            &
#ifdef r16p
                   ctor_R16P, &
#endif
                   ctor_R8P,  &
                   ctor_R4P,  &
                   ctoi_I8P,  &
                   ctoi_I4P,  &
                   ctoi_I2P,  &
                   ctoi_I1P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
interface bstr
  !< Convert number (real and integer) to bit-string (number to bit-string type casting).
  module procedure            &
#ifdef r16p
                   bstr_R16P, &
#endif
                   bstr_R8P,  &
                   bstr_R4P,  &
                   bstr_I8P,  &
                   bstr_I4P,  &
                   bstr_I2P,  &
                   bstr_I1P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
interface bcton
  !< Convert bit-string to number (real and integer, bit-string to number type casting).
  module procedure             &
#ifdef r16p
                   bctor_R16P, &
#endif
                   bctor_R8P,  &
                   bctor_R4P,  &
                   bctoi_I8P,  &
                   bctoi_I4P,  &
                   bctoi_I2P,  &
                   bctoi_I1P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function strf_R16P(fm, n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: fm  !< Format different from the standard for the kind.
  real(R16P),   intent(in) :: n   !< Real to be converted.
  character(DR16P)         :: str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, trim(fm)) n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_R16P

  elemental function strf_R8P(fm, n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: fm  !< Format different from the standard for the kind.
  real(R8P),    intent(in) :: n   !< Real to be converted.
  character(DR8P)          :: str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, trim(fm)) n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_R8P

  elemental function strf_R4P(fm, n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: fm  !< Format different from the standard for the kind.
  real(R4P),    intent(in) :: n   !< Real to be converted.
  character(DR4P)          :: str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, trim(fm)) n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_R4P

  elemental function strf_I8P(fm, n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: fm  !< Format different from the standard for the kind.
  integer(I8P), intent(in) :: n   !< Integer to be converted.
  character(DI8P)          :: str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, trim(fm)) n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_I8P

  elemental function strf_I4P(fm, n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: fm  !< Format different from the standard for the kind.
  integer(I4P), intent(in) :: n   !< Integer to be converted.
  character(DI4P)          :: str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, trim(fm)) n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_I4P

  elemental function strf_I2P(fm, n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: fm  !< Format different from the standard for the kind.
  integer(I2P), intent(in) :: n   !< Integer to be converted.
  character(DI2P)          :: str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, trim(fm)) n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_I2P

  elemental function strf_I1P(fm, n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: fm  !< Format different from the standard for the kind.
  integer(I1P), intent(in) :: n   !< Integer to be converted.
  character(DI1P)          :: str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, trim(fm)) n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strf_I1P

  elemental function str_R16P(n, no_sign) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R16P), intent(in)           :: n       !< Real to be converted.
  logical,    intent(in), optional :: no_sign !< Flag for leaving out the sign.
  character(DR16P)                 :: str     !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FR16P) n               ! Casting of n to string.
  if (n>0._R16P) str(1:1)='+'       ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_R16P

  elemental function str_R8P(n, no_sign) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R8P), intent(in)           :: n       !< Real to be converted.
  logical,   intent(in), optional :: no_sign !< Flag for leaving out the sign.
  character(DR8P)                 :: str     !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FR8P) n                ! Casting of n to string.
  if (n>0._R8P) str(1:1)='+'        ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_R8P

  elemental function str_R4P(n, no_sign) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R4P), intent(in)           :: n       !< Real to be converted.
  logical,   intent(in), optional :: no_sign !< Flag for leaving out the sign.
  character(DR4P)                 :: str     !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FR4P) n                ! Casting of n to string.
  if (n>0._R4P) str(1:1)='+'        ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_R4P

  elemental function str_I8P(n, no_sign) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(in)           :: n       !< Integer to be converted.
  logical,      intent(in), optional :: no_sign !< Flag for leaving out the sign.
  character(DI8P)                    :: str     !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FI8P) n                ! Casting of n to string.
  str = adjustl(trim(str))          ! Removing white spaces.
  if (n>=0_I8P) str='+'//trim(str)  ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_I8P

  elemental function str_I4P(n, no_sign) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Converting integer to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in)           :: n       !< Integer to be converted.
  logical,      intent(in), optional :: no_sign !< Flag for leaving out the sign.
  character(DI4P)                    :: str     !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FI4P) n                ! Casting of n to string.
  str = adjustl(trim(str))          ! Removing white spaces.
  if (n>=0_I4P) str='+'//trim(str)  ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_I4P

  elemental function str_I2P(n, no_sign) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I2P), intent(in)           :: n       !< Integer to be converted.
  logical,      intent(in), optional :: no_sign !< Flag for leaving out the sign.
  character(DI2P)                    :: str     !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FI2P) n                ! Casting of n to string.
  str = adjustl(trim(str))          ! Removing white spaces.
  if (n>=0_I2P) str='+'//trim(str)  ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_I2P

  elemental function str_I1P(n, no_sign) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I1P), intent(in)           :: n       !< Integer to be converted.
  logical,      intent(in), optional :: no_sign !< Flag for leaving out the sign.
  character(DI1P)                    :: str     !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FI1P) n                ! Casting of n to string.
  str = adjustl(trim(str))          ! Removing white spaces.
  if (n>=0_I1P) str='+'//trim(str)  ! Prefixing plus if n>0.
  if (present(no_sign)) str=str(2:) ! Leaving out the sign.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_I1P

  elemental function str_bol(n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert logical to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  logical, intent(in):: n   !< Logical to be converted.
  character(1)::        str !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, '(L1)') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_bol

  pure function str_a_R16P(n, no_sign, delimiters) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Converting real array to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R16P),   intent(in)           :: n(:)            !< Real array to be converted.
  logical,      intent(in), optional :: no_sign         !< Flag for leaving out the sign.
  character(*), intent(in), optional :: delimiters(1:2) !< Eventual delimiters of array values.
  character(len=:), allocatable      :: str             !< Returned string containing input number.
  character(DR16P)                   :: strn            !< String containing of element of input array number.
  integer                            :: i               !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(no_sign)) then
    str = ''
    do i=1,size(n)
      strn = str_R16P(no_sign=no_sign, n=n(i))
      str = str//','//trim(strn)
    enddo
  else
    str = ''
    do i=1,size(n)
      strn = str_R16P(n=n(i))
      str = str//','//trim(strn)
    enddo
  endif
  str = trim(str(2:))
  if (present(delimiters)) str = delimiters(1)//str//delimiters(2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_a_R16P

  pure function str_a_R8P(n, no_sign, delimiters) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real array to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R8P),    intent(in)           :: n(:)            !< Real array to be converted.
  logical,      intent(in), optional :: no_sign         !< Flag for leaving out the sign.
  character(*), intent(in), optional :: delimiters(1:2) !< Eventual delimiters of array values.
  character(len=:), allocatable      :: str             !< Returned string containing input number.
  character(DR8P)                    :: strn            !< String containing of element of input array number.
  integer                            :: i               !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(no_sign)) then
    str = ''
    do i=1,size(n)
      strn = str_R8P(no_sign=no_sign, n=n(i))
      str = str//','//trim(strn)
    enddo
  else
    str = ''
    do i=1,size(n)
      strn = str_R8P(n=n(i))
      str = str//','//trim(strn)
    enddo
  endif
  str = trim(str(2:))
  if (present(delimiters)) str = delimiters(1)//str//delimiters(2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_a_R8P

  pure function str_a_R4P(n, no_sign, delimiters) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real array to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R4P),    intent(in)           :: n(:)            !< Real array to be converted.
  logical,      intent(in), optional :: no_sign         !< Flag for leaving out the sign.
  character(*), intent(in), optional :: delimiters(1:2) !< Eventual delimiters of array values.
  character(len=:), allocatable      :: str             !< Returned string containing input number.
  character(DR4P)                    :: strn            !< String containing of element of input array number.
  integer                            :: i               !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(no_sign)) then
    str = ''
    do i=1,size(n)
      strn = str_R4P(no_sign=no_sign, n=n(i))
      str = str//','//trim(strn)
    enddo
  else
    str = ''
    do i=1,size(n)
      strn = str_R4P(n=n(i))
      str = str//','//trim(strn)
    enddo
  endif
  str = trim(str(2:))
  if (present(delimiters)) str = delimiters(1)//str//delimiters(2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_a_R4P

  pure function str_a_I8P(n, no_sign, delimiters) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer array to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(in)           :: n(:)            !< Integer array to be converted.
  logical,      intent(in), optional :: no_sign         !< Flag for leaving out the sign.
  character(*), intent(in), optional :: delimiters(1:2) !< Eventual delimiters of array values.
  character(len=:), allocatable      :: str             !< Returned string containing input number.
  character(DI8P)                    :: strn            !< String containing of element of input array number.
  integer                            :: i               !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(no_sign)) then
    str = ''
    do i=1,size(n)
      strn = str_I8P(no_sign=no_sign, n=n(i))
      str = str//','//trim(strn)
    enddo
  else
    str = ''
    do i=1,size(n)
      strn = str_I8P(n=n(i))
      str = str//','//trim(strn)
    enddo
  endif
  str = trim(str(2:))
  if (present(delimiters)) str = delimiters(1)//str//delimiters(2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_a_I8P

  pure function str_a_I4P(n, no_sign, delimiters) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer array to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in)           :: n(:)            !< Integer array to be converted.
  logical,      intent(in), optional :: no_sign         !< Flag for leaving out the sign.
  character(*), intent(in), optional :: delimiters(1:2) !< Eventual delimiters of array values.
  character(len=:), allocatable      :: str             !< Returned string containing input number.
  character(DI4P)                    :: strn            !< String containing of element of input array number.
  integer                            :: i               !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(no_sign)) then
    str = ''
    do i=1,size(n)
      strn = str_I4P(no_sign=no_sign, n=n(i))
      str = str//','//trim(strn)
    enddo
  else
    str = ''
    do i=1,size(n)
      strn = str_I4P(n=n(i))
      str = str//','//trim(strn)
    enddo
  endif
  str = trim(str(2:))
  if (present(delimiters)) str = delimiters(1)//str//delimiters(2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_a_I4P

  pure function str_a_I2P(n, no_sign, delimiters) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer array to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I2P), intent(in)           :: n(:)            !< Integer array to be converted.
  logical,      intent(in), optional :: no_sign         !< Flag for leaving out the sign.
  character(*), intent(in), optional :: delimiters(1:2) !< Eventual delimiters of array values.
  character(len=:), allocatable      :: str             !< Returned string containing input number.
  character(DI2P)                    :: strn            !< String containing of element of input array number.
  integer                            :: i               !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(no_sign)) then
    str = ''
    do i=1,size(n)
      strn = str_I2P(no_sign=no_sign, n=n(i))
      str = str//','//trim(strn)
    enddo
  else
    str = ''
    do i=1,size(n)
      strn = str_I2P(n=n(i))
      str = str//','//trim(strn)
    enddo
  endif
  str = trim(str(2:))
  if (present(delimiters)) str = delimiters(1)//str//delimiters(2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_a_I2P

  pure function str_a_I1P(n, no_sign, delimiters) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer array to string.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I1P), intent(in)           :: n(:)            !< Integer array to be converted.
  logical,      intent(in), optional :: no_sign         !< Flag for leaving out the sign.
  character(*), intent(in), optional :: delimiters(1:2) !< Eventual delimiters of array values.
  character(len=:), allocatable      :: str             !< Returned string containing input number.
  character(DI1P)                    :: strn            !< String containing of element of input array number.
  integer                            :: i               !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(no_sign)) then
    str = ''
    do i=1,size(n)
      strn = str_I1P(no_sign=no_sign, n=n(i))
      str = str//','//trim(strn)
    enddo
  else
    str = ''
    do i=1,size(n)
      strn = str_I1P(n=n(i))
      str = str//','//trim(strn)
    enddo
  endif
  str = trim(str(2:))
  if (present(delimiters)) str = delimiters(1)//str//delimiters(2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_a_I1P

  elemental function strz_I8P(n, nz_pad) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Converting integer to string, prefixing with the right number of zeros.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(in)           :: n      !< Integer to be converted.
  integer(I4P), intent(in), optional :: nz_pad !< Number of zeros padding.
  character(DI8P)                    :: str    !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI8PZP) n                              ! Casting of n to string.
  str=str(2:)                                      ! Leaving out the sign.
  if (present(nz_pad)) str=str(DI8P-nz_pad:DI8P-1) ! Leaving out the extra zeros padding
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strz_I8P

  elemental function strz_I4P(n, nz_pad) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string, prefixing with the right number of zeros.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in)           :: n      !< Integer to be converted.
  integer(I4P), intent(in), optional :: nz_pad !< Number of zeros padding.
  character(DI4P)                    :: str    !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI4PZP) n                              ! Casting of n to string.
  str=str(2:)                                      ! Leaving out the sign.
  if (present(nz_pad)) str=str(DI4P-nz_pad:DI4P-1) ! Leaving out the extra zeros padding
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strz_I4P

  elemental function strz_I2P(n, nz_pad) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string, prefixing with the right number of zeros.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I2P), intent(in)           :: n      !< Integer to be converted.
  integer(I4P), intent(in), optional :: nz_pad !< Number of zeros padding.
  character(DI2P)                    :: str    !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI2PZP) n                              ! Casting of n to string.
  str=str(2:)                                      ! Leaving out the sign.
  if (present(nz_pad)) str=str(DI2P-nz_pad:DI2P-1) ! Leaving out the extra zeros padding
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strz_I2P

  elemental function strz_I1P(n, nz_pad) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string, prefixing with the right number of zeros.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I1P), intent(in)           :: n      !< Integer to be converted.
  integer(I4P), intent(in), optional :: nz_pad !< Number of zeros padding.
  character(DI1P)                    :: str    !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI1PZP) n                              ! Casting of n to string.
  str=str(2:)                                      ! Leaving out the sign.
  if (present(nz_pad)) str=str(DI1P-nz_pad:DI1P-1) ! Leaving out the extra zeros padding
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strz_I1P

  function ctor_R16P(str, knd, pref, error) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert string to real.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*),           intent(in)  :: str   !< String containing input number.
  real(R16P),             intent(in)  :: knd   !< Number kind.
  character(*), optional, intent(in)  :: pref  !< Prefixing string.
  integer(I4P), optional, intent(out) :: error !< Error trapping flag: 0 no errors, >0 error occurs.
  real(R16P)                          :: n     !< Number returned.
  integer(I4P)                        :: err   !< Error trapping flag: 0 no errors, >0 error occurs.
  character(len=:), allocatable       :: prefd !< Prefixing string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str, *, iostat=err) n ! Casting of str to n.
  if (err/=0) then
    prefd = '' ; if (present(pref)) prefd = pref
    write(stderr, '(A,I1,A)') prefd//' Error: conversion of string "'//str//'" to real failed! real(', kind(knd), ')'
  endif
  if (present(error)) error = err
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctor_R16P

  function ctor_R8P(str, knd, pref, error) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert string to real.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*),           intent(in)  :: str   !< String containing input number.
  real(R8P),              intent(in)  :: knd   !< Number kind.
  character(*), optional, intent(in)  :: pref  !< Prefixing string.
  integer(I4P), optional, intent(out) :: error !< Error trapping flag: 0 no errors, >0 error occurs.
  real(R8P)                           :: n     !< Number returned.
  integer(I4P)                        :: err   !< Error trapping flag: 0 no errors, >0 error occurs.
  character(len=:), allocatable       :: prefd !< Prefixing string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str, *, iostat=err) n ! Casting of str to n.
  if (err/=0) then
    prefd = '' ; if (present(pref)) prefd = pref
    write(stderr, '(A,I1,A)') prefd//' Error: conversion of string "'//str//'" to real failed! real(', kind(knd), ')'
  endif
  if (present(error)) error = err
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctor_R8P

  function ctor_R4P(str, knd, pref, error) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert string to real.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*),           intent(in)  :: str   !< String containing input number.
  real(R4P),              intent(in)  :: knd   !< Number kind.
  character(*), optional, intent(in)  :: pref  !< Prefixing string.
  integer(I4P), optional, intent(out) :: error !< Error trapping flag: 0 no errors, >0 error occurs.
  real(R4P)                           :: n     !< Number returned.
  integer(I4P)                        :: err   !< Error trapping flag: 0 no errors, >0 error occurs.
  character(len=:), allocatable       :: prefd !< Prefixing string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str, *, iostat=err) n ! Casting of str to n.
  if (err/=0) then
    prefd = '' ; if (present(pref)) prefd = pref
    write(stderr, '(A,I1,A)') prefd//' Error: conversion of string "'//str//'" to real failed! real(', kind(knd), ')'
  endif
  if (present(error)) error = err
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctor_R4P

  function ctoi_I8P(str, knd, pref, error) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert string to integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*),           intent(in)  :: str   !< String containing input number.
  integer(I8P),           intent(in)  :: knd   !< Number kind.
  character(*), optional, intent(in)  :: pref  !< Prefixing string.
  integer(I4P), optional, intent(out) :: error !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I8P)                        :: n     !< Number returned.
  integer(I4P)                        :: err   !< Error trapping flag: 0 no errors, >0 error occurs.
  character(len=:), allocatable       :: prefd !< Prefixing string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str, *, iostat=err) n ! Casting of str to n.
  if (err/=0) then
    prefd = '' ; if (present(pref)) prefd = pref
    write(stderr, '(A,I1,A)') prefd//' Error: conversion of string "'//str//'" to integer failed! integer(', kind(knd), ')'
  endif
  if (present(error)) error = err
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctoi_I8P

  function ctoi_I4P(str, knd, pref, error) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert string to integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*),           intent(in)  :: str   !< String containing input number.
  integer(I4P),           intent(in)  :: knd   !< Number kind.
  character(*), optional, intent(in)  :: pref  !< Prefixing string.
  integer(I4P), optional, intent(out) :: error !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I4P)                        :: n     !< Number returned.
  integer(I4P)                        :: err   !< Error trapping flag: 0 no errors, >0 error occurs.
  character(len=:), allocatable       :: prefd !< Prefixing string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str, *, iostat=err) n ! Casting of str to n.
  if (err/=0) then
    prefd = '' ; if (present(pref)) prefd = pref
    write(stderr, '(A,I1,A)') prefd//' Error: conversion of string "'//str//'" to integer failed! integer(', kind(knd), ')'
  endif
  if (present(error)) error = err
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctoi_I4P

  function ctoi_I2P(str, knd, pref, error) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert string to integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*),           intent(in)  :: str   !< String containing input number.
  integer(I2P),           intent(in)  :: knd   !< Number kind.
  character(*), optional, intent(in)  :: pref  !< Prefixing string.
  integer(I4P), optional, intent(out) :: error !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I2P)                        :: n     !< Number returned.
  integer(I4P)                        :: err   !< Error trapping flag: 0 no errors, >0 error occurs.
  character(len=:), allocatable       :: prefd !< Prefixing string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str, *, iostat=err) n ! Casting of str to n.
  if (err/=0) then
    prefd = '' ; if (present(pref)) prefd = pref
    write(stderr, '(A,I1,A)') prefd//' Error: conversion of string "'//str//'" to integer failed! integer(', kind(knd), ')'
  endif
  if (present(error)) error = err
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctoi_I2P

  function ctoi_I1P(str, knd, pref, error) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert string to integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*),           intent(in)  :: str   !< String containing input number.
  integer(I1P),           intent(in)  :: knd   !< Number kind.
  character(*), optional, intent(in)  :: pref  !< Prefixing string.
  integer(I4P), optional, intent(out) :: error !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I1P)                        :: n     !< Number returned.
  integer(I4P)                        :: err   !< Error trapping flag: 0 no errors, >0 error occurs.
  character(len=:), allocatable       :: prefd !< Prefixing string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(str, *, iostat=err) n ! Casting of str to n.
  if (err/=0) then
    prefd = '' ; if (present(pref)) prefd = pref
    write(stderr, '(A,I1,A)') prefd//' Error: conversion of string "'//str//'" to integer failed! integer(', kind(knd), ')'
  endif
  if (present(error)) error = err
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ctoi_I1P

  elemental function bstr_R16P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real to string of bits.
  !<
  !< @note It is assumed that R16P is represented by means of 128 bits, but this is not ensured in all architectures.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R16P), intent(in) :: n    !< Real to be converted.
  character(128)         :: bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr, '(B128.128)') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_R16P

  elemental function bstr_R8P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real to string of bits.
  !<
  !< @note It is assumed that R8P is represented by means of 64 bits, but this is not ensured in all architectures.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R8P), intent(in) :: n    !< Real to be converted.
  character(64)         :: bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr, '(B64.64)') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_R8P

  elemental function bstr_R4P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert real to string of bits.
  !<
  !< @note It is assumed that R4P is represented by means of 32 bits, but this is not ensured in all architectures.
  !---------------------------------------------------------------------------------------------------------------------------------
  real(R4P), intent(in) :: n    !< Real to be converted.
  character(32)         :: bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr, '(B32.32)') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_R4P

  elemental function bstr_I8P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string of bits.
  !<
  !< @note It is assumed that I8P is represented by means of 64 bits, but this is not ensured in all architectures.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(in) :: n    !< Real to be converted.
  character(64)            :: bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr, '(B64.64)') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_I8P

  elemental function bstr_I4P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string of bits.
  !<
  !< @note It is assumed that I4P is represented by means of 32 bits, but this is not ensured in all architectures.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in) :: n    !< Real to be converted.
  character(32)            :: bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr, '(B32.32)') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_I4P

  elemental function bstr_I2P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string of bits.
  !<
  !< @note It is assumed that I2P is represented by means of 16 bits, but this is not ensured in all architectures.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I2P), intent(in) :: n    !< Real to be converted.
  character(16)            :: bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr, '(B16.16)') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_I2P

  elemental function bstr_I1P(n) result(bstr)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert integer to string of bits.
  !<
  !< @note It is assumed that I1P is represented by means of 8 bits, but this is not ensured in all architectures.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I1P), intent(in) :: n    !< Real to be converted.
  character(8)             :: bstr !< Returned bit-string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(bstr, '(B8.8)') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bstr_I1P

  elemental function bctor_R16P(bstr, knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert bit-string to real.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: bstr !< String containing input number.
  real(R16P),   intent(in) :: knd  !< Number kind.
  real(R16P)               :: n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr, '(B'//trim(str(bit_size(knd), .true.))//'.'//trim(str(bit_size(knd), .true.))//')') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctor_R16P

  elemental function bctor_R8P(bstr, knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert bit-string to real.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: bstr !< String containing input number.
  real(R8P),    intent(in) :: knd  !< Number kind.
  real(R8P)                :: n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr, '(B'//trim(str(bit_size(knd), .true.))//'.'//trim(str(bit_size(knd), .true.))//')') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctor_R8P

  elemental function bctor_R4P(bstr, knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert bit-string to real.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: bstr !< String containing input number.
  real(R4P),    intent(in) :: knd  !< Number kind.
  real(R4P)                :: n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(bit_size(knd), .true.))//'.'//trim(str(bit_size(knd), .true.))//')') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctor_R4P

  elemental function bctoi_I8P(bstr, knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert bit-string to integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: bstr !< String containing input number.
  integer(I8P), intent(in) :: knd  !< Number kind.
  integer(I8P)             :: n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(bit_size(knd), .true.))//'.'//trim(str(bit_size(knd), .true.))//')') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctoi_I8P

  elemental function bctoi_I4P(bstr, knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert bit-string to integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: bstr !< String containing input number.
  integer(I4P), intent(in) :: knd  !< Number kind.
  integer(I4P)             :: n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(bit_size(knd), .true.))//'.'//trim(str(bit_size(knd), .true.))//')') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctoi_I4P

  elemental function bctoi_I2P(bstr, knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert bit-string to integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: bstr !< String containing input number.
  integer(I2P), intent(in) :: knd  !< Number kind.
  integer(I2P)             :: n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(bit_size(knd), .true.))//'.'//trim(str(bit_size(knd), .true.))//')') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctoi_I2P

  elemental function bctoi_I1P(bstr, knd) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Convert bit-string to integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(in) :: bstr !< String containing input number.
  integer(I1P), intent(in) :: knd  !< Number kind.
  integer(I1P)             :: n    !< Number returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(bstr,'(B'//trim(str(bit_size(knd), .true.))//'.'//trim(str(bit_size(knd), .true.))//')') n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bctoi_I1P
endmodule penf_stringify

module penf
!-----------------------------------------------------------------------------------------------------------------------------------
!< Portability Environment for Fortran poor people.
!<{!README-PENF.md!}
!-----------------------------------------------------------------------------------------------------------------------------------
use penf_global_parameters_variables
#ifdef __GFORTRAN__
use penf_b_size, only : bit_size, byte_size
#else
use penf_b_size
#endif
use penf_stringify, only : str, strz, cton, bstr, bcton
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
! Global parameters and variables
public :: endianL, endianB, endian, is_initialized
public :: R16P, FR16P, DR16P, MinR16P, MaxR16P, BIR16P, BYR16P, smallR16P, ZeroR16
public :: R8P,  FR8P,  DR8P,  MinR8P,  MaxR8P,  BIR8P,  BYR8P,  smallR8P,  ZeroR8
public :: R4P,  FR4P,  DR4P,  MinR4P,  MaxR4P,  BIR4P,  BYR4P,  smallR4P,  ZeroR4
public :: R_P,  FR_P,  DR_P,  MinR_P,  MaxR_P,  BIR_P,  BYR_P,  smallR_P,  Zero
public :: I8P,  FI8P,  DI8P,  MinI8P,  MaxI8P,  BII8P,  BYI8P
public :: I4P,  FI4P,  DI4P,  MinI4P,  MaxI4P,  BII4P,  BYI4P
public :: I2P,  FI2P,  DI2P,  MinI2P,  MaxI2P,  BII2P,  BYI2P
public :: I1P,  FI1P,  DI1P,  MinI1P,  MaxI1P,  BII1P,  BYI1P
public :: I_P,  FI_P,  DI_P,  MinI_P,  MaxI_P,  BII_P,  BYI_P
public :: REAL_KINDS_LIST, REAL_FORMATS_LIST
public :: INTEGER_KINDS_LIST, INTEGER_FORMATS_LIST
! Bit/byte size functions
public :: bit_size, byte_size
! Stringify facility
public :: str, strz, cton
public :: bstr, bcton
! Miscellanea facility
public :: check_endian
public :: digit
public :: penf_Init
public :: penf_print
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef __GFORTRAN__
! work-around for strange gfortran bug...
interface bit_size
  !< Overloading of the intrinsic *bit_size* function for computing the number of bits of (also) real and character variables.
endinterface
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

interface digit
  !< Compute the number of digits in decimal base of the input integer.
  module procedure digit_I8, digit_I4, digit_I2, digit_I1
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine check_endian()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Check the type of bit ordering (big or little endian) of the running architecture.
  !<
  !> @note The result is stored into the *endian* global variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (is_little_endian()) then
    endian = endianL
  else
    endian = endianB
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure function is_little_endian() result(is_little)
    !-------------------------------------------------------------------------------------------------------------------------------
    !< Check if the type of the bit ordering of the running architecture is little endian.
    !-------------------------------------------------------------------------------------------------------------------------------
    logical      :: is_little !< Logical output: true is the running architecture uses little endian ordering, false otherwise.
    integer(I1P) :: int1(1:4) !< One byte integer array for casting 4 bytes integer.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    int1 = transfer(1_I4P, int1)
    is_little = (int1(1)==1_I1P)
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction is_little_endian
  endsubroutine check_endian

  elemental function digit_I8(n) result(digit)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of digits in decimal base of the input integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(in) :: n     !< Input integer.
  character(DI8P)          :: str   !< Returned string containing input number plus padding zeros.
  integer(I4P)             :: digit !< Number of digits.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FI8P) abs(n)        ! Casting of n to string.
  digit = len_trim(adjustl(str)) ! Calculating the digits number of n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction digit_I8

  elemental function digit_I4(n) result(digit)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of digits in decimal base of the input integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in) :: n     !< Input integer.
  character(DI4P)          :: str   !< Returned string containing input number plus padding zeros.
  integer(I4P)             :: digit !< Number of digits.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FI4P) abs(n)        ! Casting of n to string.
  digit = len_trim(adjustl(str)) ! Calculating the digits number of n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction digit_I4

  elemental function digit_I2(n) result(digit)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of digits in decimal base of the input integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I2P), intent(in) :: n     !< Input integer.
  character(DI2P)          :: str   !< Returned string containing input number plus padding zeros.
  integer(I4P)             :: digit !< Number of digits.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FI2P) abs(n)        ! Casting of n to string.
  digit = len_trim(adjustl(str)) ! Calculating the digits number of n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction digit_I2

  elemental function digit_I1(n) result(digit)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the number of digits in decimal base of the input integer.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I1P), intent(in) :: n     !< Input integer.
  character(DI1P)          :: str   !< Returned string containing input number plus padding zeros.
  integer(I4P)             :: digit !< Number of digits.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str, FI1P) abs(n)        ! Casting of n to string.
  digit = len_trim(adjustl(str)) ! Calculating the digits number of n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction digit_I1

  subroutine penf_init()
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Initialize PENF's variables that are not initialized into the definition specification.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call check_endian
  BIR8P  = bit_size(MaxR8P)     ; BYR8P  = BIR8P/8_I1P
  BIR4P  = bit_size(MaxR4P)     ; BYR4P  = BIR4P/8_I1P
  BIR_P  = bit_size(MaxR_P)     ; BYR_P  = BIR_P/8_I1P
#ifdef r16p
  BIR16P = bit_size(MaxR16P)    ; BYR16P = BIR16P/8_I2P
#else
  BIR16P = int(BIR8P, kind=I2P) ; BYR16P = BIR16P/8_I2P
#endif
  is_initialized = .true.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine penf_init

  subroutine penf_print(unit, pref, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Print to the specified unit the PENF's environment data.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in)            :: unit    !< Logic unit.
  character(*), intent(in),  optional :: pref    !< Prefixing string.
  integer(I4P), intent(out), optional :: iostat  !< IO error.
  character(*), intent(out), optional :: iomsg   !< IO error message.
  character(len=:), allocatable       :: prefd   !< Prefixing string.
  integer(I4P)                        :: iostatd !< IO error.
  character(500)                      :: iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.is_initialized) call penf_init
  prefd = '' ; if (present(pref)) prefd = pref
  if (endian==endianL) then
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' This architecture has LITTLE Endian bit ordering'
  else
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' This architecture has BIG Endian bit ordering'
  endif
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//' Reals kind, format and characters number:'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   R16P: '//str(n=R16P)//','//FR16P//','//str(n=DR16P)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   R8P:  '//str(n=R8P )//','//FR8P //','//str(n=DR8P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   R4P:  '//str(n=R4P )//','//FR4P //','//str(n=DR4P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//' Integers kind, format and characters number:'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I8P:  '//str(n=I8P )//','//FI8P //','//str(n=DI8P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I4P:  '//str(n=I4P )//','//FI4P //','//str(n=DI4P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I2P:  '//str(n=I2P )//','//FI2P //','//str(n=DI2P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I1P:  '//str(n=I1P )//','//FI1P //','//str(n=DI1P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//' Reals minimum and maximum values:'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   R16P: '//str(n=MinR16P)//','//str(n=MaxR16P)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   R8P:  '//str(n=MinR8P )//','//str(n=MaxR8P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   R4P:  '//str(n=MinR4P )//','//str(n=MaxR4P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//' Integergs minimum and maximum values:'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I8P:  '//str(n=MinI8P )//','//str(n=MaxI8P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I4P:  '//str(n=MinI4P )//','//str(n=MaxI4P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I2P:  '//str(n=MinI2P )//','//str(n=MaxI2P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I1P:  '//str(n=MinI1P )//','//str(n=MaxI1P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//' Reals bits/bytes sizes:'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   R16P: '//str(n=BIR16P)//'/'//str(n=BYR16P)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   R8P:  '//str(n=BIR8P )//'/'//str(n=BYR8P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   R4P:  '//str(n=BIR4P )//'/'//str(n=BYR4P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//' Integers bits/bytes sizes:'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I8P:  '//str(n=BII8P )//'/'//str(n=BYI8P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I4P:  '//str(n=BII4P )//'/'//str(n=BYI4P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I2P:  '//str(n=BII2P )//'/'//str(n=BYI2P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   I1P:  '//str(n=BII1P )//'/'//str(n=BYI1P )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//' Machine precisions'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   ZeroR16: '//str(ZeroR16,.true.)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   ZeroR8:  '//str(ZeroR8 ,.true.)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)  prefd//'   ZeroR4:  '//str(ZeroR4 ,.true.)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine penf_print
endmodule penf
