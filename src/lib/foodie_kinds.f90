!< FOODIE kinds: definition of reals and integer kind parameters of FOODIE library.

module foodie_kinds
!< FOODIE kinds: definition of reals and integer kind parameters of FOODIE library.

implicit none
private
public :: R8P
public :: R4P
public :: R_P
public :: I8P
public :: I4P
public :: I2P
public :: I1P
public :: I_P

integer, parameter :: R8P  = selected_real_kind(15,307) !< 15  digits, range \([10^{-307} , 10^{+307}  - 1]\); 64 bits.
integer, parameter :: R4P  = selected_real_kind(6,37)   !< 6   digits, range \([10^{-37}  , 10^{+37}   - 1]\); 32 bits.
integer, parameter :: R_P  = R8P                        !< Default real precision.
integer, parameter :: I8P  = selected_int_kind(18)      !< Range \([-2^{63},+2^{63} - 1]\), 19 digits plus sign; 64 bits.
integer, parameter :: I4P  = selected_int_kind(9)       !< Range \([-2^{31},+2^{31} - 1]\), 10 digits plus sign; 32 bits.
integer, parameter :: I2P  = selected_int_kind(4)       !< Range \([-2^{15},+2^{15} - 1]\), 5  digits plus sign; 16 bits.
integer, parameter :: I1P  = selected_int_kind(2)       !< Range \([-2^{7} ,+2^{7}  - 1]\), 3  digits plus sign; 8  bits.
integer, parameter :: I_P  = I4P                        !< Default integer precision.
endmodule foodie_kinds
