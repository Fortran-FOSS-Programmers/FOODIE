!< FOODIE error codes list.

module foodie_error_codes
!< FOODIE error codes list.

use foodie_kinds, only : I_P

implicit none
private
public :: ERROR_MISSING_STEPS_NUMBER
public :: ERROR_BAD_STEPS_NUMBER

integer(I_P), parameter :: ERROR_MISSING_STEPS_NUMBER=1 !< Error code: missing steps number (for FOODIE factory).
integer(I_P), parameter :: ERROR_BAD_STEPS_NUMBER=2     !< Error code: bad (unsupported) steps number.
endmodule foodie_error_codes
