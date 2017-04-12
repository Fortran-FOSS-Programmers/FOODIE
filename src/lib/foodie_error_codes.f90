!< FOODIE error codes list.

module foodie_error_codes
!< FOODIE error codes list.

use penf, only : I_P

implicit none
private
public :: ERROR_UNSUPPORTED_SCHEME

integer(I_P), parameter :: ERROR_UNSUPPORTED_SCHEME = 1 !< Error unsupported scheme.
endmodule foodie_error_codes
