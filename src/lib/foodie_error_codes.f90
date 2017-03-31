!< FOODIE error codes list.

module foodie_error_codes
!< FOODIE error codes list.

use foodie_kinds, only : I_P

implicit none
private
public :: ERROR_MISSING_STEPS_NUMBER
public :: ERROR_BAD_STEPS_NUMBER
public :: ERROR_MISSING_STAGES_NUMBER
public :: ERROR_BAD_STAGES_NUMBER
public :: ERROR_INTEGRATOR_INIT_FAIL

integer(I_P), parameter :: ERROR_MISSING_STEPS_NUMBER  = 1 !< Error missing steps number (for FOODIE factory).
integer(I_P), parameter :: ERROR_BAD_STEPS_NUMBER      = 2 !< Error bad (unsupported) steps number.
integer(I_P), parameter :: ERROR_MISSING_STAGES_NUMBER = 3 !< Error missing stages number (for FOODIE factory).
integer(I_P), parameter :: ERROR_BAD_STAGES_NUMBER     = 4 !< Error bad (unsupported) stages number.
integer(I_P), parameter :: ERROR_INTEGRATOR_INIT_FAIL  = 5 !< Error integrator initialization failed.
endmodule foodie_error_codes
