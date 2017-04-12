!< FOODIE utils: module of (possible) unrelated utilities of FOODIE library.

module foodie_utils
!< FOODIE utils: module of (possible) unrelated utilities of FOODIE library.

use penf, only : I_P

implicit none
private
public :: is_admissible

contains
  ! public procedures
  elemental function is_admissible(n, adm_range)
  !< Check if the queried number *n* is admitted by the *admissible* range list *adm_range*.
  !<
  !< The admissible range list must be formatted as string containing admissible numbers; valid list are:
  !<+ `adm_range = '1-5'` => 1, 2, 3, 4, 5 are admissible numbers;
  !<+ `adm_range = '1,3,5,10-12'` => 1, 3, 5, 10, 11, 12 are admissible numbers;
  !<+ `adm_range = '1-4,8,21-22'` => 1, 2, 3, 4, 8, 21, 22 are admissible numbers;
  !<
  !< You can mix any number of range (`min-max` format) and/or single number (`,` comma separated) entries.
  integer(I_P), intent(IN)               :: n             !< Number queried.
  character(*), intent(IN)               :: adm_range     !< Admissible range string.
  logical                                :: is_admissible !< Is true is the number is in *adm_range*.
  character(len(adm_range)), allocatable :: tokens(:)     !< Tokens for parsing *adm_range* string.
  character(len(adm_range)), allocatable :: subtokens(:)  !< Tokens for parsing *adm_range* string.
  integer(I_P)                           :: t             !< Counter.
  integer(I_P)                           :: n_parsed(1:2) !< Values parsed from *adm_range*..

  is_admissible = .false.
  call tokenize(string=adm_range, delimiter=',', toks=tokens)
  search_me : do t=1, size(tokens)
    if (index(tokens(t), '-')>0) then
      call tokenize(string=tokens(t), delimiter='-', toks=subtokens)
      read(subtokens(1), *) n_parsed(1)
      read(subtokens(2), *) n_parsed(2)
      is_admissible = (n_parsed(1)<=n.and.n<=n_parsed(2))
    else
      read(tokens(t), *) n_parsed(1)
      is_admissible = (n_parsed(1)==n)
    endif
    if (is_admissible) exit search_me
  enddo search_me
  endfunction is_admissible

  ! private procedures
  pure subroutine tokenize(string, delimiter, toks, Nt)
  !< Tokenize a string in order to parse it.
  !<
  !< @note The dummy array containing tokens must be allocatable and its character elements must have the same length of the input
  !< string. If the length of the delimiter is higher than the input string one then the output tokens array is allocated with
  !< only one element set to char(0).
  character(len=*),                        intent(IN)  :: string    !< String to be tokenized.
  character(len=*),                        intent(IN)  :: delimiter !< Delimiter of tokens.
  character(len=len(string)), allocatable, intent(OUT) :: toks(:)   !< Tokens.
  integer(I_P), optional,                  intent(OUT) :: Nt        !< Number of tokens.
  character(len=len(string))                           :: strsub    !< Temporary string.
  integer(I_P)                                         :: dlen      !< Delimiter length.
  integer(I_P)                                         :: c         !< Counter.
  integer(I_P)                                         :: n         !< Counter.
  integer(I_P)                                         :: t         !< Counter.

  ! initialize
  if (allocated(toks)) deallocate(toks)
  strsub = string
  dlen = len(delimiter)
  if (dlen>len(string)) then
    allocate(toks(1:1))
    toks(1) = char(0)
    if (present(Nt)) Nt = 1
    return
  endif

  ! compute the number of tokens
  n = 1
  do c=1, len(strsub) - dlen ! loop over string characters
    if (strsub(c:c+dlen-1)==delimiter) n = n + 1
  enddo
  allocate(toks(1:n))

  ! tokenize
  do t=1, n ! loop over tokens
    c = index(strsub, delimiter)
    if (c>0) then
      toks(t) = strsub(1:c-1)
      strsub = strsub(c+dlen:)
    else
      toks(t) = strsub
    endif
  enddo
  if (present(Nt)) Nt = n
  endsubroutine tokenize
endmodule foodie_utils
