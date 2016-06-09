!< StringiFor, definition of `string` type.
module stringifor_string_t
!-----------------------------------------------------------------------------------------------------------------------------------
!< StringiFor, definition of `string` type.
!-----------------------------------------------------------------------------------------------------------------------------------
use befor64, only : b64_decode, b64_encode
use penf, only : I1P, I2P, I4P, I8P, R4P, R8P, R16P, str
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: CK
public :: sadjustl_character, sadjustr_character,                                 &
          sindex_string_string, sindex_string_character, sindex_character_string, &
          slen, slen_trim,                                                        &
          srepeat_string_string,                                                  &
          sscan_string_string, sscan_string_character, sscan_character_string,    &
          strim
public :: string
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
integer, parameter :: CK = selected_char_kind('DEFAULT') !< Default character kind.

type :: string
  !< OOP designed string class.
  private
  character(kind=CK, len=:), allocatable :: raw !< Raw data.
  contains
    ! public methods
    ! builtins replacements
    procedure, pass(self) :: adjustl  => sadjustl                 !< Adjustl replacement.
    procedure, pass(self) :: adjustr  => sadjustr                 !< Adjustr replacement.
    procedure, pass(self) :: count    => scount                   !< Count replacement.
    generic               :: index    => sindex_string_string, &
                                         sindex_string_character  !< Index replacement.
    procedure, pass(self) :: len      => slen                     !< Len replacement.
    procedure, pass(self) :: len_trim => slen_trim                !< Len_trim replacement.
    generic               :: repeat   => srepeat_string_string, &
                                         srepeat_character_string !< Repeat replacement.
    generic               :: scan     => sscan_string_string,    &
                                         sscan_string_character   !< Scan replacement.
    procedure, pass(self) :: trim     => strim                    !< Trim replacement.
    procedure, pass(self) :: verify   => sverify                  !< Verify replacement.
    ! auxiliary methods
    procedure, pass(self) :: basedir          !< Return the base directory name of a string containing a file name.
    procedure, pass(self) :: basename         !< Return the base file name of a string containing a file name.
    procedure, pass(self) :: camelcase        !< Return a string with all words capitalized without spaces.
    procedure, pass(self) :: capitalize       !< Return a string with its first character capitalized and the rest lowercased.
    procedure, pass(self) :: chars            !< Return the raw characters data.
    procedure, pass(self) :: decode           !< Decode string.
    procedure, pass(self) :: encode           !< Encode string.
    procedure, pass(self) :: escape           !< Escape backslashes (or custom escape character).
    procedure, pass(self) :: extension        !< Return the extension of a string containing a file name.
    procedure, pass(self) :: fill             !< Pad string on the left (or right) with zeros (or other char) to fill width.
    procedure, pass(self) :: free             !< Free dynamic memory.
    generic               :: insert =>      &
                             insert_string, &
                             insert_character !< Insert substring into string at a specified position.
    generic               :: join =>       &
                             join_strings, &
                             join_characters  !< Return a string that is a join of an array of strings or characters.
    procedure, pass(self) :: lower            !< Return a string with all lowercase characters.
    procedure, pass(self) :: partition        !< Split string at separator and return the 3 parts (before, the separator and after).
    procedure, pass(self) :: read_file        !< Read a file a single string stream.
    procedure, pass(self) :: read_line        !< Read line (record) from a connected unit.
    procedure, pass(self) :: read_lines       !< Read (all) lines (records) from a connected unit as a single ascii stream.
    procedure, pass(self) :: replace          !< Return a string with all occurrences of substring old replaced by new.
    procedure, pass(self) :: reverse          !< Return a reversed string.
    procedure, pass(self) :: search           !< Search for *tagged* record into string.
    procedure, pass(self) :: slice            !< Return the raw characters data sliced.
    procedure, pass(self) :: snakecase        !< Return a string with all words lowercase separated by "_".
    procedure, pass(self) :: split            !< Return a list of substring in the string, using sep as the delimiter string.
    procedure, pass(self) :: startcase        !< Return a string with all words capitalized, e.g. title case.
    procedure, pass(self) :: strip            !< Return a string with the leading and trailing characters removed.
    procedure, pass(self) :: swapcase         !< Return a string with uppercase chars converted to lowercase and vice versa.
    generic               :: to_number =>   &
                             to_integer_I1P,&
                             to_integer_I2P,&
                             to_integer_I4P,&
                             to_integer_I8P,&
                             to_real_R4P,   &
#ifdef r16p
                             to_real_R8P,   &
                             to_real_R16P     !< Cast string to number.
#else
                             to_real_R8P      !< Cast string to number.
#endif
    procedure, pass(self) :: unescape         !< Unescape double backslashes (or custom escaped character).
    procedure, pass(self) :: unique           !< Reduce to one (unique) multiple occurrences of a substring into a string.
    procedure, pass(self) :: upper            !< Return a string with all uppercase characters.
    procedure, pass(self) :: write_file       !< Write a single string stream into file.
    procedure, pass(self) :: write_line       !< Write line (record) to a connected unit.
    procedure, pass(self) :: write_lines      !< Write lines (records) to a connected unit.
    ! inquire methods
    procedure, pass(self) :: end_with     !< Return true if a string ends with a specified suffix.
    procedure, pass(self) :: is_allocated !< Return true if the string is allocated.
    procedure, pass(self) :: is_digit     !< Return true if all characters in the string are digits.
    procedure, pass(self) :: is_integer   !< Return true if the string contains an integer.
    procedure, pass(self) :: is_lower     !< Return true if all characters in the string are lowercase.
    procedure, pass(self) :: is_number    !< Return true if the string contains a number (real or integer).
    procedure, pass(self) :: is_real      !< Return true if the string contains an real.
    procedure, pass(self) :: is_upper     !< Return true if all characters in the string are uppercase.
    procedure, pass(self) :: start_with   !< Return true if a string starts with a specified prefix.
    ! operators
    generic :: assignment(=) => string_assign_string,      &
                                string_assign_character,   &
                                string_assign_integer_I1P, &
                                string_assign_integer_I2P, &
                                string_assign_integer_I4P, &
                                string_assign_integer_I8P, &
                                string_assign_real_R4P,    &
#ifdef r16p
                                string_assign_real_R8P,    &
                                string_assign_real_R16P             !< Assignment operator overloading.
#else
                                string_assign_real_R8P              !< Assignment operator overloading.
#endif
    generic :: operator(//) => string_concat_string,    &
                               string_concat_character, &
                               character_concat_string              !< Concatenation operator overloading.
    generic :: operator(.cat.) => string_concat_string_string,    &
                                  string_concat_character_string, &
                                  character_concat_string_string    !< Concatenation operator (string output) overloading.
    generic :: operator(==) => string_eq_string,    &
                               string_eq_character, &
                               character_eq_string                  !< Equal operator overloading.
    generic :: operator(/=) => string_ne_string,    &
                               string_ne_character, &
                               character_ne_string                  !< Not equal operator overloading.
    generic :: operator(<) => string_lt_string,    &
                              string_lt_character, &
                              character_lt_string                   !< Lower than operator overloading.
    generic :: operator(<=) => string_le_string,    &
                               string_le_character, &
                               character_le_string                  !< Lower equal than operator overloading.
    generic :: operator(>=) => string_ge_string,    &
                               string_ge_character, &
                               character_ge_string                  !< Greater equal than operator overloading.
    generic :: operator(>) => string_gt_string,    &
                              string_gt_character, &
                              character_gt_string                   !< Greater than operator overloading.
    ! IO
#ifndef __GFORTRAN__
    generic :: read(formatted) => read_formatted       !< Formatted input.
    generic :: write(formatted) => write_formatted     !< Formatted output.
    generic :: read(unformatted) => read_unformatted   !< Unformatted input.
    generic :: write(unformatted) => write_unformatted !< Unformatted output.
#endif
    ! private methods
    ! builtins replacements
    procedure, private, pass(self) :: sindex_string_string     !< Index replacement.
    procedure, private, pass(self) :: sindex_string_character  !< Index replacement.
    procedure, private, pass(self) :: srepeat_string_string    !< Repeat replacement.
    procedure, private, pass(self) :: srepeat_character_string !< Repeat replacement.
    procedure, private, pass(self) :: sscan_string_string      !< Scan replacement.
    procedure, private, pass(self) :: sscan_string_character   !< Scan replacement.
    ! auxiliary methods
    procedure, private, pass(self) :: insert_string    !< Insert substring into string at a specified position.
    procedure, private, pass(self) :: insert_character !< Insert substring into string at a specified position.
    procedure, private, pass(self) :: join_strings     !< Return join string of an array of strings.
    procedure, private, pass(self) :: join_characters  !< Return join string of an array of characters.
    procedure, private, pass(self) :: to_integer_I1P   !< Cast string to integer.
    procedure, private, pass(self) :: to_integer_I2P   !< Cast string to integer.
    procedure, private, pass(self) :: to_integer_I4P   !< Cast string to integer.
    procedure, private, pass(self) :: to_integer_I8P   !< Cast string to integer.
    procedure, private, pass(self) :: to_real_R4P      !< Cast string to real.
    procedure, private, pass(self) :: to_real_R8P      !< Cast string to real.
    procedure, private, pass(self) :: to_real_R16P     !< Cast string to real.
    ! assignments
    procedure, private, pass(lhs) :: string_assign_string      !< Assignment operator from string input.
    procedure, private, pass(lhs) :: string_assign_character   !< Assignment operator from character input.
    procedure, private, pass(lhs) :: string_assign_integer_I1P !< Assignment operator from integer input.
    procedure, private, pass(lhs) :: string_assign_integer_I2P !< Assignment operator from integer input.
    procedure, private, pass(lhs) :: string_assign_integer_I4P !< Assignment operator from integer input.
    procedure, private, pass(lhs) :: string_assign_integer_I8P !< Assignment operator from integer input.
    procedure, private, pass(lhs) :: string_assign_real_R4P    !< Assignment operator from real input.
    procedure, private, pass(lhs) :: string_assign_real_R8P    !< Assignment operator from real input.
    procedure, private, pass(lhs) :: string_assign_real_R16P   !< Assignment operator from real input.
    ! concatenation operators
    procedure, private, pass(lhs) :: string_concat_string           !< Concatenation with string.
    procedure, private, pass(lhs) :: string_concat_character        !< Concatenation with character.
    procedure, private, pass(rhs) :: character_concat_string        !< Concatenation with character (inverted).
    procedure, private, pass(lhs) :: string_concat_string_string    !< Concatenation with string (string output).
    procedure, private, pass(lhs) :: string_concat_character_string !< Concatenation with character (string output).
    procedure, private, pass(rhs) :: character_concat_string_string !< Concatenation with character (inverted, string output).
    ! logical operators
    procedure, private, pass(lhs) :: string_eq_string    !< Equal to string logical operator.
    procedure, private, pass(lhs) :: string_eq_character !< Equal to character logical operator.
    procedure, private, pass(rhs) :: character_eq_string !< Equal to character (inverted) logical operator.
    procedure, private, pass(lhs) :: string_ne_string    !< Not equal to string logical operator.
    procedure, private, pass(lhs) :: string_ne_character !< Not equal to character logical operator.
    procedure, private, pass(rhs) :: character_ne_string !< Not equal to character (inverted) logical operator.
    procedure, private, pass(lhs) :: string_lt_string    !< Lower than to string logical operator.
    procedure, private, pass(lhs) :: string_lt_character !< Lower than to character logical operator.
    procedure, private, pass(rhs) :: character_lt_string !< Lower than to character (inverted) logical operator.
    procedure, private, pass(lhs) :: string_le_string    !< Lower equal than to string logical operator.
    procedure, private, pass(lhs) :: string_le_character !< Lower equal than to character logical operator.
    procedure, private, pass(rhs) :: character_le_string !< Lower equal than to character (inverted) logical operator.
    procedure, private, pass(lhs) :: string_ge_string    !< Greater equal than to string logical operator.
    procedure, private, pass(lhs) :: string_ge_character !< Greater equal than to character logical operator.
    procedure, private, pass(rhs) :: character_ge_string !< Greater equal than to character (inverted) logical operator.
    procedure, private, pass(lhs) :: string_gt_string    !< Greater than to string logical operator.
    procedure, private, pass(lhs) :: string_gt_character !< Greater than to character logical operator.
    procedure, private, pass(rhs) :: character_gt_string !< Greater than to character (inverted) logical operator.
    ! IO
    procedure, private, pass(dtv) :: read_formatted                !< Formatted input.
    procedure, private, pass(dtv) :: read_delimited                !< Read a delimited input.
    procedure, private, pass(dtv) :: read_undelimited              !< Read an undelimited input.
    procedure, private, pass(dtv) :: read_undelimited_listdirected !< Read an undelimited list directed input.
    procedure, private, pass(dtv) :: write_formatted               !< Formatted output.
    procedure, private, pass(dtv) :: read_unformatted              !< Unformatted input.
    procedure, private, pass(dtv) :: write_unformatted             !< Unformatted output.
    ! miscellanea
    procedure, private, pass(self) :: replace_one_occurrence !< Replace the first occurrence of substring old by new.
endtype string

! internal parameters
character(kind=CK, len=26), parameter :: UPPER_ALPHABET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' !< Upper case alphabet.
character(kind=CK, len=26), parameter :: LOWER_ALPHABET = 'abcdefghijklmnopqrstuvwxyz' !< Lower case alphabet.
character(kind=CK, len=1),  parameter :: SPACE          = ' '                          !< Space character.
character(kind=CK, len=1),  parameter :: TAB            = achar(9)                     !< Tab character.
character(kind=CK, len=1),  parameter :: UIX_DIR_SEP    = char(47)                     !< Unix/Linux directories separator (/).
character(kind=CK, len=1),  parameter :: BACKSLASH      = char(92)                     !< Backslash character.
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods

  ! builtins replacements
  elemental function sadjustl(self) result(adjusted)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Left adjust a string by removing leading spaces.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self     !< The string.
  type(string)              :: adjusted !< Adjusted string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  adjusted = self
  if (allocated(adjusted%raw)) adjusted%raw = adjustl(adjusted%raw)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sadjustl

  pure function sadjustl_character(self) result(adjusted)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Left adjust a string by removing leading spaces (character output).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)             :: self     !< The string.
  character(kind=CK, len=len(self%raw)) :: adjusted !< Adjusted string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) adjusted = adjustl(self%raw)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sadjustl_character

  elemental function sadjustr(self) result(adjusted)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Right adjust a string by removing leading spaces.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self     !< The string.
  type(string)              :: adjusted !< Adjusted string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  adjusted = self
  if (allocated(adjusted%raw)) adjusted%raw = adjustr(adjusted%raw)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sadjustr

  pure function sadjustr_character(self) result(adjusted)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Right adjust a string by removing leading spaces (character output).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)             :: self     !< The string.
  character(kind=CK, len=len(self%raw)) :: adjusted !< Adjusted string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) adjusted = adjustr(self%raw)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sadjustr_character

  elemental function scount(self, substring, ignore_isolated) result(No)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Count the number of occurences of a substring into a string.
  !<
  !< @note If `ignore_isolated` is set to true the eventual "isolated" occurences are ignored: an isolated occurrences are those
  !< occurrences happening at the start of string (thus not having a left companion) or at the end of the string (thus not having a
  !< right companion).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)              :: self             !< The string.
  character(*),  intent(in)              :: substring        !< Substring.
  logical,       intent(in), optional    :: ignore_isolated  !< Ignore "isolated" occurrences.
  integer                                :: No               !< Number of occurrences.
  logical                                :: ignore_isolated_ !< Ignore "isolated" occurrences, local variable.
  integer                                :: c1               !< Counter.
  integer                                :: c2               !< Counter.
#ifdef __GFORTRAN__
  character(kind=CK, len=:), allocatable :: temporary        !< Temporary storage, workaround for GNU bug.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  No = 0
  if (allocated(self%raw)) then
    if (len(substring)>len(self%raw)) return
    ignore_isolated_ = .false. ; if (present(ignore_isolated)) ignore_isolated_ = ignore_isolated
#ifdef __GFORTRAN__
    temporary = self%raw
#endif
    c1 = 1
    do
#ifdef __GFORTRAN__
      c2 = index(string=temporary(c1:), substring=substring)
#else
      c2 = index(string=self%raw(c1:), substring=substring)
#endif
      if (c2==0) return
      if (.not.(ignore_isolated_.and.(c1==1.or.c1+c2-1==len(self%raw)-len(substring)+1))) then
        No = No + 1
      endif
      c1 = c1 + c2 - 1 + len(substring)
    enddo
  endif
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction scount

  elemental function sindex_string_string(self, substring, back) result(i)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the position of the start of the first occurrence of string `substring` as a substring in `string`, counting from one.
  !< If `substring` is not present in `string`, zero is returned. If the back argument is present and true, the return value is
  !< the start of the last occurrence rather than the first.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)           :: self      !< The string.
  type(string),  intent(in)           :: substring !< Searched substring.
  logical,       intent(in), optional :: back      !< Start of the last occurrence rather than the first.
  integer                             :: i         !< Result of the search.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    i = index(string=self%raw, substring=substring%raw, back=back)
  else
    i = 0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sindex_string_string

  elemental function sindex_string_character(self, substring, back) result(i)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the position of the start of the first occurrence of string `substring` as a substring in `string`, counting from one.
  !< If `substring` is not present in `string`, zero is returned. If the back argument is present and true, the return value is
  !< the start of the last occurrence rather than the first.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self      !< The string.
  character(kind=CK, len=*), intent(in)           :: substring !< Searched substring.
  logical,                   intent(in), optional :: back      !< Start of the last occurrence rather than the first.
  integer                                         :: i         !< Result of the search.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    i = index(string=self%raw, substring=substring, back=back)
  else
    i = 0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sindex_string_character

  elemental function sindex_character_string(string_, substring, back) result(i)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the position of the start of the first occurrence of string `substring` as a substring in `string`, counting from one.
  !< If `substring` is not present in `string`, zero is returned. If the back argument is present and true, the return value is
  !< the start of the last occurrence rather than the first.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in)           :: string_   !< The string.
  type(string),              intent(in)           :: substring !< Searched substring.
  logical,                   intent(in), optional :: back      !< Start of the last occurrence rather than the first.
  integer                                         :: i         !< Result of the search.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(substring%raw)) then
    i = index(string=string_, substring=substring%raw, back=back)
  else
    i = 0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sindex_character_string

  elemental function slen(self) result(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the length of a string.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self !< The string.
  integer                   :: l    !< String length.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    l = len(string=self%raw)
  else
    l = 0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction slen

  elemental function slen_trim(self) result(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the length of a string, ignoring any trailing blanks.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self !< The string.
  integer                   :: l    !< String length.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    l = len_trim(string=self%raw)
  else
    l = 0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction slen_trim

  elemental function srepeat_string_string(self, ncopies) result(repeated)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Concatenates several copies of an input string.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self     !< String to be repeated.
  integer,       intent(in) :: ncopies  !< Number of string copies.
  type(string)              :: repeated !< Repeated string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  repeated%raw = repeat(string=self%raw, ncopies=ncopies)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction srepeat_string_string

  elemental function srepeat_character_string(self, rstring, ncopies) result(repeated)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Concatenates several copies of an input string.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in) :: self     !< String to be repeated.
  character(kind=CK, len=*), intent(in) :: rstring  !< String to be repeated.
  integer,                   intent(in) :: ncopies  !< Number of string copies.
  type(string)                          :: repeated !< Repeated string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  repeated%raw = repeat(string=rstring, ncopies=ncopies)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction srepeat_character_string

  elemental function sscan_string_string(self, set, back) result(i)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the leftmost (if `back` is either absent or equals false, otherwise the rightmost) character of string that is in `set`.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)           :: self  !< The string.
  type(string),  intent(in)           :: set   !< Searched set.
  logical,       intent(in), optional :: back  !< Start of the last occurrence rather than the first.
  integer                             :: i     !< Result of the search.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw).and.allocated(set%raw)) then
    i = scan(string=self%raw, set=set%raw, back=back)
  else
    i = 0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sscan_string_string

  elemental function sscan_string_character(self, set, back) result(i)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the leftmost (if `back` is either absent or equals false, otherwise the rightmost) character of string that is in `set`.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self  !< The string.
  character(kind=CK, len=*), intent(in)           :: set   !< Searched set.
  logical,                   intent(in), optional :: back  !< Start of the last occurrence rather than the first.
  integer                                         :: i     !< Result of the search.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    i = scan(string=self%raw, set=set, back=back)
  else
    i = 0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sscan_string_character

  elemental function sscan_character_string(sstring, set, back) result(i)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the leftmost (if `back` is either absent or equals false, otherwise the rightmost) character of string that is in `set`.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in)           :: sstring !< The string.
  type(string),              intent(in)           :: set     !< Searched set.
  logical,                   intent(in), optional :: back    !< Start of the last occurrence rather than the first.
  integer                                         :: i       !< Result of the search.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(set%raw)) then
    i = scan(string=sstring, set=set%raw, back=back)
  else
    i = 0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sscan_character_string

  elemental function strim(self) result(trimmed)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Remove leading spaces.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self    !< The string.
  type(string)              :: trimmed !< Trimmed string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  trimmed = self
  if (allocated(trimmed%raw)) trimmed%raw = trim(trimmed%raw)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strim

  elemental function sverify(self, set, back) result(i)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the leftmost (if `back` is either absent or equals false, otherwise the rightmost) character of string that is not
  !< in `set`. If all characters of `string` are found in `set`, the result is zero.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self  !< The string.
  character(kind=CK, len=*), intent(in)           :: set   !< Searched set.
  logical,                   intent(in), optional :: back  !< Start of the last occurrence rather than the first.
  integer                                         :: i     !< Result of the search.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    i = verify(string=self%raw, set=set, back=back)
  else
    i = 0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sverify

  ! auxiliary methods
  elemental function basedir(self, sep)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the base directory name of a string containing a file name.
  !<
  !<### Example
  !<
  !<```fortran
  !< type(string) :: astring
  !< astring = '/bar/foo.tar.bz2'
  !< print '(A)', astring%basedir()//'' ! print "/bar"
  !<```
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self    !< The string.
  character(kind=CK, len=*), intent(in), optional :: sep     !< Directory separator.
  type(string)                                    :: basedir !< Base directory name.
  character(kind=CK, len=:), allocatable          :: sep_    !< Separator, default value.
  integer                                         :: pos     !< Character position.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    sep_ = UIX_DIR_SEP ; if (present(sep)) sep_ = sep
    basedir = self
    pos = index(self%raw, sep_, back=.true.)
    if (pos>0) basedir%raw = self%raw(1:pos-1)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction basedir

  elemental function basename(self, sep, extension, strip_last_extension)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the base file name of a string containing a file name.
  !<
  !< Optionally, the extension is also stripped if provided or the last one if required, e.g.
  !<
  !<### Example
  !<
  !<```fortran
  !< type(string) :: astring
  !< astring = 'bar/foo.tar.bz2'
  !< print '(A)', astring%basename(extension='.tar.bz2')//''        ! print "foo"
  !< print '(A)', astring%basename(strip_last_extension=.true.)//'' ! print "foo.tar"
  !<```
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self                 !< The string.
  character(kind=CK, len=*), intent(in), optional :: sep                  !< Directory separator.
  character(kind=CK, len=*), intent(in), optional :: extension            !< File extension.
  logical,                   intent(in), optional :: strip_last_extension !< Flag to enable the stripping of last extension.
  type(string)                                    :: basename             !< Base file name.
  character(kind=CK, len=:), allocatable          :: sep_                 !< Separator, default value.
  integer                                         :: pos                  !< Character position.
#ifdef __GFORTRAN__
  character(kind=CK, len=:), allocatable          :: temporary            !< Temporary storage, workaround for GNU bug.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    sep_ = UIX_DIR_SEP ; if (present(sep)) sep_ = sep
    basename = self
#ifdef __GFORTRAN__
    temporary = basename%raw
    pos = index(temporary, sep_, back=.true.)
    if (pos>0) basename%raw = temporary(pos+1:)
#else
    pos = index(basename%raw, sep_, back=.true.)
    if (pos>0) basename%raw = self%raw(pos+1:)
#endif
    if (present(extension)) then
#ifdef __GFORTRAN__
      temporary = basename%raw
      pos = index(temporary, extension, back=.true.)
      if (pos>0) basename%raw = temporary(1:pos-1)
#else
      pos = index(basename%raw, extension, back=.true.)
      if (pos>0) basename%raw = basename%raw(1:pos-1)
#endif
    elseif (present(strip_last_extension)) then
      if (strip_last_extension) then
#ifdef __GFORTRAN__
        temporary = basename%raw
        pos = index(temporary, '.', back=.true.)
        basename%raw = temporary(1:pos-1)
#else
        pos = index(basename%raw, '.', back=.true.)
        basename%raw = basename%raw(1:pos-1)
#endif
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction basename

  elemental function camelcase(self, sep)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string with all words capitalized without spaces.
  !<
  !< @note Multiple subsequent separators are collapsed to one occurence.
  !<
  !<### Example
  !<
  !<```fortran
  !< type(string) :: astring
  !< astring = 'caMeL caSe var'
  !< print '(A)', astring%camelcase()//'' ! print "CamelCaseVar"
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self      !< The string.
  character(kind=CK, len=*), intent(in), optional :: sep       !< Separator.
  type(string)                                    :: camelcase !< Camel case string.
  type(string), allocatable                       :: tokens(:) !< String tokens.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    call self%split(tokens=tokens, sep=sep)
    tokens = tokens%capitalize()
    camelcase = camelcase%join(array=tokens)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction camelcase

  elemental function capitalize(self) result(capitalized)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string with its first character capitalized and the rest lowercased.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self        !< The string.
  type(string)              :: capitalized !< Upper case string.
  integer                   :: c           !< Character counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    capitalized = self%lower()
    c = index(LOWER_ALPHABET, capitalized%raw(1:1))
    if (c>0) capitalized%raw(1:1) = UPPER_ALPHABET(c:c)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction capitalize

  pure function chars(self) result(raw)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the raw characters data.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)              :: self !< The string.
  character(kind=CK, len=:), allocatable :: raw  !< Raw characters data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    raw = self%raw
  else
    raw = ''
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction chars

  elemental function decode(self, codec) result(decoded)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string decoded accordingly the codec.
  !<
  !< @note Only BASE64 codec is currently available.
  !<
  !<### Example
  !<
  !<```fortran
  !< type(string) :: astring
  !< astring = 'SG93IGFyZSB5b3U/'
  !< print '(A)', astring%decode(codec='base64')//'' ! print "How are you?"
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in) :: self    !< The string.
  character(kind=CK, len=*), intent(in) :: codec   !< Encoding codec.
  type(string)                          :: decoded !< Decoded string.
  type(string)                          :: codec_u !< Encoding codec in upper case string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    decoded = self
    codec_u = codec
    select case(codec_u%upper()//'')
    case('BASE64')
      call b64_decode(code=self%raw, s=decoded%raw)
    endselect
    decoded = decoded%strip(remove_nulls=.true.)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction decode

  elemental function encode(self, codec) result(encoded)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string encoded accordingly the codec.
  !<
  !< @note Only BASE64 codec is currently available.
  !<
  !<### Example
  !<
  !<```fortran
  !< type(string) :: astring
  !< astring = 'How are you?'
  !< print '(A)', astring%encode(codec='base64')//'' ! print "SG93IGFyZSB5b3U/"
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in) :: self    !< The string.
  character(kind=CK, len=*), intent(in) :: codec   !< Encoding codec.
  type(string)                          :: encoded !< Encoded string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    encoded = codec
    select case(encoded%upper()//'')
    case('BASE64')
      call b64_encode(s=self%raw, code=encoded%raw)
    endselect
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction encode

  elemental function escape(self, to_escape, esc) result(escaped)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Escape backslashes (or custom escape character).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self      !< The string.
  character(kind=CK, len=1), intent(in)           :: to_escape !< Character to be escaped.
  character(kind=CK, len=*), intent(in), optional :: esc       !< Character used to escape.
  type(string)                                    :: escaped   !< Escaped string.
  character(kind=CK, len=:), allocatable          :: esc_      !< Character to escape, local variable.
  integer                                         :: c         !< Character counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    esc_ = BACKSLASH ; if (present(esc)) esc_ = esc
    escaped%raw = ''
    do c=1, len(self%raw)
      if (self%raw(c:c)==to_escape) then
        escaped%raw = escaped%raw//esc_//to_escape
      else
        escaped%raw = escaped%raw//self%raw(c:c)
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction escape

  elemental function extension(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the extension of a string containing a file name.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)              :: self      !< The string.
  type(string)                           :: extension !< Extension file name.
  integer                                :: pos       !< Character position.
#ifdef __GFORTRAN__
  character(kind=CK, len=:), allocatable :: temporary !< Temporary storage, workaround for GNU bug.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    extension = ''
    pos = index(self%raw, '.', back=.true.)
#ifdef __GFORTRAN__
    temporary = self%raw
    if (pos>0) extension%raw = temporary(pos:)
#else
    if (pos>0) extension%raw = self%raw(pos:)
#endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction extension

  elemental function fill(self, width, right, filling_char) result(filled)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Pad string on the left (or right) with zeros (or other char) to fill width.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self          !< The string.
  integer,                   intent(in)           :: width         !< Final width of filled string.
  logical,                   intent(in), optional :: right         !< Fill on the right instead of left.
  character(kind=CK, len=1), intent(in), optional :: filling_char  !< Filling character (default "0").
  type(string)                                    :: filled        !< Filled string.
  logical                                         :: right_        !< Fill on the right instead of left, local variable.
  character(kind=CK, len=1)                       :: filling_char_ !< Filling character (default "0"), local variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    if (width>len(self%raw)) then
      right_ = .false. ; if (present(right)) right_ = right
      filling_char_ = '0' ; if (present(filling_char)) filling_char_ = filling_char
      if (.not.right_) then
        filled%raw = repeat(filling_char_, width-len(self%raw))//self%raw
      else
        filled%raw = self%raw//repeat(filling_char_, width-len(self%raw))
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction fill

  elemental subroutine free(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Free dynamic memory.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(inout) :: self !< The string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) deallocate(self%raw)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  elemental function insert_character(self, substring, pos) result(inserted)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Insert substring into string at a specified position.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),    intent(in) :: self      !< The string.
  character(len=*), intent(in) :: substring !< Substring.
  integer,          intent(in) :: pos       !< Position from which insert substring.
  type(string)                 :: inserted  !< Inserted string.
  integer                      :: safepos   !< Safe position from which insert substring.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    inserted = self
    safepos = min(max(1, pos), len(self%raw))
    if (safepos==1) then
      inserted%raw = substring//self%raw
    elseif (safepos==len(self%raw)) then
      inserted%raw = self%raw//substring
    else
      inserted%raw = self%raw(1:safepos-1)//substring//self%raw(safepos:)
    endif
  else
    inserted%raw = substring
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction insert_character

  elemental function insert_string(self, substring, pos) result(inserted)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Insert substring into string at a specified position.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self      !< The string.
  type(string),  intent(in) :: substring !< Substring.
  integer,       intent(in) :: pos       !< Position from which insert substring.
  type(string)              :: inserted  !< Inserted string.
  integer                   :: safepos   !< Safe position from which insert substring.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    inserted = self
    if (allocated(substring%raw)) then
      safepos = min(max(1, pos), len(self%raw))
      if (safepos==1) then
        inserted%raw = substring%raw//self%raw
      elseif (safepos==len(self%raw)) then
        inserted%raw = self%raw//substring%raw
      else
        inserted%raw = self%raw(1:safepos-1)//substring%raw//self%raw(safepos:)
      endif
    endif
  else
    if (allocated(substring%raw)) inserted%raw = substring%raw
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction insert_string

  pure function join_strings(self, array, sep) result(join)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string that is a join of an array of strings.
  !<
  !< The join-separator is set equals to self if self has a value or it is set to a null string ''. This value can be overridden
  !< passing a custom separator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self      !< The string.
  type(string),              intent(in)           :: array(1:) !< Array to be joined.
  character(kind=CK, len=*), intent(in), optional :: sep       !< Separator.
  type(string)                                    :: join      !< The join of array.
  character(kind=CK, len=:), allocatable          :: sep_      !< Separator, default value.
  integer                                         :: a         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    sep_ = self%raw
  else
    sep_ = ''
  endif
  if (present(sep)) sep_ = sep
  join = ''
  do a=2, size(array, dim=1)
    if (allocated(array(a)%raw)) join%raw = join%raw//sep_//array(a)%raw
  enddo
  if (allocated(array(1)%raw)) then
    join%raw = array(1)%raw//join%raw
  else
    join%raw = join%raw(len(sep_)+1:len(join%raw))
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction join_strings

  pure function join_characters(self, array, sep) result(join)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string that is a join of an array of characters.
  !<
  !< The join-separator is set equals to self if self has a value or it is set to a null string ''. This value can be overridden
  !< passing a custom separator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self      !< The string.
  character(kind=CK, len=*), intent(in)           :: array(1:) !< Array to be joined.
  character(kind=CK, len=*), intent(in), optional :: sep       !< Separator.
  type(string)                                    :: join      !< The join of array.
  character(kind=CK, len=:), allocatable          :: sep_      !< Separator, default value.
  integer                                         :: a         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    sep_ = self%raw
  else
    sep_ = ''
  endif
  if (present(sep)) sep_ = sep
  join = ''
  do a=2, size(array, dim=1)
    if (array(a)/='') join%raw = join%raw//sep_//array(a)
  enddo
  if (array(1)/='') then
    join%raw = array(1)//join%raw
  else
    join%raw = join%raw(len(sep_)+1:len(join%raw))
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction join_characters

  elemental function lower(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string with all lowercase characters.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self  !< The string.
  type(string)              :: lower !< Upper case string.
  integer                   :: n1    !< Characters counter.
  integer                   :: n2    !< Characters counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    lower = self
    do n1=1, len(self%raw)
      n2 = index(UPPER_ALPHABET, self%raw(n1:n1))
      if (n2>0) lower%raw(n1:n1) = LOWER_ALPHABET(n2:n2)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction lower

  pure function partition(self, sep) result(partitions)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Split string at separator and return the 3 parts (before, the separator and after).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self            !< The string.
  character(kind=CK, len=*), intent(in), optional :: sep             !< Separator.
  type(string)                                    :: partitions(1:3) !< Partions: before the separator, the separator itsels and
                                                                     !< after the separator.
  character(kind=CK, len=:), allocatable          :: sep_            !< Separator, default value.
  integer                                         :: c               !< Character counter.
#ifdef __GFORTRAN__
  character(kind=CK, len=:), allocatable          :: temporary       !< Temporary storage, workaround for GNU bug.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    sep_ = SPACE ; if (present(sep)) sep_ = sep

    partitions(1) = self
    partitions(2) = sep_
    partitions(3) = ''
    if (len(sep_)>=len(self%raw)) return
    c = index(self%raw, sep_)
    if (c>0) then
#ifdef __GFORTRAN__
      temporary = self%raw
      partitions(1)%raw = temporary(1:c-1)
      partitions(2)%raw = temporary(c:c+len(sep_)-1)
      partitions(3)%raw = temporary(c+len(sep_):)
#else
      partitions(1)%raw = self%raw(1:c-1)
      partitions(2)%raw = self%raw(c:c+len(sep_)-1)
      partitions(3)%raw = self%raw(c+len(sep_):)
#endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction partition

  subroutine read_file(self, file, form, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Read a file as a single string stream.
  !<
  !< @note All the lines are stored into the string self as a single ascii stream. Each line (record) is separated by a `new_line`
  !< character.
  !<
  !< @note For unformatted read only `access='stream'` is supported with new_line as line terminator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),    intent(inout)           :: self       !< The string.
  character(len=*), intent(in)              :: file       !< File name.
  character(len=*), intent(in),    optional :: form       !< Format of unit.
  integer,          intent(out),   optional :: iostat     !< IO status code.
  character(len=*), intent(inout), optional :: iomsg      !< IO status message.
  type(string)                              :: form_      !< Format of unit, local variable.
  integer                                   :: iostat_    !< IO status code, local variable.
  character(len=:), allocatable             :: iomsg_     !< IO status message, local variable.
  integer                                   :: unit       !< Logical unit.
  logical                                   :: does_exist !< Check if file exist.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  iomsg_ = repeat(' ', 99) ; if (present(iomsg)) iomsg_ = iomsg
  inquire(file=file, iomsg=iomsg_, iostat=iostat_, exist=does_exist)
  if (does_exist) then
    form_ = 'FORMATTED' ; if (present(form)) form_ = form ; form_ = form_%upper()
    select case(form_%chars())
    case('FORMATTED')
      open(newunit=unit, file=file, status='OLD', action='READ', iomsg=iomsg_, iostat=iostat_, err=10)
    case('UNFORMATTED')
      open(newunit=unit, file=file, status='OLD', action='READ', form='UNFORMATTED', access='STREAM', &
           iomsg=iomsg_, iostat=iostat_, err=10)
    endselect
    call self%read_lines(unit=unit, form=form, iomsg=iomsg_, iostat=iostat_)
    10 close(unit)
  endif
  if (present(iostat)) iostat = iostat_
  if (present(iomsg)) iomsg = iomsg_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_file

  subroutine read_line(self, unit, form, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Read line (record) from a connected unit.
  !<
  !< The line is read as an ascii stream read until the eor is reached.
  !<
  !< @note For unformatted read only `access='stream'` is supported with new_line as line terminator.
  !---------------------------------------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only : iostat_eor
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),    intent(inout)           :: self    !< The string.
  integer,          intent(in)              :: unit    !< Logical unit.
  character(len=*), intent(in),    optional :: form    !< Format of unit.
  integer,          intent(out),   optional :: iostat  !< IO status code.
  character(len=*), intent(inout), optional :: iomsg   !< IO status message.
  type(string)                              :: form_   !< Format of unit, local variable.
  integer                                   :: iostat_ !< IO status code, local variable.
  character(len=:),          allocatable    :: iomsg_  !< IO status message, local variable.
  character(kind=CK, len=:), allocatable    :: line    !< Line storage.
  character(kind=CK, len=1)                 :: ch      !< Character storage.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  form_ = 'FORMATTED' ; if (present(form)) form_ = form ; form_ = form_%upper()
  iomsg_ = repeat(' ', 99) ; if (present(iomsg)) iomsg_ = iomsg
  line = ''
  select case(form_%chars())
  case('FORMATTED')
    do
      read(unit, "(A)", advance='no', iostat=iostat_, iomsg=iomsg_, err=10, end=10, eor=10) ch
      line = line//ch
    enddo
  case('UNFORMATTED')
    do
      read(unit, iostat=iostat_, iomsg=iomsg_, err=10, end=10) ch
      if (ch==new_line('a')) then
        iostat_ = iostat_eor
        exit
      endif
      line = line//ch
    enddo
  endselect
  10 if (line/='') self%raw = line
  if (present(iostat)) iostat = iostat_
  if (present(iomsg)) iomsg = iomsg_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_line

  subroutine read_lines(self, unit, form, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Read (all) lines (records) from a connected unit as a single ascii stream.
  !<
  !< @note All the lines are stored into the string self as a single ascii stream. Each line (record) is separated by a `new_line`
  !< character. The line is read as an ascii stream read until the eor is reached.
  !<
  !< @note The connected unit is rewinded. At a successful exit current record is at eof, at the beginning otherwise.
  !<
  !< @note For unformatted read only `access='stream'` is supported with new_line as line terminator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),    intent(inout)           :: self    !< The string.
  integer,          intent(in)              :: unit    !< Logical unit.
  character(len=*), intent(in),    optional :: form    !< Format of unit.
  integer,          intent(out),   optional :: iostat  !< IO status code.
  character(len=*), intent(inout), optional :: iomsg   !< IO status message.
  integer                                   :: iostat_ !< IO status code, local variable.
  character(len=:), allocatable             :: iomsg_  !< IO status message, local variable.
  type(string)                              :: lines   !< Lines storage.
  type(string)                              :: line    !< Line storage.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  iomsg_ = repeat(' ', 99) ; if (present(iomsg)) iomsg_ = iomsg
  rewind(unit)
  iostat_ = 0
  lines%raw = ''
  do
    line%raw = ''
    call line%read_line(unit=unit, form=form, iostat=iostat_, iomsg=iomsg_)
    if (iostat_/=0.and..not.is_iostat_eor(iostat_)) then
      exit
    elseif (line/='') then
      lines%raw = lines%raw//line%raw//new_line('a')
    endif
  enddo
  if (lines%raw/='') self%raw = lines%raw
  if (present(iostat)) iostat = iostat_
  if (present(iomsg)) iomsg = iomsg_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_lines

  elemental function replace(self, old, new, count) result(replaced)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string with all occurrences of substring old replaced by new.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self      !< The string.
  character(kind=CK, len=*), intent(in)           :: old       !< Old substring.
  character(kind=CK, len=*), intent(in)           :: new       !< New substring.
  integer,                   intent(in), optional :: count     !< Number of old occurences to be replaced.
  type(string)                                    :: replaced  !< The string with old replaced by new.
  integer                                         :: r         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    replaced = self
    r = 0
    do
      if (index(replaced%raw, old)>0) then
        replaced = replaced%replace_one_occurrence(old=old, new=new)
        r = r + 1
        if (present(count)) then
          if (r>=count) exit
        endif
      else
        exit
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction replace

  elemental function reverse(self) result(reversed)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a reversed string.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self     !< The string.
  type(string)              :: reversed !< The reversed string.
  integer                   :: length   !< Length of the string.
  integer                   :: c        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    reversed = self
    length = len(self%raw)
    do c=1, length
      reversed%raw(c:c) = self%raw(length-c+1:length-c+1)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction reverse

  function search(self, tag_start, tag_end, in_string, in_character, istart, iend) result(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Search for *tagged* record into string, return the first record found (if any) matching the tags.
  !<
  !< Optionally, returns the indexes of tag start/end, thus this is not an `elemental` function.
  !<
  !< @note The tagged record is searched into self if allocated otherwise into `in_string` if passed or, eventually, into
  !< `in_character` is passed. If tag is not found the return string is not allocated and the start/end indexes (if requested) are
  !< zero.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)            :: self         !< The string.
  character(kind=CK, len=*), intent(in)            :: tag_start    !< Start tag.
  character(kind=CK, len=*), intent(in)            :: tag_end      !< End tag.
  type(string),              intent(in),  optional :: in_string    !< Search into this string.
  character(kind=CK, len=*), intent(in),  optional :: in_character !< Search into this character string.
  integer,                   intent(out), optional :: istart       !< Starting index of tag inside the string.
  integer,                   intent(out), optional :: iend         !< Ending index of tag inside the string.
  type(string)                                     :: tag          !< First tag found.
  character(kind=CK, len=:), allocatable           :: raw          !< Raw string into which search the tag.
  integer                                          :: istart_      !< Starting index of tag inside the string, local variable.
  integer                                          :: iend_        !< Ending index of tag inside the string, local variable.
  logical                                          :: found        !< Flag for inquiring search result.
  integer                                          :: nested_tags  !< Number of nested tags inside tag.
  integer                                          :: t            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  raw = ''
  if (present(in_string)) then
    raw = in_string%raw
  elseif (present(in_character)) then
    raw = in_character
  else
    if (allocated(self%raw)) then
      raw = self%raw
    endif
  endif
  istart_ = 0
  iend_ = 0
  if (raw/='') then
    found = .false.
    istart_ = index(raw, tag_start)
    iend_ = index(raw, tag_end)
    if (istart_>0.and.iend_>0) then
      iend_ = iend_ + len(tag_end) - 1
      tag%raw = raw(istart_:iend_)
      nested_tags = tag%count(tag_start)
      if (nested_tags>1) then
        do t=2, nested_tags
          iend_ = iend_ + len(tag_end) - 1 + index(raw(iend_+1:), tag_end)
        enddo
        tag%raw = raw(istart_:iend_)
      endif
    endif
  endif
  if (present(istart)) istart = istart_
  if (present(iend)) iend = iend_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction search

  pure function slice(self, istart, iend) result(raw)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return the raw characters data sliced.
  !<
  !<### Example
  !<
  !<```fortran
  !< type(string) :: astring        !< A string.
  !< astring = 'the Quick Brown fox Jumps over the Lazy Dog.'
  !< print "(A)", astring%slice(11,25) ! print "Brown fox Jumps"
  !<```
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)              :: self   !< The string.
  integer,       intent(in)              :: istart !< Slice start index.
  integer,       intent(in)              :: iend   !< Slice end   index.
  character(kind=CK, len=:), allocatable :: raw    !< Raw characters data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    raw = self%raw(istart:iend)
  else
    raw = ''
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction slice

  elemental function snakecase(self, sep)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string with all words lowercase separated by "_".
  !<
  !< @note Multiple subsequent separators are collapsed to one occurence.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self      !< The string.
  character(kind=CK, len=*), intent(in), optional :: sep       !< Separator.
  type(string)                                    :: snakecase !< Snake case string.
  type(string), allocatable                       :: tokens(:) !< String tokens.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    call self%split(tokens=tokens, sep=sep)
    tokens = tokens%lower()
    snakecase = snakecase%join(array=tokens, sep='_')
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction snakecase

  pure subroutine split(self, tokens, sep)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a list of substring in the string, using sep as the delimiter string.
  !<
  !< @note Multiple subsequent separators are collapsed to one occurence.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self           !< The string.
  character(kind=CK, len=*), intent(in), optional :: sep            !< Separator.
  type(string), allocatable, intent(out)          :: tokens(:)      !< Tokens substring.
  character(kind=CK, len=:), allocatable          :: sep_           !< Separator, default value.
  integer                                         :: No             !< Number of occurrences of sep.
  integer                                         :: t              !< Character counter.
  type(string)                                    :: temporary      !< Temporary storage.
  type(string), allocatable                       :: temp_toks(:,:) !< Temporary tokens substring.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    sep_ = SPACE ; if (present(sep)) sep_ = sep

    temporary = self%unique(sep_)
    No = temporary%count(sep_)
    allocate(temp_toks(3, No))
    temp_toks(:, 1) = temporary%partition(sep_)
    if (No>1) then
      do t=2, No
        temp_toks(:, t) = temp_toks(3, t-1)%partition(sep_)
      enddo
    endif
    if (temp_toks(1, 1)%raw/=''.and.temp_toks(3, No)%raw/='') then
      allocate(tokens(No+1))
      do t=1, No
        if (t==No) then
          tokens(t  ) = temp_toks(1, t)
          tokens(t+1) = temp_toks(3, t)
        else
          tokens(t) = temp_toks(1, t)
        endif
      enddo
    elseif (temp_toks(1, 1)%raw/='') then
      allocate(tokens(No))
      do t=1, No
        tokens(t) = temp_toks(1, t)
      enddo
    elseif (temp_toks(3, No)%raw/='') then
      allocate(tokens(No))
      do t=2, No
        if (t==No) then
          tokens(t-1) = temp_toks(1, t)
          tokens(t  ) = temp_toks(3, t)
        else
          tokens(t-1) = temp_toks(1, t)
        endif
      enddo
    else
      allocate(tokens(No-1))
      do t=2, No
        tokens(t-1) = temp_toks(1, t)
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine split

  elemental function startcase(self, sep)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string with all words capitalized, e.g. title case.
  !<
  !< @note Multiple subsequent separators are collapsed to one occurence.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self      !< The string.
  character(kind=CK, len=*), intent(in), optional :: sep       !< Separator.
  type(string)                                    :: startcase !< Start case string.
  character(kind=CK, len=:), allocatable          :: sep_      !< Separator, default value.
  type(string), allocatable                       :: tokens(:) !< String tokens.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    sep_ = SPACE ; if (present(sep)) sep_ = sep
    call self%split(tokens=tokens, sep=sep_)
    tokens = tokens%capitalize()
    startcase = startcase%join(array=tokens, sep=sep_)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction startcase

  elemental function strip(self, remove_nulls)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a copy of the string with the leading and trailing characters removed.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)           :: self         !< The string.
  logical,       intent(in), optional :: remove_nulls !< Remove null characters at the end.
  type(string)                        :: strip        !< The stripped string.
  integer                             :: c            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    strip = self%adjustl()
    strip = strip%trim()
    if (present(remove_nulls)) then
      if (remove_nulls) then
        c = index(self%raw, char(0))
        if (c>0) strip%raw = strip%raw(1:c-1)
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction strip

  elemental function swapcase(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a copy of the string with uppercase characters converted to lowercase and vice versa.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self     !< The string.
  type(string)              :: swapcase !< Upper case string.
  integer                   :: n1       !< Characters counter.
  integer                   :: n2       !< Characters counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    swapcase = self
    do n1=1, len(self%raw)
      n2 = index(UPPER_ALPHABET, self%raw(n1:n1))
      if (n2>0) then
        swapcase%raw(n1:n1) = LOWER_ALPHABET(n2:n2)
      else
        n2 = index(LOWER_ALPHABET, self%raw(n1:n1))
        if (n2>0) swapcase%raw(n1:n1) = UPPER_ALPHABET(n2:n2)
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction swapcase

  elemental function to_integer_I1P(self, kind) result(to_number)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cast string to integer (I1P).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self      !< The string.
  integer(I1P),  intent(in) :: kind      !< Mold parameter for kind detection.
  integer(I1P)              :: to_number !< The number into the string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    if (self%is_integer()) read(self%raw, *) to_number
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction to_integer_I1P

  elemental function to_integer_I2P(self, kind) result(to_number)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cast string to integer (I2P).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self      !< The string.
  integer(I2P),  intent(in) :: kind      !< Mold parameter for kind detection.
  integer(I2P)              :: to_number !< The number into the string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    if (self%is_integer()) read(self%raw, *) to_number
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction to_integer_I2P

  elemental function to_integer_I4P(self, kind) result(to_number)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cast string to integer (I4P).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self      !< The string.
  integer(I4P),  intent(in) :: kind      !< Mold parameter for kind detection.
  integer(I4P)              :: to_number !< The number into the string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    if (self%is_integer()) read(self%raw, *) to_number
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction to_integer_I4P

  elemental function to_integer_I8P(self, kind) result(to_number)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cast string to integer (I8P).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self      !< The string.
  integer(I8P),  intent(in) :: kind      !< Mold parameter for kind detection.
  integer(I8P)              :: to_number !< The number into the string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    if (self%is_integer()) read(self%raw, *) to_number
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction to_integer_I8P

  elemental function to_real_R4P(self, kind) result(to_number)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cast string to real (R4P).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self      !< The string.
  real(R4P),     intent(in) :: kind      !< Mold parameter for kind detection.
  real(R4P)                 :: to_number !< The number into the string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    if (self%is_real()) read(self%raw, *) to_number
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction to_real_R4P

  elemental function to_real_R8P(self, kind) result(to_number)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cast string to real (R8P).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self      !< The string.
  real(R8P),     intent(in) :: kind      !< Mold parameter for kind detection.
  real(R8P)                 :: to_number !< The number into the string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    if (self%is_real()) read(self%raw, *) to_number
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction to_real_R8P

  elemental function to_real_R16P(self, kind) result(to_number)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Cast string to real (R16P).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self      !< The string.
  real(R16P),    intent(in) :: kind      !< Mold parameter for kind detection.
  real(R16P)                :: to_number !< The number into the string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    if (self%is_real()) read(self%raw, *) to_number
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction to_real_R16P

  elemental function unescape(self, to_unescape, unesc) result(unescaped)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Unescape double backslashes (or custom escaped character).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self        !< The string.
  character(kind=CK, len=1), intent(in)           :: to_unescape !< Character to be unescaped.
  character(kind=CK, len=*), intent(in), optional :: unesc       !< Character used to unescape.
  type(string)                                    :: unescaped   !< Escaped string.
  character(kind=CK, len=:), allocatable          :: unesc_      !< Character to unescape, local variable.
  integer                                         :: c           !< Character counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    unesc_ = '' ; if (present(unesc)) unesc_ = unesc
    unescaped%raw = ''
    c = 1
    do
      if (c>len(self%raw)) exit
      if (c==len(self%raw)) then
        unescaped%raw = unescaped%raw//self%raw(c:c)
        exit
      else
        if (self%raw(c:c+1)==BACKSLASH//to_unescape) then
          unescaped%raw = unescaped%raw//to_unescape
          c = c + 2
        else
          unescaped%raw = unescaped%raw//self%raw(c:c)
          c = c + 1
        endif
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction unescape

  elemental function unique(self, substring) result(uniq)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Reduce to one (unique) multiple (sequential) occurrences of a substring into a string.
  !<
  !< For example the string ' ab-cre-cre-ab' is reduce to 'ab-cre-ab' if the substring is '-cre'.
  !< @note Eventual multiple trailing white space are not reduced to one occurrence.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self       !< The string.
  character(kind=CK, len=*), intent(in), optional :: substring  !< Substring which multiple occurences must be reduced to one.
  character(kind=CK, len=:), allocatable          :: substring_ !< Substring, default value.
  type(string)                                    :: uniq       !< String parsed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    substring_ = SPACE ; if (present(substring)) substring_ = substring

    uniq = self
    do
      if (.not.uniq%index(repeat(substring_, 2))>0) exit
      uniq = uniq%replace(old=repeat(substring_, 2), new=substring_)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction unique

  elemental function upper(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string with all uppercase characters.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self  !< The string.
  type(string)              :: upper !< Upper case string.
  integer                   :: n1    !< Characters counter.
  integer                   :: n2    !< Characters counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    upper = self
    do n1=1, len(self%raw)
      n2 = index(LOWER_ALPHABET, self%raw(n1:n1))
      if (n2>0) upper%raw(n1:n1) = UPPER_ALPHABET(n2:n2)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction upper

  subroutine write_file(self, file, form, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Write a single string stream into file.
  !<
  !< @note For unformatted read only `access='stream'` is supported with new_line as line terminator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),    intent(in)              :: self    !< The string.
  character(len=*), intent(in)              :: file    !< File name.
  character(len=*), intent(in),    optional :: form    !< Format of unit.
  integer,          intent(out),   optional :: iostat  !< IO status code.
  character(len=*), intent(inout), optional :: iomsg   !< IO status message.
  type(string)                              :: form_   !< Format of unit, local variable.
  integer                                   :: iostat_ !< IO status code, local variable.
  character(len=:), allocatable             :: iomsg_  !< IO status message, local variable.
  integer                                   :: unit    !< Logical unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  iomsg_ = repeat(' ', 99) ; if (present(iomsg)) iomsg_ = iomsg
  form_ = 'FORMATTED' ; if (present(form)) form_ = form ; form_ = form_%upper()
  select case(form_%chars())
  case('FORMATTED')
    open(newunit=unit, file=file, action='WRITE', iomsg=iomsg_, iostat=iostat_, err=10)
  case('UNFORMATTED')
    open(newunit=unit, file=file, action='WRITE', form='UNFORMATTED', access='STREAM', iomsg=iomsg_, iostat=iostat_, err=10)
  endselect
  call self%write_lines(unit=unit, form=form, iomsg=iomsg_, iostat=iostat_)
  10 close(unit)
  if (present(iostat)) iostat = iostat_
  if (present(iomsg)) iomsg = iomsg_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine write_file

  subroutine write_line(self, unit, form, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Write line (record) to a connected unit.
  !<
  !< @note If the connected unit is unformatted a `new_line()` character is added at the end (if necessary) to mark the end of line.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),    intent(in)              :: self    !< The string.
  integer,          intent(in)              :: unit    !< Logical unit.
  character(len=*), intent(in),    optional :: form    !< Format of unit.
  integer,          intent(out),   optional :: iostat  !< IO status code.
  character(len=*), intent(inout), optional :: iomsg   !< IO status message.
  type(string)                              :: form_   !< Format of unit, local variable.
  integer                                   :: iostat_ !< IO status code, local variable.
  character(len=:), allocatable             :: iomsg_  !< IO status message, local variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  iostat_ = 0
  iomsg_ = repeat(' ', 99) ; if (present(iomsg)) iomsg_ = iomsg
  if (allocated(self%raw)) then
    form_ = 'FORMATTED' ; if (present(form)) form_ = form ; form_ = form_%upper()
    select case(form_%chars())
    case('FORMATTED')
      write(unit, "(A)", iostat=iostat_, iomsg=iomsg_) self%raw
    case('UNFORMATTED')
      if (self%end_with(new_line('a'))) then
        write(unit, iostat=iostat_, iomsg=iomsg_) self%raw
      else
        write(unit, iostat=iostat_, iomsg=iomsg_) self%raw//new_line('a')
      endif
    endselect
  endif
  if (present(iostat)) iostat = iostat_
  if (present(iomsg)) iomsg = iomsg_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine write_line

  subroutine write_lines(self, unit, form, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Write lines (records) to a connected unit.
  !<
  !< This method checks if self contains more than one line (records) and writes them as lines (records).
  !<
  !< @note If the connected unit is unformatted a `new_line()` character is added at the end (if necessary) to mark the end of line.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),    intent(in)              :: self     !< The string.
  integer,          intent(in)              :: unit     !< Logical unit.
  character(len=*), intent(in),    optional :: form     !< Format of unit.
  integer,          intent(out),   optional :: iostat   !< IO status code.
  character(len=*), intent(inout), optional :: iomsg    !< IO status message.
  type(string), allocatable                 :: lines(:) !< Lines.
  integer                                   :: l        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    call self%split(tokens=lines, sep=new_line('a'))
    do l=1, size(lines, dim=1)
       call lines(l)%write_line(unit=unit, form=form, iostat=iostat, iomsg=iomsg)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine write_lines

  ! inquire
  elemental function end_with(self, suffix, start, end)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return true if a string ends with a specified suffix.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self     !< The string.
  character(kind=CK, len=*), intent(in)           :: suffix   !< Searched suffix.
  integer,                   intent(in), optional :: start    !< Start position into the string.
  integer,                   intent(in), optional :: end      !< End position into the string.
  logical                                         :: end_with !< Result of the test.
  integer                                         :: start_   !< Start position into the string, local variable.
  integer                                         :: end_     !< End position into the string, local variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  end_with = .false.
  if (allocated(self%raw)) then
    start_ = 1             ; if (present(start)) start_ = start
    end_   = len(self%raw) ; if (present(end))   end_   = end
    if (len(suffix)<=len(self%raw(start_:end_))) then
      end_with = index(self%raw(start_:end_), suffix)==(len(self%raw(start_:end_)) - len(suffix) + 1)
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction end_with

  elemental function is_allocated(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return true if the string is allocated.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self     !< The string.
  logical                   :: is_allocated !< Result of the test.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_allocated = allocated(self%raw)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_allocated

  elemental function is_digit(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return true if all characters in the string are digits.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self     !< The string.
  logical                   :: is_digit !< Result of the test.
  integer                   :: c        !< Character counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_digit = .false.
  if (allocated(self%raw)) then
    do c=1, len(self%raw)
      select case (self%raw(c:c))
      case ('0':'9')
        is_digit = .true.
      case default
        is_digit = .false.
        exit
      end select
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_digit

  elemental function is_integer(self, allow_spaces)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return true if the string contains an integer.
  !<
  !< The regular expression is `\s*[\+\-]?\d+([eE]\+?\d+)?\s*`. The parse algorithm is done in stages:
  !<
  !< | S0  | S1      | S2  | S3   | S4  | S5  | S6  |
  !< |-----|---------|-----|------|-----|-----|-----|
  !< |`\s*`|`[\+\-]?`|`\d+`|`[eE]`|`\+?`|`\d+`|`\s*`|
  !<
  !< Exit on stages-parsing results in:
  !<
  !< | S0 | S1 | S2 | S3 | S4 | S5 | S6 |
  !< |----|----|----|----|----|----|----|
  !< |  F |  F |  T |  F |  F |  T |  T |
  !<
  !< @note This implementation is courtesy of
  !< [tomedunn](https://github.com/tomedunn/fortran-string-utility-module/blob/master/src/string_utility_module.f90#L294)
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)           :: self          !< The string.
  logical,       intent(in), optional :: allow_spaces  !< Allow leading-trailing spaces.
  logical                             :: is_integer    !< Result of the test.
  logical                             :: allow_spaces_ !< Allow leading-trailing spaces, local variable.
  integer                             :: stage         !< Stages counter.
  integer                             :: c             !< Character counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    allow_spaces_ = .true. ; if (present(allow_spaces)) allow_spaces_ = allow_spaces
    stage = 0
    is_integer = .true.
    do c=1, len(self%raw)
      select case(self%raw(c:c))
      case(SPACE, TAB)
        select case(stage)
        case(0, 6)
          is_integer = allow_spaces_
        case(2, 5)
          is_integer = allow_spaces_
          stage = 6
        case default
          is_integer = .false.
        endselect
      case('-')
        select case(stage)
        case(0)
          stage = 1
        case default
          is_integer = .false.
        end select
      case('+')
        select case(stage)
        case(0)
          stage = 1
        case(3)
          stage = 4
        case default
          is_integer = .false.
        endselect
      case('0':'9')
        select case(stage)
        case(0:1)
          stage = 2
        case(3:4)
          stage = 5
        case default
          continue
        endselect
      case ('e','E')
        select case(stage)
        case(2)
          stage = 3
        case default
          is_integer = .false.
        endselect
      case default
        is_integer = .false.
      endselect
      if (.not.is_integer) exit
    enddo
  endif
  if (is_integer) then
    select case(stage)
    case(2, 5, 6)
      is_integer = .true.
    case default
      is_integer = .false.
    end select
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_integer

  elemental function is_lower(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return true if all characters in the string are lowercase.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self     !< The string.
  logical                   :: is_lower !< Result of the test.
  integer                   :: c        !< Character counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_lower = .false.
  if (allocated(self%raw)) then
    is_lower = .true.
    do c=1, len(self%raw)
      if (index(UPPER_ALPHABET, self%raw(c:c))>0) then
        is_lower = .false.
        exit
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_lower

  elemental function is_number(self, allow_spaces)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return true if the string contains a number (real or integer).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)           :: self         !< The string.
  logical,       intent(in), optional :: allow_spaces !< Allow leading-trailing spaces.
  logical                             :: is_number    !< Result of the test.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_number = (self%is_integer(allow_spaces=allow_spaces).or.self%is_real(allow_spaces=allow_spaces))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_number

  elemental function is_real(self, allow_spaces)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return true if the string contains a real.
  !<
  !< The regular expression is `\s*[\+\-]?\d*(|\.?\d*([deDE][\+\-]?\d+)?)\s*`. The parse algorithm is done in stages:
  !<
  !< | S0  | S1      | S2  | S3  | S4  | S5     | S6      | S7  | S8  |
  !< |-----|---------|-----|-----|-----|--------|---------|-----|-----|
  !< |`\s*`|`[\+\-]?`|`\d*`|`\.?`|`\d*`|`[deDE]`|`[\+\-]?`|`\d*`|`\s*`|
  !<
  !< Exit on stages-parsing results in:
  !<
  !< | S0 | S1 | S2 | S3 | S4 | S5 | S6 | S7 | S8 |
  !< |----|----|----|----|----|----|----|----|----|
  !  |  F |  F |  T |  T |  T |  F |  F |  T |  T |
  !<
  !< @note This implementation is courtesy of
  !< [tomedunn](https://github.com/tomedunn/fortran-string-utility-module/blob/master/src/string_utility_module.f90#L614)
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)           :: self              !< The string.
  logical,       intent(in), optional :: allow_spaces      !< Allow leading-trailing spaces.
  logical                             :: is_real           !< Result of the test.
  logical                             :: allow_spaces_     !< Allow leading-trailing spaces, local variable.
  logical                             :: has_leading_digit !< Check the presence of leading digits.
  integer                             :: stage             !< Stages counter.
  integer                             :: c                 !< Character counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    allow_spaces_ = .true. ; if (present(allow_spaces)) allow_spaces_ = allow_spaces
    stage = 0
    is_real = .true.
    has_leading_digit = .false.
    do c=1, len(self%raw)
      select case(self%raw(c:c))
      case(SPACE, TAB)
        select case(stage)
        case(0, 8)
          is_real = allow_spaces_
          continue
        case(2:4, 7)
          is_real = allow_spaces_
          stage = 8
        case default
          is_real = .false.
        endselect
      case('+', '-')
        select case(stage)
        case(0)
          stage = 1
        case(5)
          stage = 6
        case default
          is_real = .false.
        endselect
      case('0':'9')
        select case(stage)
        case(0:1)
          stage = 2
          has_leading_digit = .true.
        case(3)
          stage = 4
        case(5:6)
          stage = 7
        case default
          continue
        endselect
      case('.')
        select case(stage)
        case(0:2)
          stage = 3
        case default
          is_real = .false.
        endselect
      case('e','E','d','D')
        select case(stage)
        case(2:4)
          stage = 5
        case default
          is_real = .false.
        endselect
      case default
        is_real = .false.
      endselect
      if (.not.is_real) exit
    enddo
  endif
  if (is_real) then
    select case(stage)
    case(2, 4, 7, 8)
      is_real = .true.
    case(3)
      is_real = has_leading_digit
    case default
      is_real = .false.
    endselect
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_real

  elemental function is_upper(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return true if all characters in the string are uppercase.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: self     !< The string.
  logical                   :: is_upper !< Result of the test.
  integer                   :: c        !< Character counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_upper = .false.
  if (allocated(self%raw)) then
    is_upper = .true.
    do c=1, len(self%raw)
      if (index(LOWER_ALPHABET, self%raw(c:c))>0) then
        is_upper = .false.
        exit
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_upper

  elemental function start_with(self, prefix, start, end)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return true if a string starts with a specified prefix.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)           :: self       !< The string.
  character(kind=CK, len=*), intent(in)           :: prefix     !< Searched prefix.
  integer,                   intent(in), optional :: start      !< Start position into the string.
  integer,                   intent(in), optional :: end        !< End position into the string.
  logical                                         :: start_with !< Result of the test.
  integer                                         :: start_     !< Start position into the string, local variable.
  integer                                         :: end_       !< End position into the string, local variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  start_with = .false.
  if (allocated(self%raw)) then
    start_ = 1             ; if (present(start)) start_ = start
    end_   = len(self%raw) ; if (present(end))   end_   = end
    if (len(prefix)<=len(self%raw(start_:end_))) then
      start_with = index(self%raw(start_:end_), prefix)==1
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction start_with

  ! private methods

  ! assignments
  elemental subroutine string_assign_string(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assignment operator from string input.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(inout) :: lhs !< Left hand side.
  type(string),  intent(in)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(rhs%raw)) lhs%raw = rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine string_assign_string

  elemental subroutine string_assign_character(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assignment operator from character input.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(inout) :: lhs !< Left hand side.
  character(kind=CK, len=*), intent(in)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lhs%raw = rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine string_assign_character

  elemental subroutine string_assign_integer_I1P(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assignment operator from real input.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(inout) :: lhs !< Left hand side.
  integer(I1P),  intent(in)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lhs%raw = trim(str(rhs))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine string_assign_integer_I1P

  elemental subroutine string_assign_integer_I2P(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assignment operator from real input.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(inout) :: lhs !< Left hand side.
  integer(I2P),  intent(in)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lhs%raw = trim(str(rhs))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine string_assign_integer_I2P

  elemental subroutine string_assign_integer_I4P(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assignment operator from real input.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(inout) :: lhs !< Left hand side.
  integer(I4P),  intent(in)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lhs%raw = trim(str(rhs))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine string_assign_integer_I4P

  elemental subroutine string_assign_integer_I8P(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assignment operator from real input.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(inout) :: lhs !< Left hand side.
  integer(I8P),  intent(in)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lhs%raw = trim(str(rhs))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine string_assign_integer_I8P

  elemental subroutine string_assign_real_R4P(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assignment operator from real input.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(inout) :: lhs !< Left hand side.
  real(R4P),     intent(in)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lhs%raw = trim(str(rhs))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine string_assign_real_R4P

  elemental subroutine string_assign_real_R8P(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assignment operator from real input.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(inout) :: lhs !< Left hand side.
  real(R8P),     intent(in)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lhs%raw = trim(str(rhs))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine string_assign_real_R8P

  elemental subroutine string_assign_real_R16P(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Assignment operator from real input.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(inout) :: lhs !< Left hand side.
  real(R16P),    intent(in)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lhs%raw = trim(str(rhs))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine string_assign_real_R16P

  ! contatenation operators
  pure function string_concat_string(lhs, rhs) result(concat)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Concatenation with string.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)              :: lhs    !< Left hand side.
  type(string),  intent(in)              :: rhs    !< Right hand side.
  character(kind=CK, len=:), allocatable :: concat !< Concatenated string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  concat = ''
  if (allocated(lhs%raw)) concat = lhs%raw
  if (allocated(rhs%raw)) concat = concat//rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_concat_string

  pure function string_concat_character(lhs, rhs) result(concat)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Concatenation with character.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)  :: lhs    !< Left hand side.
  character(kind=CK, len=*), intent(in)  :: rhs    !< Right hand side.
  character(kind=CK, len=:), allocatable :: concat !< Concatenated string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(lhs%raw)) then
    concat = lhs%raw//rhs
  else
    concat = rhs
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_concat_character

  pure function character_concat_string(lhs, rhs) result(concat)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Concatenation with character (inverted).
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in)  :: lhs    !< Left hand side.
  class(string),             intent(in)  :: rhs    !< Right hand side.
  character(kind=CK, len=:), allocatable :: concat !< Concatenated string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(rhs%raw)) then
    concat = lhs//rhs%raw
  else
    concat = lhs
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction character_concat_string

  elemental function string_concat_string_string(lhs, rhs) result(concat)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Concatenation with string.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in)              :: lhs       !< Left hand side.
  type(string),  intent(in)              :: rhs       !< Right hand side.
  type(string)                           :: concat    !< Concatenated string.
  character(kind=CK, len=:), allocatable :: temporary !< Temporary concatenated string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  temporary = ''
  if (allocated(lhs%raw)) temporary = lhs%raw
  if (allocated(rhs%raw)) temporary = temporary//rhs%raw
  if (temporary/='') concat%raw = temporary
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_concat_string_string

  elemental function string_concat_character_string(lhs, rhs) result(concat)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Concatenation with character.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)  :: lhs    !< Left hand side.
  character(kind=CK, len=*), intent(in)  :: rhs    !< Right hand side.
  type(string)                           :: concat !< Concatenated string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(lhs%raw)) then
    concat%raw = lhs%raw//rhs
  else
    concat%raw = rhs
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_concat_character_string

  elemental function character_concat_string_string(lhs, rhs) result(concat)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Concatenation with character (inverted).
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in)  :: lhs    !< Left hand side.
  class(string),             intent(in)  :: rhs    !< Right hand side.
  type(string)                           :: concat !< Concatenated string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(rhs%raw)) then
    concat%raw = lhs//rhs%raw
  else
    concat%raw = lhs
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction character_concat_string_string

  ! logical operators
  elemental function string_eq_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Equal to string logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: lhs   !< Left hand side.
  type(string),  intent(in) :: rhs   !< Right hand side.
  logical                   :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw == rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_eq_string

  elemental function string_eq_character(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Equal to character logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in) :: lhs   !< Left hand side.
  character(kind=CK, len=*), intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw == rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_eq_character

  elemental function character_eq_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Equal to character (inverted) logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in) :: lhs   !< Left hand side.
  class(string),             intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = rhs%raw == lhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction character_eq_string

  elemental function string_ne_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Not equal to string logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: lhs   !< Left hand side.
  type(string),  intent(in) :: rhs   !< Right hand side.
  logical                   :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw /= rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_ne_string

  elemental function string_ne_character(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Not equal to character logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in) :: lhs   !< Left hand side.
  character(kind=CK, len=*), intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw /= rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_ne_character

  elemental function character_ne_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Not equal to character (inverted) logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in) :: lhs   !< Left hand side.
  class(string),             intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = rhs%raw /= lhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction character_ne_string

  elemental function string_lt_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Lower than to string logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: lhs   !< Left hand side.
  type(string),  intent(in) :: rhs   !< Right hand side.
  logical                   :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw < rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_lt_string

  elemental function string_lt_character(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Lower than to character logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in) :: lhs   !< Left hand side.
  character(kind=CK, len=*), intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw < rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_lt_character

  elemental function character_lt_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Lower than to character (inverted) logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in) :: lhs   !< Left hand side.
  class(string),             intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs < rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction character_lt_string

  elemental function string_le_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Lower equal than to string logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: lhs   !< Left hand side.
  type(string),  intent(in) :: rhs   !< Right hand side.
  logical                   :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw <= rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_le_string

  elemental function string_le_character(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Lower equal than to character logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in) :: lhs   !< Left hand side.
  character(kind=CK, len=*), intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw <= rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_le_character

  elemental function character_le_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Lower equal than to character (inverted) logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in) :: lhs   !< Left hand side.
  class(string),             intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs <= rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction character_le_string

  elemental function string_ge_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Greater equal than to string logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: lhs   !< Left hand side.
  type(string),  intent(in) :: rhs   !< Right hand side.
  logical                   :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw >= rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_ge_string

  elemental function string_ge_character(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Greater equal than to character logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in) :: lhs   !< Left hand side.
  character(kind=CK, len=*), intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw >= rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_ge_character

  elemental function character_ge_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Greater equal than to character (inverted) logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in) :: lhs   !< Left hand side.
  class(string),             intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs >= rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction character_ge_string

  elemental function string_gt_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Greater than to string logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string), intent(in) :: lhs   !< Left hand side.
  type(string),  intent(in) :: rhs   !< Right hand side.
  logical                   :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw > rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_gt_string

  elemental function string_gt_character(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Greater than to character logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in) :: lhs   !< Left hand side.
  character(kind=CK, len=*), intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs%raw > rhs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_gt_character

  elemental function character_gt_string(lhs, rhs) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Greater than to character (inverted) logical operator.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(kind=CK, len=*), intent(in) :: lhs   !< Left hand side.
  class(string),             intent(in) :: rhs   !< Right hand side.
  logical                               :: is_it !< Opreator test result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = lhs > rhs%raw
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction character_gt_string

  ! IO
  subroutine read_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Formatted input.
  !<
  !< @bug Change temporary acks: find a more precise length of the input string and avoid the trimming!
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(inout) :: dtv         !< The string.
  integer,                   intent(in)    :: unit        !< Logical unit.
  character(len=*),          intent(in)    :: iotype      !< Edit descriptor.
  integer,                   intent(in)    :: v_list(:)   !< Edit descriptor list.
  integer,                   intent(out)   :: iostat      !< IO status code.
  character(len=*),          intent(inout) :: iomsg       !< IO status message.
  character(len=len(iomsg))                :: local_iomsg !< Local variant of iomsg, so it doesn't get inappropriately redefined.
  character(kind=CK, len=1)                :: delim       !< String delimiter, if any.
  character(kind=CK, len=100)              :: temporary   !< Temporary storage string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (iotype == 'LISTDIRECTED') then
    call get_next_non_blank_character_any_record(unit=unit, ch=delim, iostat=iostat, iomsg=iomsg)
    if (iostat/=0) return
    if (delim=='"'.OR.delim=="'") then
      call dtv%read_delimited(unit=unit, delim=delim, iostat=iostat, iomsg=local_iomsg)
    else
      ! step back before the non-blank
      read(unit, "(TL1)", iostat=iostat, iomsg=iomsg)
      if (iostat /= 0) return
      call dtv%read_undelimited_listdirected(unit=unit, iostat=iostat, iomsg=local_iomsg)
    endif
    if (is_iostat_eor(iostat)) then
      ! suppress IOSTAT_EOR
      iostat = 0
    elseif (iostat /= 0) then
      iomsg = local_iomsg
    endif
    return
  else
    read(unit, "(A)", iostat=iostat, iomsg=iomsg)temporary
    dtv%raw = trim(temporary)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_formatted

  subroutine read_delimited(dtv, unit, delim, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Read a delimited string from a unit connected for formatted input.
  !<
  !< If the closing delimiter is followed by end of record, then we return end of record.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(out)   :: dtv       !< The string.
  integer,                   intent(in)    :: unit      !< Logical unit.
  character(kind=CK, len=1), intent(in)    :: delim     !< String delimiter.
  integer,                   intent(out)   :: iostat    !< IO status code.
  character(kind=CK, len=*), intent(inout) :: iomsg     !< IO status message.
  character(kind=CK, len=1)                :: ch        !< A character read.
  logical                                  :: was_delim !< Indicates that the last character read was a delimiter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  was_delim = .false.
  dtv%raw = ''
  do
    read(unit, "(A)", iostat=iostat, iomsg=iomsg) ch
    if (is_iostat_eor(iostat)) then
      if (was_delim) then
        ! end of delimited string followed by end of record is end of the string. Pass back the end of record condition to the
        ! caller
        return
      else
        ! end of record without terminating delimiter - move along
        cycle
      endif
    elseif (iostat /= 0) THEN
      return
    endif
    if (ch == delim) then
      if (was_delim) then
        ! doubled delimiter is one delimiter in the value
        dtv%raw = dtv%raw // ch
        was_delim = .false.
      else
        ! need to test next character to see what is happening
        was_delim = .true.
      endif
    elseif (was_delim) then
      ! the previous character was actually the delimiter for the end of the string. Put back this character
      read(unit, "(TL1)", iostat=iostat, iomsg=iomsg)
      return
    else
      dtv%raw = dtv%raw // ch
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_delimited

  subroutine read_undelimited_listdirected(dtv, unit, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Read an undelimited (no leading apostrophe or double quote) character value according to the rules for list directed input.
  !<
  !< A blank, comma/semicolon (depending on the decimal mode), slash or end of record terminates the string.
  !<
  !< If input is terminated by end of record, then this procedure returns an end-of-record condition.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),    intent(inout) :: dtv           !< The string.
  integer,          intent(in)    :: unit          !< Logical unit.
  integer,          intent(out)   :: iostat        !< IO status code.
  character(len=*), intent(inout) :: iomsg         !< IO status message.
  logical                         :: decimal_point !<True if DECIMAL=POINT in effect.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call get_decimal_mode(unit=unit, decimal_point=decimal_point, iostat=iostat, iomsg=iomsg)
  if (iostat /= 0) return
  call dtv%read_undelimited(unit=unit, terminators=' '//'/'//merge(CK_',', CK_';', decimal_point), iostat=iostat, iomsg=iomsg)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_undelimited_listdirected

  subroutine read_undelimited(dtv, unit, terminators, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Read an undelimited string up until end of record or a character from a set of terminators is encountered.
  !<
  !< If a terminator is encountered, the file position will be at that terminating character. If end of record is encountered, the
  !< file remains at end of record.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(inout) :: dtv         !< The string.
  integer,                   intent(in)    :: unit        !< Logical unit.
  character(kind=CK, len=*), intent(in)    :: terminators !< Characters that are considered to terminate the string.
                                                          !< Blanks in this string are meaningful.
  integer,                   intent(out)   :: iostat      !< IO status code.
  character(len=*),          intent(inout) :: iomsg       !< IO status message.
  character(kind=CK, len=1)                :: ch          !< A character read.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  dtv%raw = ''
  do
    read(unit, "(A)", iostat=iostat, iomsg=iomsg) ch
    if (is_iostat_eor(iostat)) then
      ! end of record just means end of string. We pass on the condition
      return
    elseif (iostat /= 0) then
      ! something odd happened
      return
    endif
    if (scan(ch, terminators) /= 0) then
      ! change the file position so that the next read sees the terminator
      read(unit, "(TL1)", iostat=iostat, iomsg=iomsg)
      if (iostat /= 0) return
      iostat = 0
      return
    endif
    ! we got a character - append it
    dtv%raw = dtv%raw // ch
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_undelimited

  subroutine write_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Formatted output.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)    :: dtv       !< The string.
  integer,                   intent(in)    :: unit      !< Logical unit.
  character(kind=CK, len=*), intent(in)    :: iotype    !< Edit descriptor.
  integer,                   intent(in)    :: v_list(:) !< Edit descriptor list.
  integer,                   intent(out)   :: iostat    !< IO status code.
  character(kind=CK, len=*), intent(inout) :: iomsg     !< IO status message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(dtv%raw)) then
    write(unit, "(A)", iostat=iostat, iomsg=iomsg)dtv%raw
  else
    write(unit, "(A)", iostat=iostat, iomsg=iomsg)''
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine write_formatted

  subroutine read_unformatted(dtv, unit, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Unformatted input.
  !<
  !< @bug Change temporary acks: find a more precise length of the input string and avoid the trimming!
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(inout) :: dtv       !< The string.
  integer,                   intent(in)    :: unit      !< Logical unit.
  integer,                   intent(out)   :: iostat    !< IO status code.
  character(kind=CK, len=*), intent(inout) :: iomsg     !< IO status message.
  character(kind=CK, len=100)              :: temporary !< Temporary storage string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit, iostat=iostat, iomsg=iomsg)temporary
  dtv%raw = trim(temporary)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_unformatted

  subroutine write_unformatted(dtv, unit, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Unformatted output.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)    :: dtv    !< The string.
  integer,                   intent(in)    :: unit   !< Logical unit.
  integer,                   intent(out)   :: iostat !< IO status code.
  character(kind=CK, len=*), intent(inout) :: iomsg  !< IO status message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(dtv%raw)) then
    write(unit, iostat=iostat, iomsg=iomsg)dtv%raw
  else
    write(unit, iostat=iostat, iomsg=iomsg)''
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine write_unformatted

  ! miscellanea
  elemental function replace_one_occurrence(self, old, new) result(replaced)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string with the first occurrence of substring old replaced by new.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(string),             intent(in)  :: self      !< The string.
  character(kind=CK, len=*), intent(in)  :: old       !< Old substring.
  character(kind=CK, len=*), intent(in)  :: new       !< New substring.
  type(string)                           :: replaced  !< The string with old replaced by new.
  integer                                :: pos       !< Position from which replace old.
#ifdef __GFORTRAN__
  character(kind=CK, len=:), allocatable :: temporary !< Temporary storage, workaround for GNU bug.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%raw)) then
    replaced = self
    pos = index(string=self%raw, substring=old)
    if (pos>0) then
#ifdef __GFORTRAN__
      temporary = self%raw
      if (pos==1) then
        replaced%raw = new//temporary(len(old)+1:)
      else
        replaced%raw = temporary(1:pos-1)//new//temporary(pos+len(old):)
      endif
#else
      if (pos==1) then
        replaced%raw = new//self%raw(len(old)+1:)
      else
        replaced%raw = self%raw(1:pos-1)//new//self%raw(pos+len(old):)
      endif
#endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction replace_one_occurrence

  ! non type-bound-procedures
  subroutine get_delimiter_mode(unit, delim, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Get the DELIM changeable connection mode for the given unit.
  !<
  !< If the unit is connected to an internal file, then the default value of NONE is always returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only : iostat_inquire_internal_unit
  !---------------------------------------------------------------------------------------------------------------------------------
  integer,                   intent(in)    :: unit         !< The unit for the connection.
  character(len=1, kind=CK), intent(out)   :: delim        !< Represents the value of the DELIM mode.
  integer,                   intent(out)   :: iostat       !< IOSTAT error code, non-zero on error.
  character(*),              intent(inout) :: iomsg        !< IOMSG explanatory message - only defined if iostat is non-zero.
  character(10)                            :: delim_buffer !< Buffer for INQUIRE about DELIM, sized for APOSTROHPE.
  character(len(iomsg))                    :: local_iomsg  !< Local variant of iomsg, so it doesn't get inappropriately redefined.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! get the string representation of the changeable mode
  inquire(unit, delim=delim_buffer, iostat=iostat, iomsg=local_iomsg)
  if (iostat == iostat_inquire_internal_unit) then
    ! no way of determining the DELIM mode for an internal file
    iostat = 0
    delim = ''
    return
  elseif (iostat /= 0) then
    iomsg = local_iomsg
    return
  endif
  ! interpret the DELIM string
  if (delim_buffer == 'QUOTE') then
    delim = '"'
  elseif (delim_buffer == 'APOSTROPHE') then
    delim = ''''
  else
    delim = '"'
  endif
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_delimiter_mode

  subroutine get_next_non_blank_character_this_record(unit, ch, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Get the next non-blank character in the current record.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer,                   intent(in)    :: unit   !< Logical unit.
  character(kind=CK, len=1), intent(out)   :: ch     !< The non-blank character read. Not valid if IOSTAT is non-zero.
  integer,                   intent(out)   :: iostat !< IO status code.
  character(kind=CK, len=*), intent(inout) :: iomsg  !< IO status message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do
    ! we spcify non-advancing, just in case we want this callable outside the context of a child input statement
    ! the PAD specifier simply saves the need for the READ statement to define ch if EOR is hit
    ! read(unit, "(A)", iostat=iostat, iomsg=iomsg, advance='NO') ch
    ! ...but that causes ifort to blow up at runtime
    read(unit, "(A)", iostat=iostat, iomsg=iomsg, pad='NO') ch
    if (iostat /= 0) return
    if (ch /= '') exit
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_next_non_blank_character_this_record

  subroutine get_next_non_blank_character_any_record(unit, ch, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Get the next non-blank character, advancing records if necessary.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer,                   intent(in)    :: unit        !< Logical unit.
  character(kind=CK, len=1), intent(out)   :: ch          !< The non-blank character read. Not valid if IOSTAT is non-zero.
  integer,                   intent(out)   :: iostat      !< IO status code.
  character(kind=CK, len=*), intent(inout) :: iomsg       !< IO status message.
  character(len(iomsg))                    :: local_iomsg !< Local variant of iomsg, so it doesn't get inappropriately redefined.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do
    call get_next_non_blank_character_this_record(unit=unit, ch=ch, iostat=iostat, iomsg=local_iomsg)
    if (is_iostat_eor(iostat)) then
      ! try again on the next record
      read (unit, "(/)", iostat=iostat, iomsg=iomsg)
      if (iostat /= 0) return
    elseif (iostat /= 0) then
      ! some sort of problem
      iomsg = local_iomsg
      return
    else
      ! got it
      exit
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_next_non_blank_character_any_record

  subroutine get_decimal_mode(unit, decimal_point, iostat, iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Get the DECIMAL changeable connection mode for the given unit.
  !<
  !< If the unit is connected to an internal file, then the default value of DECIMAL is always returned. This may not be the
  !< actual value in force at the time of the call to this procedure.
  !---------------------------------------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only : iostat_inquire_internal_unit
  !---------------------------------------------------------------------------------------------------------------------------------
  integer,                   intent(in)    :: unit           !< Logical unit.
  logical,                   intent(out)   :: decimal_point  !< True if the decimal mode is POINT, false otherwise.
  integer,                   intent(out)   :: iostat         !< IO status code.
  character(kind=CK, len=*), intent(inout) :: iomsg          !< IO status message.
  character(5)                             :: decimal_buffer !< Buffer for INQUIRE about DECIMAL, sized for POINT or COMMA.
  character(len(iomsg))                    :: local_iomsg    !< Local variant of iomsg, so it doesn't get inappropriately redefined.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(unit, decimal=decimal_buffer, iostat=iostat, iomsg=local_iomsg)
  if (iostat == iostat_inquire_internal_unit) then
    ! no way of determining the decimal mode for an internal file
    iostat = 0
    decimal_point = .true.
    return
  else if (iostat /= 0) then
    iomsg = local_iomsg
    return
  endif
  decimal_point = decimal_buffer == 'POINT'
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_decimal_mode
endmodule stringifor_string_t
