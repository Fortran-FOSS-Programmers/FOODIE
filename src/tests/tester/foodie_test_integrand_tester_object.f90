!< Define [[integrand_tester_object]], the abstract tester integrand.

module foodie_test_integrand_tester_object
!< Define [[integrand_tester_object]], the abstract tester integrand.

use flap, only : command_line_interface
use foodie, only : integrand_object
use penf, only : FR_P, R_P, I_P, str

implicit none
private
public :: integrand_tester_object

type, abstract, extends(integrand_object) :: integrand_tester_object
   !< The abstract tester integrand.
   !<
   !< This abstract provided some auxiliary methods useful for the tester machinery.
   contains
      procedure(description_interface),    pass(self), deferred :: description    !< Return an informative description of the test.
      procedure(exact_solution_interface), pass(self), deferred :: exact_solution !< Return exact solution.
      procedure(export_tecplot_interface), pass(self), deferred :: export_tecplot !< Export integrand to Tecplot file.
      procedure(parse_cli_interface),      pass(self), deferred :: parse_cli      !< Initialize from command line interface.
      procedure(set_cli_interface),        nopass,     deferred :: set_cli        !< Set command line interface.
endtype integrand_tester_object

abstract interface
   !< Abstract interfaces of [[integrand_tester_object]] class.
   pure function description_interface(self, prefix) result(desc)
   !< Return informative integrator description.
   import :: integrand_tester_object
   class(integrand_tester_object), intent(in)           :: self   !< Integrand.
   character(*),                   intent(in), optional :: prefix !< Prefixing string.
   character(len=:), allocatable                        :: desc   !< Description.
   endfunction description_interface

   pure function exact_solution_interface(self, t, U0) result(exact)
   !< Return exact solution.
   import :: integrand_object, integrand_tester_object, R_P
   class(integrand_tester_object), intent(in)           :: self     !< Integrand.
   real(R_P),                      intent(in)           :: t        !< Time.
   class(integrand_object),        intent(in), optional :: U0       !< Initial conditions.
   real(R_P), allocatable                               :: exact(:) !< Exact solution.
   endfunction exact_solution_interface

   subroutine export_tecplot_interface(self, file_name, t, scheme, close_file)
   !< Export integrand to Tecplot file.
   import :: integrand_tester_object, R_P
   class(integrand_tester_object), intent(in)           :: self       !< Advection field.
   character(*),                   intent(in), optional :: file_name  !< File name.
   real(R_P),                      intent(in), optional :: t          !< Time.
   character(*),                   intent(in), optional :: scheme     !< Scheme used to integrate integrand.
   logical,                        intent(in), optional :: close_file !< Flag for closing file.
   endsubroutine export_tecplot_interface

   subroutine parse_cli_interface(self, cli)
   !< Initialize from command line interface.
   import :: command_line_interface, integrand_tester_object
   class(integrand_tester_object), intent(inout) :: self !< Advection field.
   type(command_line_interface),   intent(inout) :: cli  !< Command line interface handler.
   endsubroutine parse_cli_interface

   subroutine set_cli_interface(cli)
   !< Set command line interface.
   import :: command_line_interface
   type(command_line_interface), intent(inout) :: cli !< Command line interface handler.
   endsubroutine set_cli_interface
endinterface
endmodule foodie_test_integrand_tester_object
