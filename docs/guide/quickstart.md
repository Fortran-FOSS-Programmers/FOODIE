---
title: Quick start
---

# Quick start

Using FOODIE takes two steps:

1. define your system as a **concrete integrand**, a type extending the abstract `integrand_object`;
2. create an **integrator** (directly, or by name through `foodie_integrator_factory`) and call its `integrate` method in
   your time loop.

This page integrates the linear constant coefficients equation

$$
U_t = R(U) = a U + b \quad \Rightarrow \quad U(t) = \left(U_0 + \frac{b}{a}\right) e^{a(t - t_0)} - \frac{b}{a}
$$

used by the FOODIE tests. The code is the integrand `integrand_lcce` of `src/tests/tester/foodie_test_integrand_lcce.f90`
reduced to the procedures required by `integrand_object` (the tests version extends `integrand_tester_object`, which adds
the exact solution, the I/O and the command line handling) and the time loop of
`src/tests/tester/foodie_test_integrate.f90`. The complete program, [at the end of this page](#the-complete-program), is
compiled and run against the library.

## Define the integrand

A concrete integrand holds the state and the parameters of the system, and implements the deferred procedures of
`integrand_object`:

```fortran
use foodie, only : integrand_object
use penf, only : I_P, R_P

type, extends(integrand_object) :: integrand_lcce
   !< The linear constant coefficient equation field, U_t = a*U + b.
   real(R_P) :: a=0._R_P  !< *a* constant.
   real(R_P) :: b=0._R_P  !< *b* constant.
   real(R_P) :: U=0._R_P  !< Integrand (state) variable.
   contains
      ! integrand_object deferred methods
      procedure, pass(self), public :: description         !< Return an informative description of the integrand.
      procedure, pass(self), public :: integrand_dimension !< Return integrand dimension.
      procedure, pass(self), public :: t => dU_dt          !< Time derivative, residuals.
      ! operators
      procedure, pass(lhs), public :: local_error !<`||integrand_lcce - integrand_lcce||` operator.
      ! +
      procedure, pass(lhs), public :: integrand_add_integrand !< `+` operator.
      procedure, pass(lhs), public :: integrand_add_real      !< `+ real` operator.
      procedure, pass(rhs), public :: real_add_integrand      !< `real +` operator.
      ! *
      procedure, pass(lhs), public :: integrand_multiply_integrand   !< `*` operator.
      procedure, pass(lhs), public :: integrand_multiply_real        !< `* real` operator.
      procedure, pass(rhs), public :: real_multiply_integrand        !< `real *` operator.
      procedure, pass(lhs), public :: integrand_multiply_real_scalar !< `* real_scalar` operator.
      procedure, pass(rhs), public :: real_scalar_multiply_integrand !< `real_scalar *` operator.
      ! -
      procedure, pass(lhs), public :: integrand_sub_integrand !< `-` operator.
      procedure, pass(lhs), public :: integrand_sub_real      !< `- real` operator.
      procedure, pass(rhs), public :: real_sub_integrand      !< `real -` operator.
      ! =
      procedure, pass(lhs), public :: assign_integrand !< `=` operator.
      procedure, pass(lhs), public :: assign_real      !< `= real` operator.
endtype integrand_lcce
```

The deferred procedures are:

| Binding | Interface | Meaning |
|---|---|---|
| `description` | `function (self, prefix) result(desc)`, `character(len=:), allocatable :: desc` | informative description |
| `integrand_dimension` | `function (self) result(integrand_dimension)`, `integer(I_P)` | size of the state |
| `t` | `function (self, t) result(dState_dt)`, `real(R_P), allocatable :: dState_dt(:)` | residual $R(t, U)$ |
| `local_error` (`.lterror.`) | `function (lhs, rhs) result(error)`, `real(R_P) :: error` | local truncation error estimate between two integrands |
| `integrand_add_integrand`, `integrand_sub_integrand`, `integrand_multiply_integrand` | `function (lhs, rhs) result(opr)`, both `class(integrand_object)` | `U + V`, `U - V`, `U * V` |
| `integrand_add_real`, `integrand_sub_real`, `integrand_multiply_real` | `rhs` is `real(R_P) :: rhs(1:)` | `U + r(:)`, `U - r(:)`, `U * r(:)` |
| `real_add_integrand`, `real_sub_integrand`, `real_multiply_integrand` | `lhs` is `real(R_P) :: lhs(1:)` | `r(:) + U`, `r(:) - U`, `r(:) * U` |
| `integrand_multiply_real_scalar`, `real_scalar_multiply_integrand` | the real is a scalar `real(R_P)` | `U * s`, `s * U` |
| `assign_integrand`, `assign_real` | `subroutine (lhs, rhs)`, `lhs` is `intent(inout)` | `U = V`, `U = r(:)` |

Every operator returns `real(R_P), allocatable :: opr(:)`, the state as a plain array: the expression
`U + (U%t(t=t) * Dt)` is evaluated as arrays and assigned back to `U` by `assign_real`. The operators, the assignments,
`description` and `integrand_dimension` are `pure` (see [Pure or impure integrands](./installation#pure-or-impure-integrands));
`t` and `local_error` need not be.

The residual function returns the time derivative of the state, as an array:

```fortran
pure function dU_dt(self, t) result(dState_dt)
!< Time derivative of field.
class(integrand_lcce), intent(in)           :: self         !< Integrand.
real(R_P),             intent(in), optional :: t            !< Time.
real(R_P), allocatable                      :: dState_dt(:) !< Integrand time derivative.

dState_dt = [self%a * self%U + self%b]
endfunction dU_dt
```

The symmetric operators receive the other operand as `class(integrand_object)` and must recover its dynamic type:

```fortran
pure function integrand_add_integrand(lhs, rhs) result(opr)
!< `+` operator.
class(integrand_lcce),   intent(in) :: lhs    !< Left hand side.
class(integrand_object), intent(in) :: rhs    !< Right hand side.
real(R_P), allocatable              :: opr(:) !< Operator result.

select type(rhs)
class is(integrand_lcce)
  opr = [lhs%U + rhs%U]
endselect
endfunction integrand_add_integrand
```

The two assignments copy a whole integrand, or overwrite the state with an array:

```fortran
pure subroutine assign_integrand(lhs, rhs)
!< `=` operator.
class(integrand_lcce),   intent(inout) :: lhs !< Left hand side.
class(integrand_object), intent(in)    :: rhs !< Right hand side.

select type(rhs)
class is(integrand_lcce)
   lhs%a  = rhs%a
   lhs%b  = rhs%b
   lhs%U  = rhs%U
endselect
endsubroutine assign_integrand

pure subroutine assign_real(lhs, rhs)
!< `= real` operator.
class(integrand_lcce), intent(inout) :: lhs     !< Left hand side.
real(R_P),             intent(in)    :: rhs(1:) !< Right hand side.

lhs%U = rhs(1)
endsubroutine assign_real
```

The other operators follow the same pattern; for a state of many variables the array returned by the operators is the
state flattened, and `assign_real` reshapes it back (see the [1D Euler example](./example-euler-1d)).

## Integrate with a multistage scheme

The factory returns an initialized integrator given the name of the scheme:

```fortran
subroutine foodie_integrator_factory(scheme, integrator, stages, tolerance, nu, alpha, iterations, autoupdate, U)
character(*),                          intent(in)  :: scheme     !< Selected integrator given.
class(integrator_object), allocatable, intent(out) :: integrator !< The FOODIE integrator.
integer(I_P),            optional,     intent(in)  :: stages     !< Stages of multi-stage methods.
real(R_P),               optional,     intent(in)  :: tolerance  !< Tolerance on the local truncation error.
real(R_P),               optional,     intent(in)  :: nu         !< Williams-Robert-Asselin filter coefficient.
real(R_P),               optional,     intent(in)  :: alpha      !< Robert-Asselin filter coefficient.
integer(I_P),            optional,     intent(in)  :: iterations !< Implicit iterations.
logical,                 optional,     intent(in)  :: autoupdate !< Enable cyclic autoupdate for multistep.
class(integrand_object), optional,     intent(in)  :: U          !< Integrand molding prototype.
```

`U` is a prototype used to allocate the integrator internal registers (stages, previous steps) with the dynamic type of
your integrand: always pass it. The optional arguments are used only by the classes they apply to: `stages` by the linear
SSP Runge-Kutta, `tolerance` by the embedded Runge-Kutta, `nu` and `alpha` by the filtered leapfrog, `iterations` by the
implicit schemes (fixed point iterations), `autoupdate` by the multistep schemes. An unknown scheme name stops the
program with an error message.

The returned integrator is polymorphic: a `select type` on the abstract multistage class gives access to its `integrate`
method:

```fortran
use foodie, only : foodie_integrator_factory, integrator_object, integrator_multistage_object
...
type(integrand_lcce)                  :: integrand   !< The integrand.
class(integrator_object), allocatable :: integrator  !< The integrator.
real(R_P), parameter                  :: Dt=0.01_R_P !< Time step.
real(R_P)                             :: time        !< Time.
integer(I_P)                          :: step        !< Time steps counter.

integrand%a = -1._R_P ; integrand%b = 0._R_P ; integrand%U = 1._R_P
time = 0._R_P

call foodie_integrator_factory(scheme='runge_kutta_ssp_stages_3_order_3', integrator=integrator, U=integrand)
select type(integrator)
class is(integrator_multistage_object)
   do step=1, 100
      call integrator%integrate(U=integrand, Dt=Dt, t=time)
      time = time + Dt
   enddo
endselect
```

The integrator does not compute the time step: `Dt` is always given by the caller (from a CFL condition, for a PDE). The
embedded Runge-Kutta schemes can reduce it to meet their tolerance: pass the optional `new_Dt` argument to know the step
actually used and advance `time` by it.

### Select a scheme by name

Being a string, the scheme can come from an input file or the command line: the same executable runs any of the
[schemes of the catalogue](./schemes). `foodie_integrator_schemes()` returns all their names,
`foodie_integrator_schemes(class_name='adams_bashforth')` those of one class, and `is_available(scheme)` checks a name
before using it.

### Use a concrete integrator directly

The factory is a convenience: every concrete integrator can be declared and initialized directly, e.g.

```fortran
use foodie, only : integrator_runge_kutta_ssp
...
type(integrator_runge_kutta_ssp) :: integrator

call integrator%initialize(scheme='runge_kutta_ssp_stages_3_order_3', U=integrand)
call integrator%integrate(U=integrand, Dt=Dt, t=time)
```

## Integrate with a multistep scheme

A multistep scheme of $N_s$ steps computes $U^{n+N_s}$ from the solutions at the previous steps, which the integrator
stores in its members `previous(1:steps)`, with their times `t(1:steps)` and time steps `Dt(1:steps)`. The first $N_s$
steps (the start-up) must be computed by other means, for example by a one-step scheme of adequate order; then, with
`autoupdate=.true.` (the default), every call of `integrate` shifts the stored steps by itself. The argument `t` of
`integrate` is the time of the last stored step:

```fortran
call foodie_integrator_factory(scheme='runge_kutta_ssp_stages_3_order_3', integrator=starter, U=integrand)
call foodie_integrator_factory(scheme='adams_bashforth_3', integrator=integrator, U=integrand)
select type(integrator)
class is(integrator_multistep_object)
   do step=1, 100
      if (step <= integrator%steps_number()) then
         select type(starter)
         class is(integrator_multistage_object)
            call starter%integrate(U=integrand, Dt=Dt, t=time)
         endselect
         time = time + Dt
         integrator%previous(step) = integrand
         integrator%Dt(step) = Dt
         integrator%t(step) = time
      else
         call integrator%integrate(U=integrand, Dt=Dt, t=time)
         time = time + Dt
      endif
   enddo
endselect
```

The tests suite starts the multistep schemes with the exact solution instead, so that the measured error is due to the
scheme alone. The multistage-multistep class (`integrator_multistage_multistep_object`) is started in the same way.

## Fast mode

Every integrator also provides `integrate_fast`, with the same arguments as `integrate`. It is written with in-place
procedures instead of the operators, so it does not allocate a temporary array at every operation:

```fortran
call integrator%integrate_fast(U=integrand, Dt=Dt, t=time)
```

The in-place procedures have default implementations in `integrand_object` that **do nothing**: an integrand used in
fast mode must override all of them, as `integrand_lcce` of the tests does:

| Binding | Interface | Meaning |
|---|---|---|
| `t_fast` | `subroutine (self, t)`, `self` is `intent(inout)` | `self = R(t, self)` |
| `integrand_add_integrand_fast` (`add_fast`) | `subroutine (opr, lhs, rhs)` | `opr = lhs + rhs` |
| `integrand_subtract_integrand_fast` (`subtract_fast`) | `subroutine (opr, lhs, rhs)` | `opr = lhs - rhs` |
| `integrand_multiply_integrand_fast` (`multiply_fast`) | `subroutine (opr, lhs, rhs)` | `opr = lhs * rhs` |
| `integrand_multiply_real_scalar_fast` (`multiply_fast`) | `subroutine (opr, lhs, rhs)`, `rhs` is `real(R_P)` | `opr = lhs * rhs` |

`integrator%has_fast_mode()` tells whether an integrator class provides the fast mode (all the current ones do).

## The complete program

Compiled with `gfortran -I lib/mod quickstart.f90 lib/libfoodie.a` after `fobis build --mode foodie-static-gnu`, it
prints the solution at $t=1$ of $U_t = -U$, $U(0)=1$, first by the 3rd order SSP Runge-Kutta, then by the 3 steps
Adams-Bashforth scheme, next to the exact value $e^{-1}$:

```text
U(1) =   3.678794257199937E-01, exact   3.678794411714423E-01
U(1) =   3.678793054498171E-01, exact   3.678794411714423E-01
```

<<< @/examples/quickstart.f90
