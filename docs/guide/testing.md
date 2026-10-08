---
title: Testing
---

# Testing

The order of accuracy is the property that defines a scheme: a coefficient mistyped in a tableau or a stage evaluated at
the wrong time leaves a scheme stable and plausible, but of lower order. FOODIE verifies every scheme at three independent
levels, all run by one command:

```bash
fobis rule --ex maketests
```

which

1. builds the tests suite with `fobis build --mode tests-gnu`: every program found under `src/` (excluding the
   directories and files listed in the `fobos` file) is a test, built into `exe/`;
2. runs it with `./scripts/run_tests.sh`: every executable of `exe/` must exit with status 0;
3. audits the Runge-Kutta tableaus with `python3 scripts/audit_rk_tableaus.py`: its exit status is the number of failed
   tableaus.

The tests suite is made of the programs of `src/tests/suite/`, built on the integrands and the integration driver of
`src/tests/tester/` (shared with the [CLI tester](#the-cli-tester)). Both test programs loop over all the schemes returned
by `foodie_integrator_schemes()`, so a new scheme is tested as soon as it is added to its class, and both test the
standard (`integrate`) and the fast (`integrate_fast`) mode of every scheme.

The expected orders are encoded once in `src/tests/suite/foodie_test_suite_utils.f90`: they are the orders in the
scheme names, with the documented exceptions of `runge_kutta_emd_stages_7_order_4` (5th order), `leapfrog_raw`
(1st order) and the linear SSP schemes (see [Notes on the orders](./schemes#notes-on-the-orders)).

## Layer 1: polynomial exactness

`foodie_test_polynomial_exactness` integrates the quadrature

$$ U_t = q\, t^{q-1} \quad \Rightarrow \quad U(t) = U_0 + t^q - t_0^q $$

whose exact solution is a polynomial of degree $q$, with $\Delta t = 0.1$ up to $t = 2$ (20 steps). A scheme of order
$p$ must integrate it exactly, to a relative error below $10^{-12}$, for every $q \le p$; a linear multistep scheme (Adams,
BDF, LMM-SSP, unfiltered leapfrog) must also **not** integrate it exactly for $q = p+1$ (relative error above
$10^{-10}$), its error constant being not null: this verifies that the claimed order is not underestimated.

For the linear multistep schemes polynomial exactness up to degree $p$ is equivalent to the order conditions; for the
multistage schemes it verifies the quadrature conditions and the correct use of the stage times. The residual does not
depend on $U$, so errors are not propagated (no stability constraint): this is why the test verifies also the highest
order schemes (Adams-Bashforth 16, Feagin 10), out of reach of the convergence test.

```text
PASS adams_bashforth_1        exact up to degree 1, max rel error +0.0
PASS adams_bashforth_1 [fast] exact up to degree 1, max rel error +0.0
PASS adams_bashforth_2        exact up to degree 2, max rel error +2.220446049250313E-16
...
PASS runge_kutta_emd_stages_7_order_4        exact up to degree 5, max rel error +2.220446049250313E-16
...
0 failures
```

## Layer 2: convergence order

`foodie_test_convergence_order` measures the observed order of every scheme. Each scheme integrates its convergence
problem with 11 time steps, halved from the coarsest one; the observed order is computed on the finest couple of
solutions whose errors are both in the asymptotic range $[10^{-12}, 10^{-3}]$, above round-off:

$$ p_{obs} = \frac{\log(e_k / e_{k+1})}{\log 2} $$

The test fails if $p_{obs} < p - 0.3$ or if no couple of solutions is in range. Two choices make the measure reliable:

- **the problem is nonlinear and non-autonomous.** The default problem is the Riccati equation

  $$ U_t = a\, t\, U^2 \quad \Rightarrow \quad U(t) = \frac{1}{\frac{1}{U_0} - \frac{a}{2}\left(t^2 - t_0^2\right)} $$

  with $a = -2$, $U_0 = 1$ (solution $1/(1+t^2)$), integrated up to $t = 2$ from $\Delta t = 0.5$. A linear autonomous
  problem verifies only the linear order conditions of the multistage schemes: a Runge-Kutta scheme can be 4th order on
  $U_t = aU$ and 2nd order on a nonlinear problem. The explicit dependence on $t$ also verifies the time at which every
  stage evaluates the residual;
- **the error is the maximum over all time steps** (L-inf norm in time), not the error at the final time: the global
  error of a scheme can change sign and vanish at the final time, giving a spurious order.

Two families are tested on other problems: the leapfrog schemes on the oscillation equations (the unfiltered leapfrog is
unstable on dissipative problems), the linear SSP schemes on the linear constant coefficients equation (their order holds
only for linear autonomous problems). Orders higher than 6 are reported (`SKIP`) but not asserted: in double precision
the asymptotic range is not reached before round-off. The multistep schemes are started with the exact solution, so the
measured error is due to the scheme alone.

```text
PASS adams_bashforth_4        on riccati: expected order 4, observed +3.9983684933355974 (errors: ...)
...
PASS leapfrog        on oscillation: expected order 2, observed +1.9998593342748432 (errors: ...)
PASS leapfrog_raw        on oscillation: expected order 1, observed +1.253738851620288 (errors: ...)
...
PASS runge_kutta_emd_stages_7_order_4        on riccati: expected order 5, observed +5.328893017300245 (errors: ...)
SKIP runge_kutta_emd_stages_17_order_10        on riccati: expected order 10, observed +13.591738570426527 (errors: ...)
...
0 failures
```

## Layer 3: audit of the Runge-Kutta tableaus

`scripts/audit_rk_tableaus.py` parses the coefficients of the Runge-Kutta schemes from the Fortran sources (decimal
literals, or `real(p_I8P)/real(q_I8P)` rationals) and evaluates, in 60 digits arithmetic, the order conditions of all the
rooted trees up to the claimed order $p$: for every tree $\tau$ with at most $p$ nodes

$$ \sum_i b_i\, \Phi_i(\tau) = \frac{1}{\gamma(\tau)} $$

where $\Phi_i(\tau)$ are the elementary weights of the tableau and $\gamma(\tau)$ the density of the tree. A tableau
passes if every residual is below $10^{-14}$ and the abscissae are the row sums of the matrix, $c_i = \sum_j a_{ij}$.
The low storage (Williamson 2N) schemes are converted to their Butcher form; both formulas of each embedded pair are
audited, the advancing one and the error estimator. The script uses only the Python standard library.

```text
PASS runge_kutta_ssp_stages_3_order_3                 order  3: max |order residual| = 1.00e-60, max |c - row sum| = 0.00e+00
PASS runge_kutta_emd_stages_7_order_4 [beta(:,1)]     order  5: max |order residual| = 4.10e-59, max |c - row sum| = 4.80e-59
PASS runge_kutta_emd_stages_7_order_4 [beta(:,2)]     order  4: max |order residual| = 2.40e-59, max |c - row sum| = 4.80e-59
PASS runge_kutta_emd_stages_17_order_10 [beta(:,1)]   order 10: max |order residual| = 2.00e-21, max |c - row sum| = 2.00e-21
PASS runge_kutta_emd_stages_17_order_10 [beta(:,2)]   order  8: max |order residual| = 2.00e-21, max |c - row sum| = 2.00e-21
...
0 failures
```

This is the only layer that verifies the full nonlinear order conditions of the high order Runge-Kutta schemes (Calvo 6,
Feagin 10), whose orders cannot be measured in double precision.

## The tests runner

`scripts/run_tests.sh` runs every executable of `exe/` and checks its exit status (an executable named `*_xfail_*` must
fail, one named `*mpi*` runs under `mpirun`); `--vmem KB` caps the virtual memory of every test. The test programs print
their reports and `stop 1` on failure.

## Coverage

```bash
fobis rule --ex makecoverage
```

builds the tests suite with coverage instrumentation (mode `tests-gnu-coverage`), runs it and computes the line coverage
of the FOODIE sources with `gcov` (`scripts/compute-coverage.sh`, which also writes `docs/public/coverage.json`).

## The CLI tester

The CLI tester, `src/tests/tester/foodie_tester.f90`, integrates any scheme (or all the schemes of a class, or all the
schemes) on one of five problems with the time steps given on the command line, and prints the errors and the observed
orders. It is not part of the tests suite; it is built with the `_IMPURE_` macro because of its WenOOF based linear
advection problem:

```bash
fobis build --mode tester-gnu
./build/tester/foodie_tester --help
```

| Problem | Equation | Options |
|---|---|---|
| `lcce` | linear constant coefficients, $U_t = aU + b$ | `-a`, `-b`, `-U0` |
| `linear_advection` | 1D linear advection with constant speed $a$, finite volumes with WENO reconstruction | `--w-scheme`, `--weno-order`, `--cfl`, `-a`, `--Ni`, `--initial_state` (`sin_wave`, `square_wave`), ... |
| `oscillation` | inertial oscillations, $v_{1,t} = -f v_2$, $v_{2,t} = f v_1$ | `--frequency` (default $10^{-4}$), `--U0` |
| `polynomial` | quadrature with polynomial solution, $U_t = q t^{q-1}$ | `--degree` (`-q`), `--U0` |
| `riccati` | nonlinear non-autonomous, $U_t = a t U^2$ | `-a`, `--U0` |

The general settings precede the problem name: `--scheme` (`-s`, a scheme, a class of schemes or `all`), `--time_step`
(`-Dt`, one or more time steps), `--final_time` (`-ft`), `--iterations` (implicit schemes), `--stages` (linear SSP),
`--fast`, `--save_results` (`-r`, Tecplot ASCII files), `--output`, `--save_frequency`. For example:

```bash
./build/tester/foodie_tester test -s runge_kutta_ssp_stages_3_order_3 -Dt 0.1 0.05 0.025 riccati
```

```text
runge_kutta_ssp_stages_3_order_3
riccati
 Dt =   0.10000000000000001      , error =    1.0368510165353895E-004
riccati
 Dt =    5.0000000000000003E-002 , error =    1.3015239291869207E-005
 Observed order =  2.99
riccati
 Dt =    2.5000000000000001E-002 , error =    1.6029303120390637E-006
 Observed order =  3.02
```

The final time must be an exact multiple of the time step. Each problem prints its own options with
`./build/tester/foodie_tester <problem> --help`.
