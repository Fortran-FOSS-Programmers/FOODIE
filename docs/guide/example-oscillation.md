---
title: Oscillation equations
---

# Example: oscillation equations

::: info Migrated from the wiki
This page replaces the *Oscillation equations* page of the old GitHub wiki. The wiki described an integrand written for
a former FOODIE API (`type, extends(integrand)`, with an `update_previous_steps` method and integrator-specific
`init` calls) that no longer exists. The current integrand is `integrand_oscillation` of
`src/tests/tester/foodie_test_integrand_oscillation.f90`, integrated by the CLI tester.
:::

## Mathematical model

The inertial oscillation equations are a system of two ODEs [1]:

$$
U_t = R(U), \qquad
U = \begin{bmatrix} v_1 \\ v_2 \end{bmatrix}, \qquad
R(U) = \begin{bmatrix} -f v_2 \\ f v_1 \end{bmatrix}
$$

with constant frequency $f$, by default $f = 10^{-4}$. The exact solution is a rotation in the $v_1$-$v_2$ plane:

$$
v_1(t) = v_1(0) \cos(f t) - v_2(0) \sin(f t), \qquad
v_2(t) = v_1(0) \sin(f t) + v_2(0) \cos(f t)
$$

The path of the solution must be a circle, its amplitude constant: the eigenvalues of the system, $\pm i f$, are purely
imaginary, so the problem tests the behaviour of the schemes on the imaginary axis. A scheme whose stability region does
not contain a segment of the imaginary axis, such as the forward Euler scheme, amplifies the solution at every step
whatever the step size; the leapfrog scheme, unstable on dissipative problems, is neutrally stable here.

## The FOODIE integrand

The state is an array of two reals; the residual function is, verbatim from the source:

```fortran
pure function dU_dt(self, t) result(dState_dt)
!< Time derivative of integrand_oscillation field.
class(integrand_oscillation), intent(in)           :: self         !< Integrand.
real(R_P),                    intent(in), optional :: t            !< Time.
real(R_P), allocatable                             :: dState_dt(:) !< Integrand time derivative.

dState_dt = [-self%f * self%U(2), &
              self%f * self%U(1)]
endfunction dU_dt
```

The operators and the assignments follow the pattern of the [quick start](./quickstart#define-the-integrand), acting on
the two components of the state.

## Integration

Build the CLI tester and integrate the problem with three time steps up to $t = 10^6$ (about 16 periods), from the
default initial state $(v_1, v_2) = (0, 1)$:

```bash
fobis build --mode tester-gnu
./build/tester/foodie_tester test -s runge_kutta_ssp_stages_5_order_4 -Dt 1000 500 250 -ft 1e6 oscillation
```

The tester prints the error of $v_1$ and $v_2$ at the final time and the observed orders. The errors for some schemes:

| Scheme | $\Delta t = 1000$ | $\Delta t = 500$ | $\Delta t = 250$ | Observed order |
|---|---|---|---|---|
| `euler_explicit` | 109.4 | 6.49 | 1.32 | (unstable) |
| `leapfrog` | 1.50e-1 | 3.64e-2 | 9.01e-3 | 2.01 |
| `leapfrog_raw` | 1.52e-1 | 3.70e-2 | 9.27e-3 | 2.00 |
| `adams_bashforth_3` | 2.18e-2 | 2.58e-3 | 3.10e-4 | 3.06 |
| `adams_moulton_2` | 2.25e-3 | 2.72e-4 | 3.35e-5 | 3.02 |
| `runge_kutta_ssp_stages_3_order_3` | 2.38e-3 | 2.81e-4 | 3.41e-5 | 3.05 |
| `runge_kutta_ssp_stages_5_order_4` | 3.19e-5 | 2.04e-6 | 1.29e-7 | 3.99 |

(error of $v_1$; observed order on the two finest steps; default options: 1 fixed point iteration for Adams-Moulton.)

- The forward Euler solution spirals outwards: with $f \Delta t = 0.1$ its amplitude grows by $\sqrt{1 + (f\Delta
  t)^2}$ per step, the error (109 for an amplitude 1) is not a truncation error but an instability.
- The leapfrog family is neutrally stable on this problem. At these step sizes `leapfrog_raw` shows an observed order of
  about 2: its 1st order error term is proportional to the filter coefficient $\nu = 0.01$, and dominates only for
  smaller steps. Its formal order is 1, see [Notes on the orders](./schemes#notes-on-the-orders).
- The multistep schemes are started up with the exact solution by the tester.

Add `-r` (`--save_results`) to save the solutions as Tecplot ASCII files and plot the path in the $v_1$-$v_2$ plane.

## Bibliography

[1] *Numerical Methods for Fluid Dynamics With Applications to Geophysics*, D. R. Durran, Springer, 2010.
