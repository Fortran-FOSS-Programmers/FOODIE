---
title: Multistep schemes
---

# Multistep schemes

The multistep classes compute the solution at the new step from the solutions at $N_s$ previous steps. They extend
`integrator_multistep_object` (the multistep Runge-Kutta SSP class extends `integrator_multistage_multistep_object`) and
store the previous steps in their members `previous(1:steps)`, with times `t(1:steps)` and time steps `Dt(1:steps)`. The
first $N_s$ steps must be provided by the caller (start-up), then `integrate(U, Dt, t)` shifts them by itself if
`autoupdate` is `.true.` (the default); `t` is the time of the last stored step. See the
[quick start](./quickstart#integrate-with-a-multistep-scheme) for a complete start-up loop.

Below, $U^{n+s}$ denotes the solution at the step $n+s$ and $R^{n+s} = R(t^{n+s}, U^{n+s})$. The value of $\Delta t$ is
always provided by the caller; all the schemes but the variable step size LMM-SSP require a constant $\Delta t$.

## Adams-Bashforth

Class `adams_bashforth` (`integrator_adams_bashforth`), explicit, schemes `adams_bashforth_1` ... `adams_bashforth_16`,
of order $N_s$ equal to the number of steps [7, 12]:

$$ U^{n+N_s} = U^{n+N_s-1} + \Delta t \sum_{s=1}^{N_s} b_s R^{n+s-1} $$

`adams_bashforth_1` is the forward Euler scheme. The first coefficients are:

| Scheme | $b_1$ | $b_2$ | $b_3$ | $b_4$ |
|---|---|---|---|---|
| `adams_bashforth_1` | $1$ | | | |
| `adams_bashforth_2` | $-\frac{1}{2}$ | $\frac{3}{2}$ | | |
| `adams_bashforth_3` | $\frac{5}{12}$ | $-\frac{16}{12}$ | $\frac{23}{12}$ | |
| `adams_bashforth_4` | $-\frac{9}{24}$ | $\frac{37}{24}$ | $-\frac{59}{24}$ | $\frac{55}{24}$ |

## Adams-Moulton

Class `adams_moulton` (`integrator_adams_moulton`), implicit, schemes `adams_moulton_0` ... `adams_moulton_15`, of
order $N_s + 1$ [7, 12]:

$$ U^{n+N_s} = U^{n+N_s-1} + \Delta t \left[ \sum_{s=0}^{N_s-1} b_s R^{n+s} + b_{N_s} R\left(t^{n+N_s}, U^{n+N_s}\right) \right] $$

The implicit equation is solved by fixed point iterations, their number given by the `iterations` argument (default 1):
each iteration evaluates the residual at the last iterate of $U^{n+N_s}$, starting from the value of `U` passed to
`integrate`. `adams_moulton_0` is the **implicit backward Euler** scheme, 1st order, $U^{n+1} = U^n + \Delta t\,
R(t^{n+1}, U^{n+1})$, with $U^n$ stored in the single register of the previous steps; at least one fixed point iteration
is always performed.

| Scheme | $b_0$ | $b_1$ | $b_2$ | $b_3$ |
|---|---|---|---|---|
| `adams_moulton_0` | $1$ | | | |
| `adams_moulton_1` | $\frac{1}{2}$ | $\frac{1}{2}$ | | |
| `adams_moulton_2` | $-\frac{1}{12}$ | $\frac{8}{12}$ | $\frac{5}{12}$ | |
| `adams_moulton_3` | $\frac{1}{24}$ | $-\frac{5}{24}$ | $\frac{19}{24}$ | $\frac{9}{24}$ |

Fixed point iterations converge only if $\Delta t$ is small with respect to the Lipschitz constant of $R$: on stiff
problems they impose a step restriction comparable to the one of an explicit scheme.

## Adams-Bashforth-Moulton

Class `adams_bashforth_moulton` (`integrator_adams_bashforth_moulton`), predictor-corrector, schemes
`adams_bashforth_moulton_1` ... `adams_bashforth_moulton_16`. The scheme of number $N$ predicts with the Adams-Bashforth
scheme of $N$ steps and corrects with the Adams-Moulton scheme of $N-1$ steps, both of order $N$:

$$
\begin{aligned}
U^{n+N,p} &= U^{n+N-1} + \Delta t \sum_{s=1}^{N} b_s^p R^{n+s-1} \\
U^{n+N} &= U^{n+N-1} + \Delta t \left[ \sum_{s=0}^{N-2} b_s^c R^{n+s+1} + b_{N-1}^c R\left(t^{n+N}, U^{n+N,p}\right) \right]
\end{aligned}
$$

`adams_bashforth_moulton_1` is AB(1)-AM(0), forward Euler predictor and backward Euler corrector, 1st order. The
`iterations` argument (default 1) sets the corrector iterations.

## Backward Differentiation Formula

Class `back_df` (`integrator_back_df`), implicit, schemes `back_df_1` ... `back_df_6`, of order $N_s$:

$$ U^{n+N_s} + \sum_{s=1}^{N_s} \alpha_s U^{n+N_s-s} = \Delta t\, \beta R\left(t^{n+N_s}, U^{n+N_s}\right) $$

solved by fixed point iterations (`iterations`, default 1). `back_df_1` is the backward Euler scheme.

| Steps | $\beta$ | $\alpha_1$ | $\alpha_2$ | $\alpha_3$ | $\alpha_4$ | $\alpha_5$ | $\alpha_6$ |
|---|---|---|---|---|---|---|---|
| 1 | $1$ | $-1$ | | | | | |
| 2 | $\frac{2}{3}$ | $-\frac{4}{3}$ | $\frac{1}{3}$ | | | | |
| 3 | $\frac{6}{11}$ | $-\frac{18}{11}$ | $\frac{9}{11}$ | $-\frac{2}{11}$ | | | |
| 4 | $\frac{12}{25}$ | $-\frac{48}{25}$ | $\frac{36}{25}$ | $-\frac{16}{25}$ | $\frac{3}{25}$ | | |
| 5 | $\frac{60}{137}$ | $-\frac{300}{137}$ | $\frac{300}{137}$ | $-\frac{200}{137}$ | $\frac{75}{137}$ | $-\frac{12}{137}$ | |
| 6 | $\frac{60}{147}$ | $-\frac{360}{147}$ | $\frac{450}{147}$ | $-\frac{400}{147}$ | $\frac{225}{147}$ | $-\frac{72}{147}$ | $\frac{10}{147}$ |

The BDF schemes are designed for stiff problems, but their stability advantage requires a Newton-like solver of the
implicit equation: with fixed point iterations the step is limited as for the Adams-Moulton schemes.

## Leapfrog

Class `leapfrog` (`integrator_leapfrog`), explicit, 2 steps [4]:

$$ U^{n+2} = U^n + 2\Delta t\, R^{n+1} $$

The scheme `leapfrog` is 2nd order accurate, but it is unstable on dissipative problems (its spurious computational
mode is amplified): use it for oscillatory, non dissipative problems.

The scheme `leapfrog_raw` applies the Robert-Asselin-Williams (RAW) filter [5, 6, 20, 21] after each step:

$$
\begin{aligned}
\Delta &= \frac{\nu}{2}\left(U^n - 2U^{n+1} + U^{n+2}\right) \\
U^{n+1} &\leftarrow U^{n+1} + \alpha \Delta \\
U^{n+2} &\leftarrow U^{n+2} + (\alpha - 1) \Delta
\end{aligned}
$$

with $\nu \in (0, 1]$ and $\alpha \in (0.5, 1]$, defaults $\nu = 0.01$ and $\alpha = 0.53$, set by the `nu` and `alpha`
arguments. For $\alpha = 1$ the filter reverts to the classic Robert-Asselin filter.

::: warning Order of accuracy
With a fixed filter coefficient $\nu$, `leapfrog_raw` is formally **1st order** accurate, as the tests suite measures.
The "3rd order accuracy" of Williams (2011) refers to the accuracy of the amplitude of the oscillations, not to the order
of the global error.
:::

## Linear multistep SSP

Class `lmm_ssp` (`integrator_lmm_ssp`), explicit, Strong Stability Preserving [16]:

$$ U^{n+N_s} = \sum_{s=1}^{N_s} \left[ a_s U^{n+s-1} + \Delta t\, b_s R^{n+s-1} \right] $$

| Scheme | Order | Non-zero coefficients |
|---|---|---|
| `lmm_ssp_steps_3_order_2` | 2 | $a_1 = \frac{1}{4}$, $a_3 = \frac{3}{4}$, $b_3 = \frac{3}{2}$ |
| `lmm_ssp_steps_4_order_3` | 3 | $a_1 = \frac{11}{27}$, $a_4 = \frac{16}{27}$, $b_1 = \frac{12}{27}$, $b_4 = \frac{16}{9}$ |
| `lmm_ssp_steps_5_order_3` | 3 | $a_1 = \frac{7}{32}$, $a_5 = \frac{25}{32}$, $b_1 = \frac{5}{16}$, $b_5 = \frac{25}{16}$ |

## Linear multistep SSP with variable step size

Class `lmm_ssp_vss` (`integrator_lmm_ssp_vss`), explicit, SSP, the time step can change from step to step [17]. With

$$ \omega_i = \frac{\Delta t^{n+i}}{\Delta t^{n+N_s}}, \qquad \Omega_s = \sum_{i=1}^{s} \omega_i, \quad 1 \le s \le N_s $$

the 2nd order formula is

$$ U^{n+N_s} = \frac{1}{\Omega_{N_s-1}^2} U^n + \frac{\Omega_{N_s-1}^2 - 1}{\Omega_{N_s-1}^2} U^{n+N_s-1}
   + \frac{\Omega_{N_s-1} + 1}{\Omega_{N_s-1}} \Delta t^{n+N_s} R\left(U^{n+N_s-1}\right) $$

and the 3rd order one is

$$
\begin{aligned}
U^{n+N_s} = {} & \frac{3 \Omega_{N_s-1} + 2}{\Omega_{N_s-1}^3} U^n
 + \frac{(\Omega_{N_s-1} + 1)^2 (\Omega_{N_s-1} - 2)}{\Omega_{N_s-1}^3} U^{n+N_s-1} \\
 & + \frac{\Omega_{N_s-1} + 1}{\Omega_{N_s-1}^2} \Delta t^{n+N_s} R(U^n)
 + \frac{(\Omega_{N_s-1} + 1)^2}{\Omega_{N_s-1}^2} \Delta t^{n+N_s} R\left(U^{n+N_s-1}\right)
\end{aligned}
$$

| Scheme | Order |
|---|---|
| `lmm_ssp_vss_steps_2_order_2` | 2 |
| `lmm_ssp_vss_steps_3_order_2` | 2 |
| `lmm_ssp_vss_steps_3_order_3` | 3 |
| `lmm_ssp_vss_steps_4_order_3` | 3 |
| `lmm_ssp_vss_steps_5_order_3` | 3 |

## Multistep Runge-Kutta SSP

Class `ms_runge_kutta_ssp` (`integrator_ms_runge_kutta_ssp`), explicit, multistep-multistage, SSP [18]. With $k$ steps
and $s$ stages:

$$
\begin{aligned}
y_1^n &= u^n \\
y_i^n &= \sum_{l=1}^{k} d_{il}\, u^{n-k+l} + \Delta t \sum_{l=1}^{k-1} \hat{a}_{il}\, R(u^{n-k+l})
        + \Delta t \sum_{j=1}^{i-1} a_{ij}\, R(y_j^n), \quad 2 \le i \le s \\
u^{n+1} &= \sum_{l=1}^{k} \theta_l\, u^{n-k+l} + \Delta t \sum_{l=1}^{k-1} \hat{b}_l\, R(u^{n-k+l})
          + \Delta t \sum_{j=1}^{s} b_j\, R(y_j^n)
\end{aligned}
$$

| Scheme | Steps | Stages | Order |
|---|---|---|---|
| `ms_runge_kutta_ssp_steps_2_stages_2_order_3` | 2 | 2 | 3 |
| `ms_runge_kutta_ssp_steps_3_stages_2_order_3` | 3 | 2 | 3 |
| `ms_runge_kutta_ssp_steps_4_stages_5_order_8` | 4 | 5 | 8 |

The coefficients are in the source,
[`foodie_integrator_ms_runge_kutta_ssp.F90`](https://github.com/Fortran-FOSS-Programmers/FOODIE/blob/master/src/lib/foodie_integrator_ms_runge_kutta_ssp.F90).
The 8th order of the last scheme is verified by the polynomial exactness test and observed on the Riccati problem, but
it is beyond the orders the convergence test can assert in double precision (see [Testing](./testing)).

The references are numbered as in the [bibliography of the catalogue](./schemes#bibliography).
