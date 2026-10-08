---
title: Catalogue of the schemes
---

# Catalogue of the schemes

FOODIE provides 86 schemes grouped in 13 classes. A scheme is selected by its name, passed to
`foodie_integrator_factory` or to the `initialize` method of its class (see the [quick start](./quickstart)). The names
are the authoritative list returned by `foodie_integrator_schemes()`; the CLI tester prints them for an unknown scheme:

```bash
fobis build --mode tester-gnu
./build/tester/foodie_tester test --scheme unknown -Dt 0.1 lcce
```

The formulas of each class are on the pages of the [Runge-Kutta](./schemes-runge-kutta) and of the
[multistep](./schemes-multistep) families.

## Summary

| Class (type) | Family | Schemes | Order | Ref. |
|---|---|---|---|---|
| `euler_explicit` (`integrator_euler_explicit`) | explicit, 1 stage | `euler_explicit` | 1 | |
| `runge_kutta_ssp` (`integrator_runge_kutta_ssp`) | explicit, multistage | `runge_kutta_ssp_stages_1_order_1`, `_2_order_2`, `_3_order_3`, `_5_order_4` | 1, 2, 3, 4 | [1, 19] |
| `runge_kutta_ls` (`integrator_runge_kutta_ls`) | explicit, multistage, low storage (2N) | `runge_kutta_ls_stages_1_order_1`, `_5_order_4`, `_6_order_4`, `_7_order_4`, `_12_order_4`, `_13_order_4`, `_14_order_4` | 1, then 4 | [2, 3, 9, 10] |
| `runge_kutta_emd` (`integrator_runge_kutta_emd`) | explicit, multistage, embedded pairs (adaptive step) | `runge_kutta_emd_stages_2_order_2`, `_6_order_5`, `_7_order_4`, `_9_order_6`, `_17_order_10` | 2, 5, **5**, 6, 10 | [11, 13, 14, 15] |
| `runge_kutta_lssp` (`integrator_runge_kutta_lssp`) | explicit, multistage, linear SSP, any number of stages $s$ | `runge_kutta_lssp_stages_s_order_s_1`, `runge_kutta_lssp_stages_s_order_s` | $s-1$, $s$: **linear autonomous problems only** | [16] |
| `adams_bashforth` (`integrator_adams_bashforth`) | explicit, multistep | `adams_bashforth_1` ... `adams_bashforth_16` | $N$ = steps | [7, 12] |
| `adams_moulton` (`integrator_adams_moulton`) | implicit, multistep, fixed point iterations | `adams_moulton_0` ... `adams_moulton_15` | $N+1$ | [7, 12] |
| `adams_bashforth_moulton` (`integrator_adams_bashforth_moulton`) | predictor-corrector, multistep | `adams_bashforth_moulton_1` ... `adams_bashforth_moulton_16` | $N$ | [7, 12] |
| `back_df` (`integrator_back_df`) | implicit, multistep, fixed point iterations | `back_df_1` ... `back_df_6` | $N$ | |
| `leapfrog` (`integrator_leapfrog`) | explicit, 2 steps | `leapfrog`, `leapfrog_raw` | 2, **1** | [4, 5, 6, 20, 21] |
| `lmm_ssp` (`integrator_lmm_ssp`) | explicit, multistep, SSP | `lmm_ssp_steps_3_order_2`, `lmm_ssp_steps_4_order_3`, `lmm_ssp_steps_5_order_3` | 2, 3, 3 | [16] |
| `lmm_ssp_vss` (`integrator_lmm_ssp_vss`) | explicit, multistep, SSP, variable step size | `lmm_ssp_vss_steps_2_order_2`, `_3_order_2`, `_3_order_3`, `_4_order_3`, `_5_order_3` | 2, 2, 3, 3, 3 | [17] |
| `ms_runge_kutta_ssp` (`integrator_ms_runge_kutta_ssp`) | explicit, multistep-multistage, SSP | `ms_runge_kutta_ssp_steps_2_stages_2_order_3`, `_steps_3_stages_2_order_3`, `_steps_4_stages_5_order_8` | 3, 3, 8 | [18] |

The order is the one measured or verified by the [tests suite](./testing); it coincides with the one in the scheme name
but for the cases in bold, explained below.

## Notes on the orders

These facts come from the numerical verification of the schemes (see [Testing](./testing)):

- **`runge_kutta_emd_stages_7_order_4` (Dormand-Prince) is 5th order accurate.** It advances the solution with the 5th
  order formula of the pair (local extrapolation), the 4th order formula being used only for the estimation of the
  local error. The name refers to the lower order of the pair. The other embedded pairs advance with the higher order
  formula too: `runge_kutta_emd_stages_17_order_10` (Feagin) advances with the 10th order one and estimates the error
  with an 8th order one.
- **`leapfrog_raw` is formally 1st order accurate.** It is the leapfrog scheme with the Robert-Asselin-Williams (RAW)
  filter, with a fixed filter coefficient $\nu$. The "3rd order accuracy" of Williams (2011) refers to the accuracy of the
  amplitude of the oscillations, not to the order of the global error. The unfiltered `leapfrog` is 2nd order, but it is
  unstable on dissipative problems: use it for oscillatory (non dissipative) ones.
- **`runge_kutta_lssp_*` (linear SSP) schemes are $s$ (`order_s`) and $s-1$ (`order_s_1`) order accurate only for linear
  autonomous problems**, $U_t = LU$ with $L$ constant. For nonlinear or non-autonomous problems they are at most 2nd order
  accurate. The number of stages $s$ is selected with the `stages` argument of the factory (or of `initialize`); the
  defaults are 2 stages for `order_s_1` and 1 stage for `order_s`.
- **`adams_moulton_0` is the implicit backward Euler scheme**, 1st order: $U^{n+1} = U^n + \Delta t\, R(t^{n+1},
  U^{n+1})$, solved by fixed point iterations (at least one is always performed).
- **`runge_kutta_emd_stages_2_order_2` (Heun-Euler) is supported**: a 2nd order pair with a 1st order (forward Euler)
  error estimator.
- Orders higher than 6 cannot be measured in double precision (the round-off is reached before the asymptotic range):
  they are verified by the polynomial exactness test, complete for the linear multistep schemes (Adams, BDF, LMM-SSP),
  and by the audit of the tableaus for the Runge-Kutta ones. For `ms_runge_kutta_ssp_steps_4_stages_5_order_8` the 8th
  order is verified by the polynomial exactness test (a necessary condition) and observed on the Riccati problem
  (observed order 8.4), but not asserted by the convergence test.

## Options of the classes

| Option (factory/`initialize` argument) | Classes | Default | Meaning |
|---|---|---|---|
| `U` | all | | integrand prototype: allocates the internal registers with your type |
| `stages` | `runge_kutta_lssp` | 2 (`order_s_1`), 1 (`order_s`) | number of stages $s$ |
| `tolerance` | `runge_kutta_emd` | 0.01 | tolerance on the local truncation error estimate |
| `nu`, `alpha` | `leapfrog` (`leapfrog_raw` only) | 0.01, 0.53 | RAW filter coefficients, $\nu \in (0,1]$, $\alpha \in (0.5,1]$; $\alpha=1$ gives the Robert-Asselin filter |
| `iterations` | `adams_moulton`, `adams_bashforth_moulton`, `back_df` | 1 | fixed point iterations of the implicit step (`ms_runge_kutta_ssp` accepts it, but its schemes are explicit and do not use it) |
| `autoupdate` | multistep classes | `.true.` | shift the stored previous steps after each step |

## Adaptive time step

The embedded Runge-Kutta schemes compute two solutions of different order with the same stages, and estimate the local
error with the `.lterror.` operator of the integrand. If the error exceeds `tolerance`, the step is repeated with

$$ \Delta t_{new} = 0.9\, \Delta t_{old} \left( \frac{tolerance}{error} \right)^{\frac{1}{p+1}} $$

where $p$ is the order of the error estimator (the lower order formula of the pair), until it is accepted; the step actually used is returned in the optional `new_Dt` argument of `integrate`. The step is
only reduced, never enlarged: choosing the next step is left to the caller.

## Bibliography

[1] *High Order Strong Stability Preserving Time Discretizations*, S. Gottlieb, D. I. Ketcheson, C. W. Shu, Journal of
Scientific Computing, vol. 38, N. 3, 2009, pp. 251--289.

[2] *Low-Storage Runge-Kutta Schemes*, J. H. Williamson, Journal of Computational Physics, vol. 35, 1980, pp. 48--56.

[3] *Fourth-Order 2N-Storage Runge-Kutta Schemes*, M. H. Carpenter, C. A. Kennedy, NASA Technical Memorandum 109112,
June 1994.

[4] *Numerical methods used in atmospheric models*, F. Mesinger, A. Arakawa, Global Atmospheric Research Programme
(GARP), [Technical Report](http://twister.ou.edu/CFD2003/Mesinger_ArakawaGARP.pdf), 1976.

[5] *A Proposed Modification to the Robert-Asselin Time Filter*, P. D. Williams, Monthly Weather Review, vol. 137,
pp. 2538--2546, 2009, doi:10.1175/2009MWR2724.1.

[6] *The RAW filter: An improvement to the Robert-Asselin filter in semi-implicit integrations*, P. D. Williams, Monthly
Weather Review, vol. 139(6), pp. 1996--2007, June 2011.

[7] *Linear multistep method*, [Wikipedia article](https://en.wikipedia.org/wiki/Linear_multistep_method).

[8] *Scientific Software Design: The Object-Oriented Way*, D. Rouson, J. Xia, X. Xu, Cambridge University Press, 2011,
ISBN 9780521888134.

[9] *High-accuracy large-step explicit Runge-Kutta (HALE-RK) schemes for computational aeroacoustics*, V. Allampalli,
R. Hixon, M. Nallasamy, S. D. Sawyer, Journal of Computational Physics, vol. 228, 2009, pp. 3837--3850.

[10] *Efficient low-storage Runge-Kutta schemes with optimized stability regions*, J. Niegemann, R. Diehl, K. Busch,
Journal of Computational Physics, vol. 231, 2012, pp. 364--372.

[11] *A family of embedded Runge-Kutta formulae*, J. R. Dormand, P. J. Prince, Journal of Computational and Applied
Mathematics, vol. 6 (1), 1980, pp. 19--26, doi:10.1016/0771-050X(80)90013-3.

[12] *Cowell Type Numerical Integration As Applied to Satellite Orbit Computation*, J. L. Maury Jr., G. P. Segal,
X-553-69-46, April 1969, [NASA-TM-X-63542](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19690017325.pdf).

[13] *A variable order Runge-Kutta method for initial value problems with rapidly varying right-hand sides*, J. R. Cash,
A. H. Karp, ACM Transactions on Mathematical Software, vol. 16, 1990, pp. 201--222, doi:10.1145/79505.79507.

[14] *A New Embedded Pair of Runge-Kutta Formulas of orders 5 and 6*, M. Calvo, J. I. Montijano, L. Randez, Computers &
Mathematics with Applications, vol. 20, issue 1, 1990, pp. 15--24, doi:10.1016/0898-1221(90)90064-Q.

[15] *A tenth-order Runge-Kutta method with error estimate*, T. Feagin, Proceedings of the IAENG Conference on
Scientific Computing, 2007.

[16] *Strong Stability Preserving Runge-Kutta and Multistep Time Discretizations*, S. Gottlieb, D. Ketcheson, C. W. Shu,
World Scientific Publishing, 2011, doi:10.1142/7498.

[17] *Strong Stability Preserving Explicit Linear Multistep Methods with Variable Step Size*, Y. Hadjimichael,
D. Ketcheson, L. Lóczi, A. Németh, SIAM Journal on Numerical Analysis, vol. 54 (5), 2016, pp. 2799--2832,
doi:10.1137/15M101717X.

[18] *Explicit Strong Stability Preserving Multistep Runge-Kutta Methods*, C. Bresten, S. Gottlieb, Z. Grant, D. Higgs,
D. Ketcheson, A. Németh, Mathematics of Computation, 2016, doi:10.1090/mcom/3115.

[19] *Coefficients for the study of Runge-Kutta integration processes*, J. C. Butcher, Journal of the Australian
Mathematical Society, vol. 3, 1963, pp. 185--201.

[20] *The integration of a low order spectral form of the primitive meteorological equations*, A. J. Robert, Journal of
the Meteorological Society of Japan, vol. 44, 1966, pp. 237--245.

[21] *Frequency filter for time integrations*, R. Asselin, Monthly Weather Review, vol. 100, 1972, pp. 487--490.
