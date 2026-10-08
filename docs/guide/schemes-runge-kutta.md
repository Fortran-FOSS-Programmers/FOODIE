---
title: Runge-Kutta schemes
---

# Runge-Kutta schemes

The multistage classes advance the solution from $t^n$ to $t^{n+1} = t^n + \Delta t$ using only $U^n$. All of them are
explicit and extend `integrator_multistage_object`: their `integrate(U, Dt, t, new_Dt)` takes the time step `Dt`, given
by the caller, and the time `t` $= t^n$. The orders stated here are the verified ones, see the
[catalogue](./schemes#notes-on-the-orders) and [Testing](./testing).

## Forward Euler

Class `euler_explicit`, scheme `euler_explicit`, 1st order:

$$ U^{n+1} = U^n + \Delta t\, R(t^n, U^n) $$

## SSP Runge-Kutta

Class `runge_kutta_ssp` (`integrator_runge_kutta_ssp`). The schemes have the Total Variation Diminishing (TVD) or Strong
Stability Preserving (SSP) property and are written in the Butcher form [19]:

$$ U^{n+1} = U^n + \Delta t \sum_{s=1}^{N_s} \beta^s K^s, \qquad
   K^s = R\left( t^n + \gamma^s \Delta t,\; U^n + \Delta t \sum_{i=1}^{s-1} \alpha^{s,i} K^i \right) $$

where $N_s$ is the number of stages and $\gamma^s = \sum_i \alpha^{s,i}$.

| Scheme | Order | Notes |
|---|---|---|
| `runge_kutta_ssp_stages_1_order_1` | 1 | forward Euler |
| `runge_kutta_ssp_stages_2_order_2` | 2 | optimal SSP(2,2) [1] |
| `runge_kutta_ssp_stages_3_order_3` | 3 | optimal SSP(3,3) [1] |
| `runge_kutta_ssp_stages_5_order_4` | 4 | optimal SSP(5,4) [1] |

The tableaus ($\gamma$ in the first column, $\alpha$ in the matrix, $\beta$ in the last row):

$$
\text{SSP(2,2)}\quad
\begin{array}{c|cc}
0 & & \\
1 & 1 & \\ \hline
  & \frac{1}{2} & \frac{1}{2}
\end{array}
\qquad\qquad
\text{SSP(3,3)}\quad
\begin{array}{c|ccc}
0 & & & \\
1 & 1 & & \\
\frac{1}{2} & \frac{1}{4} & \frac{1}{4} & \\ \hline
  & \frac{1}{6} & \frac{1}{6} & \frac{2}{3}
\end{array}
$$

SSP(5,4):

| $\gamma$ | $\alpha^{s,1}$ | $\alpha^{s,2}$ | $\alpha^{s,3}$ | $\alpha^{s,4}$ |
|---|---|---|---|---|
| 0 | | | | |
| 0.39175222686925 | 0.39175222686925 | | | |
| 0.58607968906690 | 0.21766909635783 | 0.36841059270907 | | |
| 0.47454236316248 | 0.08269208668309 | 0.13995850210743 | 0.25189177437196 | |
| 0.93501063109579 | 0.06796628357405 | 0.11503469845367 | 0.20703489877294 | 0.54497475029514 |

$\beta$ = 0.14681187615788, 0.24848290939132, 0.10425883027948, 0.27443890104848, 0.22600748312284.

These are the rounded values of the module documentation; the source uses more digits, refined to satisfy the order
conditions in double precision (verified by `scripts/audit_rk_tableaus.py`).

## Low storage Runge-Kutta

Class `runge_kutta_ls` (`integrator_runge_kutta_ls`). Following Williamson [2], the schemes use only two registers
$K_1, K_2$ of the size of the state (2N storage):

$$
\begin{aligned}
& K_1 = U^n, \quad K_2 = 0 \\
& \left.\begin{aligned}
  K_2 &= A_s K_2 + \Delta t\, R(t^n + C_s \Delta t, K_1) \\
  K_1 &= K_1 + B_s K_2
  \end{aligned}\right\} \quad s = 1, 2, \dots, N_s \\
& U^{n+1} = K_1
\end{aligned}
$$

with $A_1 = C_1 = 0$.

| Scheme | Order | Reference |
|---|---|---|
| `runge_kutta_ls_stages_1_order_1` | 1 | forward Euler (not a real low storage scheme, provided for completeness) |
| `runge_kutta_ls_stages_5_order_4` | 4 | LSRK(5,4)2N, solution 3 of Carpenter and Kennedy [3] |
| `runge_kutta_ls_stages_6_order_4` | 4 | RK(6,4) of Allampalli et al. [9] |
| `runge_kutta_ls_stages_7_order_4` | 4 | RK(7,4) of Allampalli et al. [9] |
| `runge_kutta_ls_stages_12_order_4` | 4 | RK(12,4) of Niegemann et al. [10] |
| `runge_kutta_ls_stages_13_order_4` | 4 | RK(13,4) of Niegemann et al. [10] |
| `runge_kutta_ls_stages_14_order_4` | 4 | RK(14,4) of Niegemann et al. [10] |

The coefficients of the 5 stages scheme are rational:

| Stage | $A$ | $B$ | $C$ |
|---|---|---|---|
| 1 | 0 | 1432997174477/9575080441755 | 0 |
| 2 | -567301805773/1357537059087 | 5161836677717/13612068292357 | 1432997174477/9575080441755 |
| 3 | -2404267990393/2016746695238 | 1720146321549/2090206949498 | 2526269341429/6820363962896 |
| 4 | -3550918686646/2091501179385 | 3134564353537/4481467310338 | 2006345519317/3224310063776 |
| 5 | -1275806237668/842570457699 | 2277821191437/14882151754819 | 2802321613138/2924317926251 |

The 6 and 7 stages coefficients published with 12 decimals violate the 4th order conditions by about $10^{-13}$: the
source uses coefficients refined to satisfy them exactly, which reproduce the published ones once rounded. The
coefficients of all the schemes are in the documentation of
[`foodie_integrator_runge_kutta_low_storage.F90`](https://github.com/Fortran-FOSS-Programmers/FOODIE/blob/master/src/lib/foodie_integrator_runge_kutta_low_storage.F90)
and in the [API](/api/).

## Embedded Runge-Kutta (adaptive)

Class `runge_kutta_emd` (`integrator_runge_kutta_emd`). An embedded pair computes, with the same stages $K^s$, two
solutions of different order:

$$ U_{high}^{n+1} = U^n + \Delta t \sum_{s=1}^{N_s} \beta_{high}^s K^s, \qquad
   U_{low}^{n+1}  = U^n + \Delta t \sum_{s=1}^{N_s} \beta_{low}^s K^s $$

The difference between the two, measured by the `.lterror.` operator of the integrand, estimates the local error and
drives the step size control (see [Adaptive time step](./schemes#adaptive-time-step)). In all the FOODIE pairs the
solution is advanced with the **higher** order formula (local extrapolation):

| Scheme | Advances with | Error estimator | Reference |
|---|---|---|---|
| `runge_kutta_emd_stages_2_order_2` | 2nd order | 1st order | Heun-Euler |
| `runge_kutta_emd_stages_6_order_5` | 5th order | 4th order | Cash-Karp [13] |
| `runge_kutta_emd_stages_7_order_4` | **5th order** | 4th order | Dormand-Prince [11] |
| `runge_kutta_emd_stages_9_order_6` | 6th order | 5th order | Calvo et al. [14] |
| `runge_kutta_emd_stages_17_order_10` | 10th order | 8th order | Feagin [15] |

The names of the schemes give the number of stages and one order of the pair: for Dormand-Prince it is the lower one,
the scheme being 5th order accurate.

Heun-Euler:

$$
\begin{array}{c|cc}
0 & & \\
1 & 1 & \\ \hline
\beta_{high} & \frac{1}{2} & \frac{1}{2} \\
\beta_{low} & 1 & 0
\end{array}
$$

Cash-Karp:

$$
\begin{array}{c|cccccc}
0 & & & & & & \\
\frac{1}{5} & \frac{1}{5} & & & & & \\
\frac{3}{10} & \frac{3}{40} & \frac{9}{40} & & & & \\
\frac{3}{5} & \frac{3}{10} & -\frac{9}{10} & \frac{6}{5} & & & \\
1 & -\frac{11}{54} & \frac{5}{2} & -\frac{70}{27} & \frac{35}{27} & & \\
\frac{7}{8} & \frac{1631}{55296} & \frac{175}{512} & \frac{575}{13824} & \frac{44275}{110592} & \frac{253}{4096} & \\ \hline
\beta_{high} & \frac{37}{378} & 0 & \frac{250}{621} & \frac{125}{594} & 0 & \frac{512}{1771} \\
\beta_{low} & \frac{2825}{27648} & 0 & \frac{18575}{48384} & \frac{13525}{55296} & \frac{277}{14336} & \frac{1}{4}
\end{array}
$$

Dormand-Prince:

$$
\begin{array}{c|ccccccc}
0 & & & & & & & \\
\frac{1}{5} & \frac{1}{5} & & & & & & \\
\frac{3}{10} & \frac{3}{40} & \frac{9}{40} & & & & & \\
\frac{4}{5} & \frac{44}{45} & -\frac{56}{15} & \frac{32}{9} & & & & \\
\frac{8}{9} & \frac{19372}{6561} & -\frac{25360}{2187} & \frac{64448}{6561} & -\frac{212}{729} & & & \\
1 & \frac{9017}{3168} & -\frac{355}{33} & \frac{46732}{5247} & \frac{49}{176} & -\frac{5103}{18656} & & \\
1 & \frac{35}{384} & 0 & \frac{500}{1113} & \frac{125}{192} & -\frac{2187}{6784} & \frac{11}{84} & \\ \hline
\beta_{high} & \frac{35}{384} & 0 & \frac{500}{1113} & \frac{125}{192} & -\frac{2187}{6784} & \frac{11}{84} & 0 \\
\beta_{low} & \frac{5179}{57600} & 0 & \frac{7571}{16695} & \frac{393}{640} & -\frac{92097}{339200} & \frac{187}{2100} & \frac{1}{40}
\end{array}
$$

The Calvo (9 stages) and Feagin (17 stages) tableaus are too large to be reported here: they are in the documentation
of
[`foodie_integrator_runge_kutta_embedded.F90`](https://github.com/Fortran-FOSS-Programmers/FOODIE/blob/master/src/lib/foodie_integrator_runge_kutta_embedded.F90)
and in the [API](/api/).

## Linear SSP Runge-Kutta

Class `runge_kutta_lssp` (`integrator_runge_kutta_lssp`), from Gottlieb, Ketcheson and Shu [16]: SSP schemes with any
number of stages $s$, selected by the `stages` argument. With $h = \frac{1}{2}$ for `runge_kutta_lssp_stages_s_order_s_1`
($s \ge 2$, default 2) and $h = 1$ for `runge_kutta_lssp_stages_s_order_s` ($s \ge 1$, default 1):

$$
\begin{aligned}
U^{(1)} &= U^n \\
U^{(k)} &= U^{(k-1)} + h \Delta t\, R\left(t^n + (k-2) h \Delta t,\; U^{(k-1)}\right), \quad k = 2, \dots, s \\
U^{(s)} &\leftarrow U^{(s)} + h \Delta t\, R\left(t^n + (s-1) h \Delta t,\; U^{(s)}\right) \\
U^{n+1} &= \sum_{k=1}^{s} \alpha_k U^{(k)}
\end{aligned}
$$

The coefficients $\alpha_k$ are computed by recursion at initialization:

- `order_s_1`: $\alpha_1 = 0$, $\alpha_2 = 1$; for $i = 3, \dots, s$: $\alpha_i = \frac{2}{i}\alpha_{i-1}$, then
  $\alpha_j = \frac{2}{j-1}\alpha_{j-1}$ for $j = i-1, \dots, 2$, then $\alpha_1 = 1 - \sum_{j=2}^{i}\alpha_j$;
- `order_s`: $\alpha_1 = 1$; for $i = 2, \dots, s$: $\alpha_i = \frac{1}{i!}$, then $\alpha_j = \frac{1}{j-1}\alpha_{j-1}$
  for $j = i-1, \dots, 2$, then $\alpha_1 = 1 - \sum_{j=2}^{i}\alpha_j$.

::: warning Order of accuracy
The schemes are $s$ and $s-1$ order accurate only for **linear autonomous** problems, $U_t = LU$ with $L$ constant. For
nonlinear or non-autonomous problems they are at most 2nd order accurate. The tests suite verifies their order on the
linear constant coefficients equation.
:::

The references are numbered as in the [bibliography of the catalogue](./schemes#bibliography).
