---
title: 1D Euler equations
---

# Example: 1D Euler equations

::: info Migrated from the wiki
This page replaces the *1D Euler equations* page of the old GitHub wiki, whose code was written for a former FOODIE API
(`type, extends(integrand)`, `update_previous_steps`, `adams_bashforth_integrator%init(steps=3)`) that no longer exists.
The current application is `src/app/euler_1D`: the program `foodie_euler_1D.f90` and the integrand
`foodie_integrand_euler_1D.f90`.
:::

## Mathematical model

The 1D Euler equations are a nonlinear hyperbolic system of conservation laws describing the dynamics of a compressible
gas when body forces, viscous stresses and heat fluxes are negligible:

$$
U_t = R(U) \Leftrightarrow U_t + F(U)_x = 0, \qquad
U = \begin{bmatrix} \rho \\ \rho u \\ \rho E \end{bmatrix}, \qquad
F(U) = \begin{bmatrix} \rho u \\ \rho u^2 + p \\ \rho u H \end{bmatrix}
$$

where $\rho$ is the density, $u$ the velocity, $p$ the pressure, $E$ the total specific energy and $H$ the total
specific enthalpy. The gas is ideal (thermally and calorically perfect):

$$
R = c_p - c_v, \quad \gamma = \frac{c_p}{c_v}, \quad e = c_v T, \quad h = c_p T, \quad
T = \frac{p}{\rho R}, \quad a = \sqrt{\frac{\gamma p}{\rho}}
$$

with $\rho E = \rho e + \frac{1}{2}\rho u^2$ and $\rho H = \rho h + \frac{1}{2}\rho u^2$.

### Multi-fluid extension

A mixture of $N_s$ gases with different properties is modelled by the standard thermodynamic model, replacing the
density with the partial densities $\rho_s$ of the species:

$$
U = \begin{bmatrix} \rho_s \\ \rho u \\ \rho E \end{bmatrix}, \quad
F(U) = \begin{bmatrix} \rho_s u \\ \rho u^2 + p \\ \rho u H \end{bmatrix}, \quad s = 1, \dots, N_s, \qquad
\rho = \sum_{s=1}^{N_s} \rho_s, \quad
c_p = \sum_{s=1}^{N_s} \frac{\rho_s}{\rho} c_{p,s}, \quad
c_v = \sum_{s=1}^{N_s} \frac{\rho_s}{\rho} c_{v,s}
$$

The system admits discontinuous solutions (shocks, contact discontinuities): it is a demanding benchmark for the
coupling of a time integrator with a non-oscillatory space discretization, while being simple enough to be controlled.

## Space discretization

The method of lines turns the PDE system into an ODE system for the cell averages. The integrand uses a finite volume,
fully conservative, Godunov-like scheme on $N_i$ uniform cells, with $N_g$ ghost cells on each side for the boundary
conditions (transmissive `TRA` or reflective `REF`):

1. the conservative variables are converted to primitive ones and the boundary conditions imposed on the ghost cells;
2. the primitive variables are reconstructed at the cell interfaces by a WENO interpolator of the
   [WenOOF](https://github.com/Fortran-FOSS-Programmers/WenOOF) library;
3. the numerical flux at each interface is computed by the local Lax-Friedrichs (Rusanov) approximate Riemann solver;
4. the residual of each cell is the flux balance, $R_i = (F_{i-1/2} - F_{i+1/2}) / \Delta x$.

The time step is computed from the CFL condition, $\Delta t = \text{CFL}\, \Delta x / \max_i (|u_i| + a_i)$.

## The FOODIE integrand

`integrand_euler_1D` extends `integrand_object` directly. Its state is the rank 2 array `U(1:Nc, 1:Ni)` of the
conservative variables ($N_c = N_s + 2$) of the physical cells; the integrand also holds the parameters of the grid,
the WENO interpolator and the integrator. The residual function returns the residuals of all the cells flattened in one
array, variables first (excerpt of `dEuler_dt`):

```fortran
function dEuler_dt(self, t) result(dState_dt)
!< Time derivative of Euler field, the residuals function.
class(integrand_euler_1D), intent(in)           :: self         !< Euler field.
real(R_P),                 intent(in), optional :: t            !< Time.
real(R_P), allocatable                          :: dState_dt(:) !< Euler field time derivative.
...
! compute residuals
allocate(dState_dt(1:self%Nc*self%Ni))
j = 0
do i=1, self%Ni
   do c=1, self%Nc
      j = j + 1
      dState_dt(j) = (F(c, i - 1) - F(c, i)) / self%dx
   enddo
enddo
endfunction dEuler_dt
```

The operators flatten the state in the same order and `assign_real` reshapes the array back into `U`:

```fortran
pure subroutine assign_real(lhs, rhs)
!< Assign one real to an advection field.
class(integrand_euler_1D), intent(inout) :: lhs     !< Left hand side.
real(R_P),                 intent(in)    :: rhs(1:) !< Right hand side.
integer(I_P)                             :: c, i, j !< Counter.

j = 0
do i=1, lhs%Ni
   do c=1, lhs%Nc
      j = j + 1
      lhs%U(c,i) = rhs(j)
   enddo
enddo
endsubroutine assign_real
```

The integrator is created by the factory from the scheme selected on the command line, and the time loop advances the
field with the CFL time step:

```fortran
do
   n = n + 1
   dt = self%dt(t=t)
   select type(integrator)
   class is(integrator_multistage_object)
      call integrator%integrate(U=self, dt=dt, t=t)
   endselect
   t = t + dt
   ...
   if (self%is_completed(n=n, t=t)) exit
enddo
```

The WenOOF interpolator is impure, so the app is built with the `_IMPURE_` macro (see
[Pure or impure integrands](./installation#pure-or-impure-integrands)).

## Build and run

The app has its own `fobos` file; build it from the repository root, after `fobis fetch`:

```bash
fobis build -f src/app/euler_1D/fobos --mode gnu
./build/app/euler_1D/foodie_euler_1D --help
```

| Option | Default | Meaning |
|---|---|---|
| `--initial_state` | `sod` | initial state: `sod` (Sod shock tube) or `bryson` |
| `--i-scheme` | `runge_kutta_ssp_stages_5_order_4` | time integrator: `runge_kutta_ssp_stages_1_order_1`, `_3_order_3` or `_5_order_4` |
| `--w-scheme` | `reconstructor-JS` | WENO scheme: `reconstructor-JS`, `-M-JS`, `-M-Z`, `-Z` |
| `--weno-order` | 3 | WENO reconstruction order |
| `--weno-eps` | 1e-6 | WENO epsilon |
| `--cfl` | 0.8 | CFL number |
| `--length` | 1.0 | domain length |
| `--Ni` | 100 | number of cells |
| `--Ns` | 1 | number of species |
| `--BC_L`, `--BC_R` | `TRA` | left and right boundary conditions, `TRA` or `REF` |
| `--Tmax` | 0.2 | final time |
| `--Nmax` | -1 | maximum number of time steps (if positive it replaces `--Tmax`) |

The Sod shock tube with the default settings:

```bash
./build/app/euler_1D/foodie_euler_1D
```

prints the settings and the progress of the integration,

```text
Integrate euler_1D-sod
Nmax              : +3
...
Integrator scheme : runge_kutta_ssp_stages_5_order_4
...
n, t, dt: +1, +0.67618841252329170E-002, +0.67618841252329170E-002
n, t, dt: +2, +0.11376653013898985E-001, +0.46147688886660678E-002
n, t, dt: +3, +0.15713893762532984E-001, +0.43372407486339986E-002
```

(here with `--Nmax 3`), and saves the solution of every step in the CSV files `euler_1D-sod-000000000.csv`,
`euler_1D-sod-000000001.csv`, ..., with columns `x,u,p,r,t` (position, velocity, pressure, density, time). The Bryson
initial state is run in a longer domain, as in the example of the program help:

```bash
./build/app/euler_1D/foodie_euler_1D --initial_state bryson --length 20 --Tmax 6
```

The app supports only the multistage schemes: the time loop selects `integrator_multistage_object`, and a multistep
scheme would need the start-up of its previous steps (see the [quick start](./quickstart#integrate-with-a-multistep-scheme)).

## Bibliography

[1] *Numerical Methods for Fluid Dynamics With Applications to Geophysics*, D. R. Durran, Springer, 2010.
