---
layout: home

hero:
  name: FOODIE
  text: Fortran Object-Oriented Differential-equations Integration Environment
  tagline: "Write your system once as an integrand type, integrate it with 86 explicit, implicit, multistage and multistep schemes selected by name. Pure Fortran 2008, built on the Abstract Data Type pattern of Rouson."
  actions:
    - theme: brand
      text: Quick start
      link: /guide/quickstart
    - theme: alt
      text: Schemes
      link: /guide/schemes
    - theme: alt
      text: Installation
      link: /guide/installation
    - theme: alt
      text: API
      link: /api/
    - theme: alt
      text: View on GitHub
      link: https://github.com/Fortran-FOSS-Programmers/FOODIE

features:
  - icon: 🧮
    title: Formulas that read like formulas
    details: "Every scheme is written against one abstract type, integrand_object: U = U + (U%t(t=t) * Dt) is the forward Euler scheme, not a pseudo-code of it."
    link: /guide/#the-abstract-integrand
    linkText: The design
  - icon: 🧩
    title: Your system, any shape
    details: "Extend integrand_object with your state (a scalar, a Lorenz attractor, a finite volume grid of a multi-fluid gas) and implement its residual function and arithmetic operators: every scheme applies to it."
    link: /guide/quickstart
    linkText: Quick start
  - icon: 🏭
    title: Schemes selected by name
    details: "foodie_integrator_factory returns a ready integrator from a string such as 'runge_kutta_ssp_stages_3_order_3' or 'adams_bashforth_4': switch scheme without recompiling."
    link: /guide/quickstart#select-a-scheme-by-name
    linkText: The factory
  - icon: 📚
    title: Thirteen classes of schemes
    details: "Forward Euler, Runge-Kutta SSP, low-storage, embedded (adaptive), linear SSP; Adams-Bashforth, Adams-Moulton, Adams-Bashforth-Moulton, BDF, leapfrog, LMM-SSP with fixed and variable step size, multistep Runge-Kutta SSP."
    link: /guide/schemes
    linkText: Catalogue
  - icon: ✅
    title: Verified orders
    details: "A tests suite measures the order of every scheme on a nonlinear non-autonomous problem, checks polynomial exactness, and audits every Runge-Kutta tableau against the order conditions of all rooted trees."
    link: /guide/testing
    linkText: Testing
  - icon: ⚡
    title: Fast mode
    details: "Every integrator also provides integrate_fast, built on in-place operators that your integrand can override to avoid the temporaries of the operator syntax."
    link: /guide/quickstart#fast-mode
    linkText: Fast mode
  - icon: 🛠️
    title: Standard Fortran, simple builds
    details: "Fortran 2008, built with FoBiS (dependencies fetched for you) or fpm; static and shared library modes for GNU and Intel compilers."
    link: /guide/installation
    linkText: Installation
  - icon: 🔓
    title: Multi-licensed
    details: "GPL v3 for FOSS projects; BSD 2-Clause, BSD 3-Clause or MIT for closed source and commercial ones: pick the license that fits."
    link: /guide/#copyrights
    linkText: Copyrights
---

## What FOODIE solves

FOODIE integrates in time the initial value problem

$$
\begin{cases}
U_t = R(t, U) \\
U(t_0) = U_0
\end{cases}
$$

where $U_t = \frac{dU}{dt}$, $U$ is the vector of state variables, function of the time-like independent variable
$t$, and $R$ is the (vectorial) residual function, possibly nonlinear in $U$. FOODIE is tailored to the systems arising
from the semi-discretization of PDEs (the method of lines), but it is not limited to them.

## A first integration

A concrete integrand, here the linear equation $U_t = aU + b$ of the tests suite, is integrated by any scheme through the
factory:

```fortran
use foodie, only : foodie_integrator_factory, integrator_object, integrator_multistage_object
...
class(integrator_object), allocatable :: integrator

call foodie_integrator_factory(scheme='runge_kutta_ssp_stages_3_order_3', integrator=integrator, U=integrand)
select type(integrator)
class is(integrator_multistage_object)
   do step=1, 100
      call integrator%integrate(U=integrand, Dt=0.01_R_P, t=time)
      time = time + 0.01_R_P
   enddo
endselect
```

The [quick start](/guide/quickstart) builds the integrand and the loop step by step; the [catalogue](/guide/schemes)
lists every scheme with its order and its references.

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

Contributions are welcome — see the [Contributing](/guide/contributing) page.

## Copyrights

This project is distributed under a multi-licensing system:

- **FOSS projects**: [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html)
- **Closed source / commercial**: [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause), [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause), or [MIT](http://opensource.org/licenses/MIT)

> Anyone interested in using, developing, or contributing to this project is welcome — pick the license that best fits your needs.
