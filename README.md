<a name="top"></a>

# FOODIE

[![GitHub tag](https://img.shields.io/github/tag/Fortran-FOSS-Programmers/FOODIE.svg)]() [![Join the chat at https://gitter.im/Fortran-FOSS-Programmers/FOODIE](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/Fortran-FOSS-Programmers/FOODIE?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)]()
[![License](https://img.shields.io/badge/license-BSD2-red.svg)]()
[![License](https://img.shields.io/badge/license-BSD3-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-red.svg)]()

[![Status](https://img.shields.io/badge/status-beta-orange.svg)]()
[![Build Status](https://travis-ci.org/Fortran-FOSS-Programmers/FOODIE.svg?branch=master)](https://travis-ci.org/Fortran-FOSS-Programmers/FOODIE)
[![Coverage Status](https://img.shields.io/codecov/c/github/Fortran-FOSS-Programmers/FOODIE.svg)](http://codecov.io/github/Fortran-FOSS-Programmers/FOODIE?branch=master)
[![GitHub issues](https://img.shields.io/github/issues/Fortran-FOSS-Programmers/FOODIE.svg)]()

### FOODIE, Fortran Object-Oriented Differential-equations Integration Environment

- FOODIE is a pure Fortran (KISS) library providing an awesome environment for the numerical integration of Differential-equations (ODE, PDE);
- FOODIE is Fortran 2008+ standard compliant;
- FOODIE is OOP designed;
- FOODIE is a Free, Open Source Project.

#### Table of Contents

+ [What is FOODIE?](#what-is-FOODIE?)
+ [Main features](#main-features)
+ [Status](#status)
+ [Copyrights](#copyrights)
+ [Documentation](#documentation)

## What is FOODIE?

Modern Fortran standards (2003+) have introduced support for Object-Oriented Programming (OOP). Exploiting new features like Abstract Data Type (ADT) is now possible to develop a KISS library providing an awesome environment for the numerical integration of Differential-equations such as Ordinary and Partial Differential Eqautions (ODE, PDE). FOODIE is tailored to the systems arising from the semi-discretization of PDEs, but it is not limited to them.

The FOODIE environment allows the (numerical) solution of general, non linear differential equations system of the form:

![IVP](https://raw.githubusercontent.com/wiki/Fortran-FOSS-Programmers/FOODIE/images/IVP.png)

where:

+ *U_t = dU/dt*;
+ *U* is the vector of *state* variables being a function of the time-like independent variable *t*;
+ *R* is the (vectorial) residual function, it could be a non linear function of the solution *U* itself;
+ *F* is the (vectorial) initial conditions function.

The FOODIE has two main purposes:

+ for developers devising new schemes for the numerical integrations of differential equations (DE):
    + provide a concise, clear, robust and comprehensive abstract environment by means of which:
        + express the solvers formulae with a very high-level language, it being close as much as possible to their *natual* mathematical formulations; this ensures:
            + clearness, conciseness and robustness;
            + fast-developing;
+ for clients that must solve a differential equations system:
    + provide a simple, standard API for many built-in DE solvers out-of-the-box, thus allowing:
        + fast-solution of new problems;
        + robustness: the same DE solver is applied to different problems, i.e. cross-validation;

Go to [Top](#top)

## Main features

FOODIE is aimed to be a KISS-pure-Fortran library for integrating Ordinary Differential Equations (ODE), it being:

+ [x] Pure Fortran implementation;
+ [x] KISS and user-friendly:
    + [x] simple API, presently based on the Rouson's Abstract Data Type Pattern [8];
    + [ ] easy building and porting on heterogeneous architectures;
+ [ ] comprehensive solvers set out-of-the-box:
    + [x] explicit schemes:
        + [x] [Adams-Bashforth](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_adams_bashforth.html) schemes see [7, 12]:
            + [x] 1 step, namely the forward explicit Euler scheme, 1st order accurate;
            + [x] 2 to 16 steps, 2nd to 16th accurate, respectively;
        + [x] [Euler](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_euler_explicit.html) (forward explicit) scheme, 1st order accurate;
        + [x] [Leapfrog](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_leapfrog.html), 2nd order accurate:
            + [x] unfiltered leapfrog, 2nd order accurate, mostly unstable, see [4];
            + [x] Robert-Asselin filtered leapfrog, 1st order accurate, see [4, 5, 6];
            + [x] Robert-Asselin-Williams filtered leapfrog, 3rd order accurate, see [5, 6];
        + [x] [Linear Multistep Methods, SSP](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_lmm_ssp.html) schemes see [16]:
            + [x] 3 steps, 2nd order accurate;
            + [x] 4 steps, 3rd order accurate;
            + [x] 5 steps, 3rd order accurate;
        + [x] [Linear Multistep Methods, SSP with Variable Step Size (VSS)](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_lmm_ssp_vss.html) schemes see [17]:
            + [x] 2 steps, 2nd order accurate;
            + [x] 3 steps, 2nd order accurate;
            + [x] 3 steps, 3rd order accurate;
            + [x] 4 steps, 3rd order accurate;
            + [x] 5 steps, 3rd order accurate;
        + [ ] [Linear Multistep Runge-Kutta Methods SSP](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_ms_runge_kutta_ssp.html) schemes see [18]:
            + [x] 2 steps, 2 stages, 3rd order accurate;
            + [x] 3 steps, 2 stages, 3rd order accurate;
            + [ ] 4 steps, 3 stages, 5th order accurate;
            + [ ] 3 steps, 6 stages, 5th order accurate;
            + [ ] 3 steps, 5 stages, 6th order accurate;
            + [ ] 3 steps, 7 stages, 7th order accurate;
            + [ ] 4 steps, 7 stages, 7th order accurate;
            + [x] 4 steps, 5 stages, 8th order accurate;
            + [ ] 5 steps, 9 stages, 8th order accurate;
            + [ ] 4 steps, 9 stages, 9th order accurate;
            + [ ] 3 steps, 20 stages, 10th order accurate;
        + [ ] Runge-Kutta schemes:
            + [ ] [low-storage](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_runge_kutta_low_storage.html) schemes, see [1, 2, 3]:
                + [x] 1 stage, namely the forward explicit Euler scheme, 1st order accurate;
                + [ ] 2 stages;
                + [ ] 3 stages;
                + [ ] 4 stages;
                + [x] 5 stages, 4th order accurate, 2N registers, see [3];
                + [x] 6 stages, 4th order accurate, 2N registers, see [9];
                + [x] 7 stages, 4th order accurate, 2N registers, see [9];
                + [x] 12 stages, 4th order accurate, 2N registers, see [10];
                + [x] 13 stages, 4th order accurate, 2N registers, see [10];
                + [x] 14 stages, 4th order accurate, 2N registers, see [10];
            + [ ] [TVD/SSP](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_kutta_tvd_runge.html) schemes, see [1]:
                + [x] 1 stage, namely the forward explicit Euler scheme, 1st order accurate;
                + [x] 2 stages, 2nd order accurate;
                + [x] 3 stages, 3rd order accurate;
                + [ ] 4 stages;
                + [x] 5 stages, 4th order accurate;
            + [ ] [embedded (adaptive)](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_runge_kutta_embeddedd.html) schemes:
                + [x] Heun-Euler, 2 stages, 2nd order accurate;
                + [ ] Runge-Kutta-Fehlberg, 5 stages, 4th order accurate;
                + [x] Runge-Kutta-Cash-Karp, 6 stages, 5th order accurate, see [13];
                + [x] Prince-Dormand, 7 stages, 4th order accurate, see [11];
                + [x] Calvo, 9 stages, 6th order accurate, see [14];
                + [x] Feagin, 17 stages, 10th order accurate, see [15];
    + [ ] implicit schemes:
        + [ ] Runge-Kutta schemes;
        + [x] [Adams-Moulton](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_adams_moulton.html) schemes:
            + [x] 0 step, 1st order accurate;
            + [x] 1 step, 2nd accurate;
            + [x] 2 steps, 3rd accurate;
            + [x] 3 steps, 4th accurate;
        + [x] [Backward Differentiation Formula](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_backward_differentiation_formula.html) schemes:
            + [x] 1 step, namely the backward implicit Euler scheme, 1st order accurate;
            + [x] 2 to 6 steps, 2nd to 6th accurate, respectively;
    + [ ] predictor-corrector schemes:
        + [x] [Adams-Bashforth-Moulton](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_adams_bashforth_moulton.html) schemes:
            + [x] 1 step, AB(1)-AM(0), 1st order accurate;
            + [x] 2 steps, AB(2)-AM(1), 2nd accurate;
            + [x] 3 steps, AB(3)-AM(2), 3rd accurate;
            + [x] 4 steps, AB(4)-AM(3), 4th accurate;
+ [ ] efficient and *non intrusive*:
    + [ ] FOODIE environment is unaware of any eventual parallel paradigms the clients used, but it is proved to preserve high scalability on parallel architectures such as:
        + [x] OpenMP directive-based codes on shared memory multi/many cores architectures;
        + [ ] CoArray Fortran (CAF) based codes for Partitioned Global Address Space (PGAS) programming model;
        + [ ] MPI based code on distributed memory clusters;
        + [ ] GPGPU/accelerators device enabled codes;
+ [x] [Tests-Driven](https://github.com/Fortran-FOSS-Programmers/FOODIE/wiki/Examples) Developed ([TDD](https://en.wikipedia.org/wiki/Test-driven_development)):
+ [x] well documented:
    + [x] clear documentation of schemes implementations, e.g. see [Adams-Bashforth API documentation](http://fortran-foss-programmers.github.io/FOODIE/module/foodie_integrator_adams_bashforth.html);
    + [x] complete [API](http://fortran-foss-programmers.github.io/FOODIE/index.html) reference;
    + [x] comprehensive [wiki](https://github.com/Fortran-FOSS-Programmers/FOODIE/wiki):
+ [x] collaborative developed on [GitHub](https://github.com/Fortran-FOSS-Programmers/FOODIE);
+ [x] [FOSS licensed](https://github.com/Fortran-FOSS-Programmers/FOODIE/wiki/Copyrights);

Any feature request is welcome.

#### Bibliography

[1] *High Order Strong Stability Preserving Time Discretizations*, Gottlieb, S., Ketcheson, D. I., Shu, C.W., Journal of Scientific Computing, vol. 38, N. 3, 2009, pp. 251--289.

[2] *Low-Storage Runge-Kutta Schemes*, J. H. Williamson, Journal of Computational Physics, vol. 35, 1980, pp. 48--56.

[3] *Fourth-Order 2N-Storage Runge-Kutta Schemes*, Mark H. Carpenter, Christopher A. Kennedy, NASA Technical Memorandum 109112, June 1994.

[4] *Numerical methods used in atmospheric models*, Mesinger F. and A. Arakawa, Global Atmospheric Research Programme (GARP), [Technical Report](http://twister.ou.edu/CFD2003/Mesinger_ArakawaGARP.pdf), 1976.

[5] *A Proposed Modification to the Robert-Asselin Time Filter*, Williams, P. D., Mon. Wea. Rev., vol. 137, pp. 2538--2546, 2009, doi: http://dx.doi.org/10.1175/2009MWR2724.1.

[6] *The RAW filter: An improvement to the Robert-Asselin filter in semi-implicit integrations*, Williams, P.D., Monthly Weather Review, vol. 139(6), pages 1996--2007, June 2011.

[7] *Linear multistep method*, [wikipedia article](https://en.wikipedia.org/wiki/Linear_multistep_method).

[8] *Scientific Software Design: The Object-Oriented Way*, Rouson, Damian and Xia, Jim and Xu, Xiaofeng, 2011, ISBN 9780521888134, Cambridge University Press, New York, NY, USA.

[9] *High-accuracy large-step explicit Runge–Kutta (HALE-RK) schemes for computational aeroacoustics*, Vasanth Allampalli and Ray Hixon and M. Nallasamy and Scott D. Sawyer, Journal of Computational Physics, vol. 228, 2009, pp. 3837--3850.

[10] *Efficient low-storage Runge–Kutta schemes with optimized stability regions*, Jens Niegemann and Richard Diehl and Kurt Busch, Journal of Computational Physics, vol. 231, 2012, pp. 364--372.

[11] *A family of embedded Runge-Kutta formulae*, Dormand, J. R.; Prince, P. J. (1980), , Journal of Computational and Applied Mathematics 6 (1): 19--26, doi:10.1016/0771-050X(80)90013-3.

[12] *Cowell Type Numerical Integration As Applied to Satellite Orbit Computation*, J. L. Maury Jr., G. P. Segal, X-553-69-46, April 1969, [NASA-TM-X-63542](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19690017325.pdf).

[13] *A variable order Runge-Kutta method for initial value problems with rapidly varying right-hand sides*, J. R. Cash, A. H. Karp, ACM Transactions on Mathematical Software, vol. 16,  pp. 201--222, 1990, doi:10.1145/79505.79507.

[14] *A New Embedded Pair of Runge-Kutta Formulas of orders 5 and 6*, M. Calvo, J.I. Montijano, L. Randez, Computers & Mathematics with Applications, Volume 20, Issue 1, 1990, Pages 15--24, ISSN 0898-1221, http://dx.doi.org/10.1016/0898-1221(90)90064-Q.

[15] *A tenth-order Runge-Kutta method with error estimate*, Feagin, T., Proceedings of the IAENG Conf. on Scientific Computing. 2007.

[16] *Strong Stability Preserving Runge-Kutta and Multistep Time Discretizations*, S. Gottlieb, D. Ketcheson, C.W. Shu, 2011, 978-981-4289-26-9, doi:10.1142/7498, World Scientific Publishing Co. Pte. Ltd.

[17] *Strong Stability Preserving Explicit Linear Multistep Methods with Variable Step Size*, Y. Hadjimichael, D. Ketcheson, L. Lóczi, A. Németh, 2016, SIAM J. Numer. Anal., 54(5), 2799–2832.DOI:10.1137/15M101717X.

[18] *Explicit Strong Stability Preserving Multistep Runge-Kutta Methods*, C. Bresten, S. Gottlieb, Z. Grant, D. Higgs, D. Ketcheson, A. Nemeth, Mathematics of Computation, 2016, http://dx.doi.org/10.1090/mcom/3115.

Go to [Top](#top)

## Status [![Status](https://img.shields.io/badge/status-beta-orange.svg)]()

FOODIE project is young, but developed with love. Many integrators have been implemented using the Rouson's Abstract Data Type Pattern and tested with complex problems, but the library API is still in beta testing status. Nevertheless, FOODIE is already proven to be able to integrate a wide range of different ODE problems, from pure ODEs (Lorenz and inertial oscillations equations) to complex PDEs (Burgers and Euler equations), see the documentation.

We are searching for Fortraners enthusiast joining our team!

## Copyrights

FOODIE is an open source project, it is distributed under a multi-licensing system:

+ for FOSS projects:
  - [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html);
+ for closed source/commercial projects:
  - [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause);
  - [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause);
  - [MIT](http://opensource.org/licenses/MIT).

Anyone is interest to use, to develop or to contribute to FOODIE is welcome, feel free to select the license that best matches your soul!

More details can be found on [wiki](https://github.com/Fortran-FOSS-Programmers/FOODIE/wiki/Copyrights).

Go to [Top](#top)

## Documentation

Besides this README file the FOODIE documentation is contained into its own [wiki](https://github.com/Fortran-FOSS-Programmers/FOODIE/wiki). Detailed documentation of the API is contained into the [GitHub Pages](http://Fortran-FOSS-Programmers.github.io/FOODIE/index.html) that can also be created locally by means of [ford tool](https://github.com/cmacmackin/ford).

Go to [Top](#top)
