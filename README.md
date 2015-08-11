<a name="top"></a>

# FOODiE

[![GitHub tag](https://img.shields.io/github/tag/Fortran-FOSS-Programmers/FOODiE.svg)]() [![Join the chat at https://gitter.im/Fortran-FOSS-Programmers/FOODiE](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/Fortran-FOSS-Programmers/FOODiE?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)]()
[![License](https://img.shields.io/badge/license-BSD2-red.svg)]()
[![License](https://img.shields.io/badge/license-BSD3-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-red.svg)]()

[![Status](https://img.shields.io/badge/status-alpha-orange.svg)]()
[![Build Status](https://travis-ci.org/Fortran-FOSS-Programmers/FOODiE.svg?branch=master)](https://travis-ci.org/Fortran-FOSS-Programmers/FOODiE)
[![Coverage Status](https://img.shields.io/codecov/c/github/Fortran-FOSS-Programmers/FOODiE.svg)](http://codecov.io/github/Fortran-FOSS-Programmers/FOODiE?branch=master)
[![GitHub issues](https://img.shields.io/github/issues/Fortran-FOSS-Programmers/FOODiE.svg)]()

### FOODiE, Fortran Object oriented Ordinary Differential Equations integration library

- FOODiE is a pure Fortran (KISS) library for integrating Ordinary Differential Equations (ODE);
- FOODiE is Fortran 2003+ standard compliant;
- FOODiE is OOP designed;
- FOODiE is a Free, Open Source Project.

#### Table of Contents

+ [What is FOODiE?](#what-is-FOODiE?)
+ [Main features](#main-features)
+ [Status](#status)
+ [Copyrights](#copyrights)
+ [Documentation](#documentation)

## What is FOODiE?

Modern Fortran standards (2003+) have introduced support for Object Oriented Programming. Exploiting new features like Abstract Data Type (ADT) is now possible to develop a KISS library for integrating Ordinary Differential Equations on ADT making the development of new numerical schemes faster, easier and clearer.

Go to [Top](#top)

## Main features

FOODiE is aimed to be a KISS-pure-Fortran library for integrating Ordinary Differential Equations (ODE), it being:

+ [x] Pure Fortran implementation;
+ [x] KISS and user-friendly:
    + [x] simple API, presently based on the Rouson's Abstract Data Type Pattern [1];
    + [ ] easy building and porting on heterogeneous architectures;
+ [ ] comprehensive:
    + [x] explicit Euler scheme, 1st order accurate;
    + [ ] explicit Runge-Kutta schemes:
        + [ ] TVD or SSP schemes, see [2]:
            + [x] 1 stage, namely the forward explicit Euler scheme, 1st order accurate;
            + [x] 2 stages, 2nd order accurate;
            + [x] 3 stages, 3rd order accurate;
            + [ ] 4 stages;
            + [x] 5 stages, 4th order accurate;
    + [ ] implicit Runge-Kutta schemes;
    + [ ] low-storage explicit schemes, see [2,3,3];
    + [ ] Leapfrog, 2nd order accurate:
        + [ ] Robert-Asselin filter, for stability and leapfrog, 1st order accurate;
        + [ ] Robert-Asselin-Williams filter, for stability and leapfrog, 3rd order accurate;
    + [ ] Adams-Bashforth variants of 2nd and 3rd order accurate;
+ [ ] efficient:
    + [ ] high scalability on parallel architectures:
        + [ ] support for shared memory multi/many cores architecture;
        + [ ] support for distributed memory cluster;
        + [ ] support for GPGPU/accelerators device;
+ [ ] well documented:
    + [ ] clear documentation of schemes implementations;
    + [x] complete API reference;
    + [ ] comprehensive wiki:
+ [x] collaborative developed;
+ [x] FOSS licensed;

Any feature request is welcome.

#### Bibliography

[1] *Scientific Software Design: The Object-Oriented Way*, Rouson, Damian and Xia, Jim and Xu, Xiaofeng, 2011, ISBN 9780521888134, Cambridge University Press, New York, NY, USA.

[2] *High Order Strong Stability Preserving Time Discretizations*, Gottlieb, S., Ketcheson, D. I., Shu, C.W., Journal of Scientific Computing, vol. 38, N. 3, 2009, pp. 251-289.

[3] *Low-Storage Runge-Kutta Schemes*, J. H. Williamson, Journal of Computational Physics, vol. 35, 1980, pp. 48--56.

[4] *Fourth-Order 2N-Storage Runge-Kutta Schemes*, Mark H. Carpenter, Christopher A. Kennedy, NASA Technical Memorandum 109112, June 1994.

Go to [Top](#top)

## Status

FOODiE project is just started. A small bunch of integrators have been implemented using the Rouson's Abstract Data Type Pattern, but the library API is not stable.

We are searching for Fortraners enthusiast joining our team.

## Copyrights

FOODiE is an open source project, it is distributed under a multi-licensing system:

+ for FOSS projects:
  - [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html);
+ for closed source/commercial projects:
  - [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause);
  - [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause);
  - [MIT](http://opensource.org/licenses/MIT).

Anyone is interest to use, to develop or to contribute to FOODiE is welcome, feel free to select the license that best matches your soul!

More details can be found on [wiki](https://github.com/Fortran-FOSS-Programmers/FOODiE/wiki/Copyrights).

Go to [Top](#top)

## Documentation

Besides this README file the FOODiE documentation is contained into its own [wiki](https://github.com/Fortran-FOSS-Programmers/FOODiE/wiki). Detailed documentation of the API is contained into the [GitHub Pages](http://Fortran-FOSS-Programmers.github.io/FOODiE/index.html) that can also be created locally by means of [ford tool](https://github.com/cmacmackin/ford).

Go to [Top](#top)
