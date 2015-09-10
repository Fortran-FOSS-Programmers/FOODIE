<a name="top"></a>

# WenOOF

[![GitHub tag](https://img.shields.io/github/tag/Fortran-FOSS-Programmers/WenOOF.svg)]() [![Join the chat at https://gitter.im/Fortran-FOSS-Programmers/WenOOF](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/Fortran-FOSS-Programmers/WenOOF?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)]()
[![License](https://img.shields.io/badge/license-BSD2-red.svg)]()
[![License](https://img.shields.io/badge/license-BSD3-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-red.svg)]()

[![Status](https://img.shields.io/badge/status-alpha-orange.svg)]()
[![Build Status](https://travis-ci.org/Fortran-FOSS-Programmers/WenOOF.svg?branch=master)](https://travis-ci.org/Fortran-FOSS-Programmers/WenOOF)
[![Coverage Status](https://img.shields.io/codecov/c/github/Fortran-FOSS-Programmers/WenOOF.svg)](http://codecov.io/github/Fortran-FOSS-Programmers/WenOOF?branch=master)
[![GitHub issues](https://img.shields.io/github/issues/Fortran-FOSS-Programmers/WenOOF.svg)]()

### WenOOF, WENO interpolation Object Oriented Fortran library

- WenOOF is a pure Fortran (KISS) library for computing WENO interpolations;
- WenOOF is Fortran 2003+ standard compliant;
- WenOOF is OOP designed;
- WenOOF is a Free, Open Source Project.

#### Table of Contents

+ [What is WenOOF?](#what-is-wenoof?)
	+ [What is WENO?](#what-is-weno?)
+ [Main features](#main-features)
+ [Status](#status)
+ [Copyrights](#copyrights)
+ [Documentation](#documentation)

## What is WenOOF?

Modern Fortran standards (2003+) have introduced support for Object Oriented Programming. Exploiting new features like Abstract Data Type (ADT) is now possible to develop a KISS library for computing Weighted Essentially Non-Oscillatory (WENO) interpolation on ADT making the development of new numerical schemes faster, easier and clearer.

### What is WENO?

Starting from the original paper of Liu, Osher and Chan<sup>1</sup> WENO interpolation schemes have gained attention mainly for solving hyperbolic partial differential equations (PDEs) whose solutions admit strong discontinuities as well complex smooth solution features. Consequently, WENO schemes are a common build-block of nonlinear, conservative finite volume methods, however they are not strictly related to the PDEs solution, they being a general, **non linear interpolation (approximation) procedure**.

A clear, yet brief online introduction of WENO schemes family can be found at Prof. Shu's WENO methods page hosted on [scholarpedia](http://www.scholarpedia.org/article/WENO_methods).

Since 1994 the WENO literature has _blowing up_, a superficial search on [sciencedirect](http://www.sciencedirect.com/) for `weno scheme` resulting in more than 1500 matches. During the last 2 decades many new WENO schemes have been proposed: the efficient implementation of Jiang and Shu<sup>2</sup>, the hybrid Compact-WENO scheme of Pirozzoli<sup>3</sup>, the bandwidth-optimized WENO scheme of Martin et al.<sup>4</sup>, the WENO-Z scheme of Borges et al.<sup>5</sup> and many others.

WenOOF is designed to provide a KISS, Object Oriented Fortran API for computing WENO interpolations accordingly the main relevant WENO schemes ever devised and with a new ones we will develop :-)

#### Cited references

[1] _Weighted Essentially Non-oscillatory Schemes_, Xu-Dong Liu, Stanley Osher, Tony Chan, JCP, 1994, vol. 115, pp. 200--212, doi:10.1006/jcph.1994.1187

[2] _Efficient Implementation of Weighted ENO Schemes_, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130

[3] _Conservative Hybrid Compact-WENO Schemes for Shock-Turbulence Interaction_, Sergio Pirozzoli, JCP, 2002, vol. 178, pp. 81--117, doi:10.1006/jcph.2002.7021

[4] _A bandwidth-optimized WENO scheme for the effective direct numerical simulation of compressible turbulence_, M.P. Martín, E.M. Taylor, M. Wu, V.G. Weirs, JCP, 2006, vol. 220, pp. 270--289, doi:10.1016/j.jcp.2006.05.009

[5] _An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws_, Rafael Borges, Monique Carmona, Bruno Costa, Wai Sun Don, JCP, 2007, vol. 227, pp. 3191--3211, doi:10.1016/j.jcp.2007.11.038

Go to [Top](#top)

## Main features

WenOOF is aimed to be a KISS-pure-Fortran library for computing WENO interpolation, it being:

+ [x] Pure Fortran implementation;
+ [x] KISS and user-friendly:
  + [x] simple API;
  + [ ] easy building and porting on heterogeneous architectures;
+ [ ] comprehensive:
  + [ ] central schemes;
  + [x] upwind biased schemes;
  + [ ] hybrid schemes;
+ [ ] efficient:
  + [ ] high scalability on parallel architectures:
    + [ ] support for shared memory multi/many cores architecture;
    + [ ] support for distributed memory cluster;
    + [ ] support for GPGPU/accelerators device;
+ [ ] well documented:
  + [x] clear documentation of schemes implementations;
  + [x] complete API reference;
  + [ ] comprehensive wiki:
+ [ ] collaborative developed;
+ [x] FOSS licensed;

Any feature request is welcome.

Go to [Top](#top)

## Copyrights

WenOOF is an open source project, it is distributed under a multi-licensing system:

+ for FOSS projects:
  - [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html);
+ for closed source/commercial projects:
  - [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause);
  - [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause);
  - [MIT](http://opensource.org/licenses/MIT).

Anyone is interest to use, to develop or to contribute to WenOOF is welcome, feel free to select the license that best matches your soul!

More details can be found on [wiki](https://github.com/Fortran-FOSS-Programmers/WenOOF/wiki/Copyrights).

### Externals libraries

WenOOF uses some external libraries (placed into the *external* subdirectory of the root project) for the testing suite. These library maybe distributed under different licensing system with respect the WenOOF one, please refer to their own licenses.

Go to [Top](#top)

## Documentation

Besides this README file the WenOOF documentation is contained into its own [wiki](https://github.com/Fortran-FOSS-Programmers/WenOOF/wiki). Detailed documentation of the API is contained into the [GitHub Pages](http://Fortran-FOSS-Programmers.github.io/WenOOF/index.html) that can also be created locally by means of [ford tool](https://github.com/cmacmackin/ford).

Go to [Top](#top)
