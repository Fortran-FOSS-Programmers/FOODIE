---
title: Installation
---

# Installation

## Prerequisites

A Fortran 2008+ compliant compiler with a C-like preprocessor (the sources ending in `.F90` are preprocessed). The
`fobos` file provides modes for GNU gfortran, Intel Fortran and OpenCoarrays-GNU. The library and the tests suite are
developed with GNU gfortran on GNU/Linux.

## Download

Clone the repository:

```bash
git clone https://github.com/Fortran-FOSS-Programmers/FOODIE
cd FOODIE
```

### Third-party dependencies

The dependencies are declared in the `[dependencies]` section of the `fobos` file and fetched by FoBiS into
`src/third_party/` (the fetched commits are pinned in `src/third_party/fobos.lock`):

```bash
fobis fetch           # clone the dependencies into src/third_party/
fobis fetch --update  # re-fetch them
```

| Library | Purpose | Needed by |
|---|---|---|
| [PENF](https://github.com/szaghi/PENF) | Portable kind parameters (`I_P`, `R_P`, ...) and number-to-string conversion | the library |
| [FLAP](https://github.com/szaghi/FLAP) | Command line interfaces | the tests suite, the CLI tester, the Euler app |
| [FACE](https://github.com/szaghi/FACE) | ANSI colours of the terminal | FLAP |
| [WenOOF](https://github.com/Fortran-FOSS-Programmers/WenOOF) | WENO interpolators | the linear advection problem of the CLI tester, the Euler app |

The library itself depends only on PENF.

## Build with FoBiS

[FoBiS](https://github.com/szaghi/FoBiS) is the build system used for the development of FOODIE: it resolves the
dependencies among modules by itself, the `fobos` file only lists the build modes.

```bash
pip install FoBiS.py
```

### List the modes

```bash
fobis build --lmodes
```

| Mode | Builds |
|---|---|
| `foodie-static-gnu`, `foodie-static-gnu-debug` | static library, GNU gfortran, release or debug |
| `foodie-shared-gnu`, `foodie-shared-gnu-debug` | shared library, GNU gfortran |
| `foodie-static-intel`, `foodie-static-intel-debug` | static library, Intel Fortran |
| `foodie-shared-intel`, `foodie-shared-intel-debug` | shared library, Intel Fortran |
| `foodie-static-opencoarrays-gnu`, `foodie-static-opencoarrays-gnu-debug` | static library, OpenCoarrays-GNU, preprocessed with `-DCAF` |
| `tests-gnu`, `tests-gnu-debug`, `tests-gnu-coverage` | tests suite into `exe/`, GNU gfortran (coverage instrumented for the last) |
| `tests-intel`, `tests-intel-debug` | tests suite into `exe/`, Intel Fortran |
| `tester-gnu`, `tester-gnu-debug` | the CLI tester into `build/tester/` |

### Build the library

```bash
fobis build --mode foodie-static-gnu   # lib/libfoodie.a
fobis build --mode foodie-shared-gnu   # lib/libfoodie.so
```

The library is placed in `./lib/`, the module files in `./lib/mod/`. The other library modes are listed above.

### Build and run the tests

```bash
fobis rule --ex maketests
```

builds the tests suite, runs it and audits the Runge-Kutta tableaus, see [Testing](./testing).

### Use FOODIE as a dependency of your FoBiS project

Declare FOODIE in the `fobos` file of your project and run `fobis fetch`:

```ini
[dependencies]
deps_dir = src/third_party
FOODIE   = https://github.com/Fortran-FOSS-Programmers/FOODIE
```

## Build with fpm

FOODIE is a library-only [fpm](https://fpm.fortran-lang.org) package depending on PENF; the tests, the tester and the
app are not wired into fpm.

```bash
fpm build
```

To use it from your fpm project add to its `fpm.toml`:

```toml
[dependencies]
FOODIE = { git = "https://github.com/Fortran-FOSS-Programmers/FOODIE" }
```

## Pure or impure integrands

The type-bound procedures of `integrand_object` (operators, assignments, description, dimension) are declared `pure` and
`elemental` by means of the macros of `src/lib/preprocessor_macros.h`:

```c
#ifdef _IMPURE_
#define _PURE_
#define _ELEMENTAL_
#else
#define _PURE_ pure
#define _ELEMENTAL_ elemental
#endif
```

By default they are pure, thus the procedures of your concrete integrand overriding them must be pure too. If your
integrand cannot be pure (it does I/O, calls impure libraries, uses an impure interpolator), build FOODIE, and your
code, with the `_IMPURE_` macro defined, which removes the `pure`/`elemental` prefixes:

```bash
fobis build --mode foodie-static-gnu --preproc "-D_IMPURE_"   # FoBiS
fpm build --flag "-D_IMPURE_"                                 # fpm
```

In a `fobos` mode the same is obtained with `preproc = -D_IMPURE_`: this is how the CLI tester (`tester-gnu`) and the
Euler app (`src/app/euler_1D/fobos`) are built, their WenOOF-based integrands being impure. The residual function `t`
and the local error operator are not constrained to be pure in either build.
