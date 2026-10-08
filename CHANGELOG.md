# Changelog

All notable changes to this project are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.4.0] — 2026-10-08
### Documentation
- Migrate documentation to VitePress, correct scheme orders


### Fixed
- **adams-moulton**: Make AM-0 backward Euler and fix ABM-1 corrector

- **runge-kutta**: Refine SSP(5,4), LS(6,4) and LS(7,4) coefficients

- **runge-kutta-emd**: Fix Heun-Euler weights and step controller exponent

- **ms-runge-kutta-ssp**: Evaluate stages at their own times, fix fast mode

- **runge-kutta-lssp**: Evaluate stages at their own times


### build
- Manage dependencies with fobis fetch, rewrite fobos ⚠ BREAKING CHANGE


## [0.0.4] — 2015-09-15
### Security
- Try to fix Travis CI token security issue



