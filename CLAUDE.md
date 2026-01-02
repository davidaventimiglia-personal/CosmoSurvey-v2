# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

This project uses GNU Autotools. The `configure` script must be generated first:

```bash
# Generate configure script (requires autoconf, automake)
autoreconf -i

# Configure, build, install
./configure
make
make install

# Clean build artifacts
make clean
make distclean  # also removes configure-generated files
```

## Dependencies

- **Cuba library** (1.5+): Numerical integration (Vegas algorithm) - http://www.feynarts.de/cuba/
- **PLplot** with Fortran bindings: Plotting library (pkg-config module: `plplotfortran0`)
- **Fortran compiler**: gfortran or similar
- **C compiler**: for Cuba library integration

## Architecture

CosmoSurvey is a Fortran 90 cosmological survey modeler that calculates expected galaxy cluster counts in redshift-luminosity bins.

### Library (`libcosmo.a`)

Core modules built into a static library:

- **constants.f90**: Physical constants (H100, c, rho_crit, M_scale, L_scale, etc.)
- **quadrature.f90**: Wrapper around QUADPACK routines for adaptive integration (`f_qag`, `f_qagi`)
- **cosmology.f90**: Cosmological calculations
  - `CosmoParams` type: stores model parameters (w, Omega_M_0, sigma_8, etc.)
  - `theta_G`: global parameter storage (thread-private for OpenMP)
  - Functions: `f_E(z)` (Hubble ratio), `f_r(z1,z2)` (comoving distance), `f_dV` (comoving volume), `f_D_L(z)` (luminosity distance)
- **structure.f90**: Structure formation model
  - Growth factor `f_G(z)` using Linder parameterization
  - Mass function `f_dN_dV_dmu` (number density per mass)
  - Mass-observable relations: `f_tau_of_mu` (M-T), `f_lambda_of_mu` (L-M)
  - Covariance matrix `f_Cov` for scatter parameters

### Executables

All executables use the kraken.f90 (M_kracken command-line parser), survey2.f90, and survey3.f90 modules:

- **cosmosurvey**: Main CLI tool - computes cluster counts in z-L bins
- **cosmoplot**: Plotting version using PLplot
- **cosmotest**: Test/diagnostic version

### Key Calculation Flow

`survey3.f90` contains the main calculation pipeline:
1. `f_N()`: Entry point - loops over redshift and luminosity bins
2. `f_dN()`: Cluster count in a single bin (integrates over solid angle)
3. `f_dN_dOmega()`: Integrates over redshift using QUADPACK
4. `f_dN_dV()`: Integrates over mass using Cuba (Cuhre algorithm)
5. `dN_dV_dl_dmu()`: Cuba integrand - mass function convolved with luminosity scatter

## CLI Usage

```bash
cosmosurvey -h                    # Show help
cosmosurvey -z1 0.01 -z2 2.5 -nz 10 -l1 44.0 -l2 46.0 -nl 10 \
            -do 12.0 -fl 1.25e-13 -w -1.0 -om 0.30 -s8 0.8 -g 0.55
```

Key parameters: redshift range (z1, z2, nz bins), luminosity range (l1, l2, nl bins), solid angle (do), flux limit (fl), cosmological parameters (w, om, s8, g).
