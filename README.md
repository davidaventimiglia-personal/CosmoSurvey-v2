# CosmoSurvey
Cosmology survey modeler

This package contains the CosmoSurvey cosmology survey modeler.

## Dependencies

- **Cuba library** (4.x): Multi-dimensional numerical integration - http://www.feynarts.de/cuba/
- **PLplot** (5.x) with Fortran bindings: Scientific plotting library
  - pkg-config module: `plplot-fortran`
  - On Debian/Ubuntu: `sudo apt install libplplot-dev plplot-driver-cairo`
- **Fortran compiler**: gfortran or similar (requires Fortran 2003 for `iso_c_binding`)
- **C compiler**: for Cuba library integration
- **GNU Autotools**: autoconf, automake (for building from source)

## Building

```bash
# Generate configure script
autoreconf -i

# Configure and build
./configure
make

# Install (optional)
make install
```

## Usage

```bash
cosmosurvey -h                    # Show help
cosmosurvey -z1 0.01 -z2 2.5 -nz 10 -l1 44.0 -l2 46.0 -nl 10 \
            -do 12.0 -fl 1.25e-13 -w -1.0 -om 0.30 -s8 0.8 -g 0.55
```

## Documentation

See the file INSTALL for generic compilation and installation
instructions.

See the section FAQ in the documentation (doc/cosmosurvey.info) for
frequently-asked questions.  The documentation is also available in
PDF and HTML, provided you have a recent version of Texinfo installed:
run "make pdf" or "make html"

Please send bug reports to <dventimi@gmail.com>.  Please include the
version number from `cosmosurvey --version', and a complete,
self-contained test case in each bug report.

If you have questions about using CosmoSurvey and the documentation
does not answer them, please send mail to <dventimi@gmail.com>.

-----

Copyright (C) 2009 David A. Ventimiglia

This file is part of CosmoSurvey, the cosmology survey modeller.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
