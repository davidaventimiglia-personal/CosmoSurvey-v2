#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([cosmosurvey], [1.0], [dventimi@gmail.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/cosmology.f90])
AC_LANG([Fortran])

# Checks for programs.
AC_PROG_CC			dnl C/C++ compiler
AC_PROG_FC			dnl Fortran compiler
AC_PROG_INSTALL			dnl install script
AC_PROG_RANLIB			dnl for working with archives
PKG_PROG_PKG_CONFIG		dnl pkg-config script

# Checks for libraries
AC_SEARCH_LIBS([Vegas], [cuba])
PKG_CHECK_MODULES([PLPLOT], [plplot-fortran])

AC_CONFIG_FILES([
	Makefile
	src/Makefile
	doc/Makefile
	doc/examples/Makefile])
AC_OUTPUT
