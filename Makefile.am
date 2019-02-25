ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libpocketfft.la

libpocketfft_la_SOURCES = pocketfft.c pocketfft.h

include_HEADERS = pocketfft.h

EXTRA_DIST = README.md LICENSE.md runtest.sh

check_PROGRAMS = ffttest
ffttest_SOURCES = ffttest.c
ffttest_LDADD = libpocketfft.la

TESTS = runtest.sh

AM_CFLAGS = @AM_CFLAGS@

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
