ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libpocketfft.la

src_pocketfft = \
  pocketfft.h \
  pocketfft.c

include_HEADERS = \
  pocketfft.h

libpocketfft_la_SOURCES = $(src_pocketfft)

check_PROGRAMS = ffttest
ffttest_SOURCES = ffttest.c
ffttest_LDADD = libpocketfft.la -lm

TESTS = ffttest

AM_CFLAGS = -I$(top_srcdir)

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
