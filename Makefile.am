ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libpocketfft.la

libpocketfft_la_SOURCES = pocketfft/pocketfft.c pocketfft/pocketfft.h

# format is "current:revision:age"
# any change: increase revision
# any interface change: increase current, revision=0
# any backward-compatible change: increase age
# any backward-incompatible change: age=0
# ==> age <= current
libpocketfft_la_LDFLAGS = -version-info 0:0:0

nobase_include_HEADERS = pocketfft/pocketfft.h

EXTRA_DIST = README.md LICENSE.md runtest.sh

check_PROGRAMS = ffttest
ffttest_SOURCES = ffttest.c
ffttest_LDADD = libpocketfft.la

TESTS = runtest.sh

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
