AC_INIT([pocketfft], [1.0.0])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
AM_MAINTAINER_MODE([enable])

dnl
dnl Needed for linking on Windows.
dnl Protect with m4_ifdef because AM_PROG_AR is required in
dnl autoconf >= 1.12 when using -Wall, but the macro is
dnl absent in old versions of autoconf.
dnl
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

LT_INIT
AC_CONFIG_MACRO_DIR([m4])

dnl
dnl By default, install the headers into a subdirectory of
dnl ${prefix}/include to avoid possible header filename collisions.
dnl
includedir="${includedir}/${PACKAGE_NAME}"

dnl
dnl Enable silent build rules if this version of Automake supports them
dnl
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_IF([test x"AS_VAR_GET(CACHEVAR)" = xyes],
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

dnl
dnl Introduce --enable-native-optimizations command line argument to turn on
dnl -march=native compiler flag, disabled by default.
dnl
AC_ARG_ENABLE(
    [native-optimizations],
    [AS_HELP_STRING([--enable-native-optimizations], [Enable non-portable optimizations for your own CPU by compiling with -march=native @<:@default=no@:>@])]
)

AC_PROG_CC_C99
AS_IF(
    [test "x$enable_native_optimizations" = "xyes"],
    [AX_CHECK_COMPILE_FLAG([-march=native],[CC="$CC -march=native"])],
    dnl
    dnl FIXME: On GCC 4.4, we hit an internal compiler error unless either
    dnl -march=native or -fno-tree-fre is specified.
    dnl
    [
        AS_IF(
            [test "x$GCC" = "xyes" -a "x`$CC -dumpversion | cut -d. -f1,2`" = "x4.4"],
            [AX_CHECK_COMPILE_FLAG([-fno-tree-fre], [CFLAGS="$CFLAGS -fno-tree-fre"])]
        )
    ]
)
AX_CHECK_COMPILE_FLAG([-fno-math-errno],[CFLAGS="$CFLAGS -fno-math-errno"])
AX_CHECK_COMPILE_FLAG([-fno-trapping-math],[CFLAGS="$CFLAGS -fno-trapping-math"])
AX_CHECK_COMPILE_FLAG([-fno-rounding-math],[CFLAGS="$CFLAGS -fno-rounding-math"])
AX_CHECK_COMPILE_FLAG([-fno-signaling-nans],[CFLAGS="$CFLAGS -fno-signaling-nans"])
AX_CHECK_COMPILE_FLAG([-fcx-limited-range],[CFLAGS="$CFLAGS -fcx-limited-range"])

AC_PROG_LIBTOOL

dnl
dnl Create pkgconfig .pc file.
dnl
AX_CREATE_PKGCONFIG_INFO(,,,,[])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
