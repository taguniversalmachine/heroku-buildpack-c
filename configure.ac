# Must init the autoconf setup
# The first parameter is project name
# second is version number
# third is bug report address
AC_INIT([mkrand], [1.0])

# Safety checks in case user overwritten --srcdir
# AC_CONFIG_SRCDIR([command.c])

# Store the auxiliary build tools (e.g., install-sh, config.sub, config.guess)
# in this dir (build-aux)
AC_CONFIG_AUX_DIR([build-aux])

# Init automake, and specify this program use relaxed structures.
# i.e. this program doesn't follow the gnu coding standards, and doesn't have
# ChangeLog, COPYING, AUTHORS, INSTALL, README etc. files.
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Check for C compiler
AC_PROG_CC
# We can add more checks in this section

AC_CHECK_FUNCS(gettimeofday)
AC_CHECK_FUNCS(memset)
AC_CHECK_FUNCS(strtol)
AC_CHECK_HEADERS([inttypes.h])
AC_CHECK_HEADERS([stdint.h])
AC_CHECK_HEADERS([sys/time.h])
AC_C_RESTRICT
AC_FUNC_MALLOC
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Tells automake to create a Makefile
# See https://www.gnu.org/software/automake/manual/html_node/Requirements.html
AC_CONFIG_FILES([Makefile])

# Generate the output
AC_OUTPUT