AC_INIT([omicsds],[0.0.2])

# Only MacOS and Linux platforms are supported
AC_CHECK_PROG([UNAME], [uname], [yes], [no])
AS_IF([test "$UNAME" = "no"],
      [AC_MSG_ERROR([Only MacOS and Linux platforms are supported])])

AC_MSG_CHECKING([for OS platform])
OS_PLATFORM=`uname -s`
case "${OS_PLATFORM}" in
     Linux*)
         BUILD_PLATFORM=Linux
         ;;
     Darwin*)
         BUILD_PLATFORM=Darwin
         ;;
     *)
         AC_MSG_ERROR([Only MacOS and Linux platforms are supported])
         ;;
esac
AC_MSG_RESULT([$BUILD_PLATFORM])

# find R home and set FLAGS
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
RBIN="${R_HOME}/bin/R"

CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`
CXX11STD=`"${RBIN}" CMD config CXX11STD`
CXXFLAGS=`"${RBIN}" CMD config CXX11FLAGS`

AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

if test -z "${OMICSDS_HOME}"; then
  if test -d "../../install"; then
    OMICSDS_HOME="../../install";
  fi
fi

# Default values
AC_SUBST([OMICSDS_INSTALL_PATH], ${OMICSDS_HOME%/})
AC_SUBST([OMICSDS_LIBS], "-lomicsds")

# Allow for overrides
AC_ARG_WITH([omicsds],
            AS_HELP_STRING([--with-omicsds=[DIR]],[OmicsDS root install path [PREFIX]]),
            [OMICSDS_INSTALL_PATH="${with_omicsds}"],
            [])

OMICSDS_LIBS="${OMICSDS_LIBS} -Wl,-rpath,${OMICSDS_INSTALL_PATH}/lib"

# Check existence of OMICSDS_INSTALL_PATH
AS_IF([test "unset$OMICSDS_INSTALL_PATH" = "unset"],
      [AC_MSG_ERROR([Omicsds installation not found.
        Invoke with either \$OMICSDS_HOME env set to the installation path or
        by explicitly specifying the path with the --with-omicsds argument])])

# Sanity check installation
AC_CHECK_FILES([${OMICSDS_INSTALL_PATH}/include/omicsds/omicsds.h],
               [],
               [AC_MSG_ERROR([OmicsDS include files not found])])

# Try compiling with Omicsds headers
CPPFLAGS="${CPPFLAGS} ${CXX11STD} ${CXXFLAGS} -I${OMICSDS_INSTALL_PATH}/include"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <omicsds/omicsds.h>]], [[]])],[],[AC_MSG_ERROR([unable to compile with omicds/omicsds.h])])

# Check Omicsds installation
LDFLAGS="${LD_FLAGS} -L${OMICSDS_INSTALL_PATH}/lib"
LIBS="${LIBS} -lomicsds"
AC_MSG_CHECKING([whether omicsds works])
AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <omicsds/omicsds.h>
int main() {
  if (OmicsDS::version().size() > 0)
    return 0;
  else
    return -1;
}]])],
               [AC_MSG_RESULT([yes])],
               [AC_MSG_ERROR([Omicsds installation at ${OMICSDS_INSTALL_PATH} is not functional])])

# create and report output
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo
echo "Final src/Makevars"
cat src/Makevars
