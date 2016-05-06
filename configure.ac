AC_INIT(qhmc, 0.1.5, [osborn at alcf.anl.gov])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR(bin/qhmc.c)
AC_CONFIG_SRCDIR(lib/base/qhmc_base.c)
AC_CONFIG_HEADER(include/qhmc_config_internal.h)
AX_PREFIX_CONFIG_H(include/qhmc_config.h)

dnl --enable-lfs
AC_ARG_ENABLE(
  lfs,
  AC_HELP_STRING([--enable-lfs],
    [include luafilesystem library]),
  AC_SUBST(LFS, [$enableval]),
  AC_SUBST(LFS, [yes])
)
AM_CONDITIONAL(HAVELFS, [test "X$LFS" = Xyes])
if test "X$LFS" = Xyes; then
  AC_DEFINE(HAVE_LFS, [], [include luafilesystem library])
fi

dnl --enable-nc=<ColorLib>
AC_ARG_ENABLE(
  nc,
  AC_HELP_STRING([--enable-nc],
    [make just the libraries necessary for the given color]),
  AC_SUBST(NC, [$enableval]),
  AC_SUBST(NC, [all])
)
if test "X$NC" = X; then NC="all"; fi 
if test "X$NC" = Xn; then NC="N"; fi
if test `expr $NC : '[0-9]*$' \& $NC \> 3` -ne 0; then NC="N"; fi
if test "X$NC" != X1 -a "X$NC" != X2 -a "X$NC" != X3 -a "X$NC" != XN; then NC="all"; fi 
AM_CONDITIONAL(HAVENC1, [test "X$NC" = "X1" -o "X$NC" = "Xall"])
AM_CONDITIONAL(HAVENC2, [test "X$NC" = "X2" -o "X$NC" = "Xall"])
AM_CONDITIONAL(HAVENC3, [test "X$NC" = "X3" -o "X$NC" = "Xall"])
AM_CONDITIONAL(HAVENCN, [test "X$NC" = "XN" -o "X$NC" = "Xall"])
if test "X$NC" = "X1" -o "X$NC" = "Xall"; then
  AC_DEFINE(HAVE_NC1, [], [include Nc=1 libraries])
  colors=1
  defaultnc=1
  #havenc="$havenc -DHAVENC1"
fi
if test "X$NC" = "X2" -o "X$NC" = "Xall"; then
  AC_DEFINE(HAVE_NC2, [], [include Nc=2 libraries])
  colors=2
  defaultnc=2
  #havenc="$havenc -DHAVENC2"
fi
if test "X$NC" = "X3" -o "X$NC" = "Xall"; then
  AC_DEFINE(HAVE_NC3, [], [include Nc=3 libraries])
  colors=3
  defaultnc=3
  #havenc="$havenc -DHAVENC3"
fi
if test "X$NC" = "XN" -o "X$NC" = "Xall"; then
  AC_DEFINE(HAVE_NCN, [], [include Nc=N libraries])
  colors="'N'"
  defaultnc=3
  #havenc="$havenc -DHAVENCN"
fi

dnl --enable-colors=<ColorLib>
AC_ARG_ENABLE(
  colors,
  AC_HELP_STRING([--enable-colors],
    [make all fields use this color library (1,2,3,n)]),
  AC_SUBST(COLORS, [$enableval]),
  : #AC_SUBST(COLORS, [])
)
if test -n "$COLORS"; then
  if test "$COLORS" = n; then
    colors="'N'"
  else
    colors="$COLORS"
    defaultnc="$COLORS"
  fi
fi
if test "X$NC" = "XN" -o "X$NC" = "Xall"; then
  AC_DEFINE(HAVE_MG, [], [include multigrid])
fi
AC_DEFINE_UNQUOTED(QOPQDP_COLORS, [$colors], [default QOPQDP color library])
AC_DEFINE_UNQUOTED(QOPQDP_DEFAULTNC, [$defaultnc], [default QOPQDP Nc])
AC_SUBST(QOPQDP_COLORS, [$colors])
#AC_SUBST(QOPQDP_DEFAULTNC, [$defaultnc])

# CUDA
AC_ARG_WITH(
  cuda,
  AC_HELP_STRING([--with-cuda=<CUDA dir>],[location of CUDA library]),
  [CUDADIR="$withval"; WITH_CUDA=1]
)
#AM_CONDITIONAL(HAVECUDA, [test -n "$WITH_CUDA"])
#if test -n "$WITH_CUDA"; then
#  AC_DEFINE(HAVE_CUDA, [], [include CUDA library])
#fi
#AC_SUBST(CUDADIR)
#AC_SUBST(CUDA_CFLAGS, "-I$CUDADIR/include")
#AC_SUBST(CUDA_LIBS, "-L$CUDADIR/lib -lcudart ")

# QUDA
AC_ARG_WITH(
  quda,
  AC_HELP_STRING([--with-quda=<QUDA dir>],[location of QUDA library]),
  [QUDADIR="$withval"; WITH_QUDA=1]
)
AM_CONDITIONAL(HAVEQUDA, [test -n "$WITH_QUDA"])
if test -n "$WITH_QUDA"; then
  if test -z "$WITH_CUDA"; then
    echo "error: specifying --with-quda requires --with-cuda"
    exit 1
  fi
  AC_DEFINE(HAVE_QUDA, [], [include QUDA library])
fi
AC_SUBST(QUDADIR)
AC_SUBST(QUDA_CFLAGS, "-I$QUDADIR/include")
AC_SUBST(QUDA_LIBS, "-L$QUDADIR/lib -lquda -L$CUDADIR/lib64 -lcudart -lstdc++")

# QOPQDP
AC_ARG_WITH(
  qopqdp,
  AC_HELP_STRING([--with-qopqdp=<QOPQDP dir>],[location of QOPQDP libraries]),
  [QOPQDPDIR="$with_qopqdp"; WITH_QOPQDP=1]
)
AC_SUBST(QOPQDPDIR)
AM_CONDITIONAL(HAVEQOPQDP, [test -n "$WITH_QOPQDP"])
if test -n "$WITH_QOPQDP"; then
  AC_DEFINE(HAVE_QOPQDP, [], [include QOPQDP library])
fi

AC_ARG_WITH(
  qdp,
  AC_HELP_STRING([--with-qdp=<QDP dir>],[location of QDP libraries]),
  [QDPDIR="$with_qdp"; WITH_QDP=1]
)
AC_SUBST(QDPDIR)

AC_ARG_WITH(
  qla,
  AC_HELP_STRING([--with-qla=<QLA dir>],[location of QLA libraries]),
  [QLADIR="$with_qla"; WITH_QLA=1]
)
AC_SUBST(QLADIR)

AC_ARG_WITH(
  qio,
  AC_HELP_STRING([--with-qio=<QIO dir>],[location of QIO libraries]),
  [QIODIR="$with_qio"; WITH_QIO=1]
)
AC_SUBST(QIODIR)

AC_ARG_WITH(
  qmp,
  AC_HELP_STRING([--with-qmp=<QMP dir>],[location of QMP libraries]),
  [QMPDIR="$with_qmp"; WITH_QMP=1]
)
AC_SUBST(QMPDIR)

qopqdpcflags="-I${QOPQDPDIR}/include"
qdpcflags="-I${QDPDIR}/include"
qlacflags="-I${QLADIR}/include"
qiocflags="-I${QIODIR}/include"
qopqdplibs="-L${QOPQDPDIR}/lib -lqopqdp"
qdplibs="-L${QDPDIR}/lib"
qlalibs="-L${QLADIR}/lib"
qiolibs="-L${QIODIR}/lib -lqio -llime"
if test "X$NC" = "X2" -o "X$NC" = "Xall"; then
  qdplibs="$qdplibs -lqdp_df2 -lqdp_f2 -lqdp_d2"
  qlalibs="$qlalibs -lqla_df2 -lqla_f2 -lqla_d2"
fi
if test "X$NC" = "X3" -o "X$NC" = "Xall"; then
  qdplibs="$qdplibs -lqdp_df3 -lqdp_f3 -lqdp_d3"
  qlalibs="$qlalibs -lqla_df3 -lqla_f3 -lqla_d3"
fi
if test "X$NC" = "XN" -o "X$NC" = "Xall"; then
  qdplibs="$qdplibs -lqdp_dfn -lqdp_fn -lqdp_dn"
  qlalibs="$qlalibs -lqla_dfn -lqla_fn -lqla_dn"
fi
qdplibs="$qdplibs -lqdp_df -lqdp_f -lqdp_d -lqdp_int -lqdp_common"
qlalibs="$qlalibs -lqla_df -lqla_f -lqla_d -lqla_int -lqla_random -lqla_cmath -lm"

# get QMP configuration
if test "X${QMPDIR}X" = "XX" ; then
  AC_PATH_PROG([QMP_CONFIG], [qmp-config], [], [${PATH}])
else
  AC_PATH_PROG([QMP_CONFIG], [qmp-config], [], [${QMPDIR}/bin ${PATH}])
fi
if test "X${QMP_CONFIG}X" != "XX" ; then
  AC_MSG_NOTICE([Found QMP configuration program ${QMP_CONFIG}])
  # use the QMP CC unless it is explicitly set by the user
  if ! set | grep -q '^CC=' ; then
    AC_SUBST(CC,  "`${QMP_CONFIG} --cc`")
    AC_MSG_NOTICE([QMP compiler: ${CC}])
  fi
  qmpcflags="`${QMP_CONFIG} --cflags`"
  AC_MSG_NOTICE([QMP compile flags: $qmpcflags])
  qmpldflags="`${QMP_CONFIG} --ldflags`"
  AC_MSG_NOTICE([QMP linking flags: $qmpldflags])
  qmplibs="`${QMP_CONFIG} --libs`"
  AC_MSG_NOTICE([QMP libraries flags: $qmplibs])
else
  AC_MSG_WARN([QMP configuration program qmp-config not found.])
  AC_MSG_WARN([Set environment variables QMP_CFLAGS QMP_LDFAGS QMP_LIBS
    before configure])
fi

AC_SUBST(QOPQDP_CFLAGS, "$qopqdpcflags $qdpcflags $qlacflags $qiocflags $qmpcflags")
AC_SUBST(QOPQDP_LIBS, "$qopqdplibs $qdplibs $qlalibs $qiolibs $qmpldflags $qmplibs")

AC_SUBST(LUA_CFLAGS, "-DLUA_USE_POSIX")

# Checks for programs.
AC_PROG_LN_S
AC_PROG_RANLIB
#AC_DISABLE_SHARED
#AC_PROG_LIBTOOL
AC_PROG_CC([c99 cc gcc])
AC_CHECK_TOOL(AR, ar, [ar])
AM_PROG_CC_C_O
#AC_PROG_CC_C99

# pass configuration to distcheck target
DISTCHECK_CONFIGURE_FLAGS="CC=\"$CC\" CFLAGS=\"$CFLAGS\" LDFLAGS=\"$LDFLAGS\" LIBS=\"$LIBS\""
if test -n "$WITH_QMP"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qmp=$QMPDIR"
fi
if test -n "$WITH_QIO"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qio=$QIODIR"
fi
if test -n "$WITH_QLA"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qla=$QLADIR"
fi
if test -n "$WITH_QDP"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qdp=$QDPDIR"
fi
if test -n "$WITH_QOPQDP"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qopqdp=$QOPQDPDIR"
fi
AC_SUBST(DISTCHECK_CONFIGURE_FLAGS)

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(lib/lua/Makefile)
AC_CONFIG_FILES(lib/lua/src/Makefile)
AC_CONFIG_FILES(lib/base/Makefile)
AC_CONFIG_FILES(lib/lbc/Makefile)
AC_CONFIG_FILES(lib/luafilesystem/Makefile)
AC_CONFIG_FILES(lib/qopqdp/Makefile)
AC_CONFIG_FILES(lib/quda/Makefile)
AC_CONFIG_FILES(bin/Makefile)
AC_CONFIG_FILES(test/Makefile)
AC_OUTPUT
