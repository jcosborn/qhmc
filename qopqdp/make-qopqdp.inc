CC = mpicc
#COPT = -O3 -std=gnu99 -Wall -Wno-strict-aliasing -fopenmp
COPT = -O3 -g -std=gnu99 -Wall -fopenmp

COLORLIB = 2
NC = 2

ifeq ($(COLORLIB),n)
  COLORS = \'N\'
else
  COLORS = $(COLORLIB)
endif

#DEFINES = -DQHMC_REPRO_UNIFORM

TAG = -mpi
SCIDAC = $(HOME)/lqcd/install
QMPDIR = $(SCIDAC)/qmp$(TAG)
QIODIR = $(SCIDAC)/qio$(TAG)
QLADIR = $(SCIDAC)/qla-cvs$(TAG)
QDPDIR = $(SCIDAC)/qdp$(TAG)
QOPDIR = $(SCIDAC)/qopqdp$(TAG)

IQMP = -I$(QMPDIR)/include
LQMP = -L$(QMPDIR)/lib -lqmp

IQIO = -I$(QIODIR)/include
LQIO = -L$(QIODIR)/lib -lqio -llime

IQLA = -I$(QLADIR)/include -DQLA_Colors=$(COLORS) -DQLA_Nc=$(NC)
LQLA = -L$(QLADIR)/lib -lqla_f -lqla_d -lqla_q -lqla_df -lqla_dq \
 -lqla_f$(COLORLIB) -lqla_d$(COLORLIB) -lqla_q$(COLORLIB) -lqla_dn \
 -lqla_df$(COLORLIB) -lqla_dq$(COLORLIB) -lqla_int -lqla_random -lqla_cmath

IQDP = -I$(QDPDIR)/include -DQDP_PROFILE -DQDP_Colors=$(COLORS) -DQDP_Nc=$(NC)
LQDP = -L$(QDPDIR)/lib -lqdp_d$(COLORLIB) -lqdp_d -lqdp_f$(COLORLIB) -lqdp_f \
 -lqdp_df$(COLORLIB) -lqdp_df -lqdp_int -lqdp_common

IQOP = -I$(QOPDIR)/include -DQOP_Colors=$(COLORS) -DQOP_Nc=$(NC)
LQOP = -L$(QOPDIR)/lib -lqopqdp

LMG = -lqdp_dfn -lqdp_dn -lqdp_fn -lqla_dfn -lqla_dn -lqla_fn

QHMCQOPQDP = ../qopqdp/libqhmc_qopqdp.a

QOPQDPINC = $(IQOP) $(IQDP) $(IQLA) $(IQIO) $(IQMP) $(DEFINES)
QOPQDPLIBS = $(QHMCQOPQDP) $(LQOP) $(LMG) $(LQDP) $(LQLA) $(LQIO) $(LQMP) -lm
QOPQDPDEPS = $(QHMCQOPQDP)
