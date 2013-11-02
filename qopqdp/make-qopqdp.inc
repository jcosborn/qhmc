#HOST := $(shell hostname)
#BGQ := $(filter cetus% vesta% mira%,$(HOST))
MACH := $(shell uname -m)
BGQ := $(filter ppc64,$(MACH))
ifneq "$(BGQ)" ""
  CC = mpixlc_r
  COPT = -O3 -qinfo=pro -g 
else
  CC = mpicc
  COPT = -O3 -g -std=gnu99 -Wall -fopenmp
endif

HAVENC1 = true
HAVENC2 = true
HAVENC3 = true
HAVENCN = true
HAVEMG = true
#DEFAULTNC = 3

#DEFINES = -DQHMC_REPRO_UNIFORM

#TAG = -bgq
SCIDAC = $(HOME)/lqcd/install
QMPDIR = $(SCIDAC)/qmp$(TAG)
QIODIR = $(SCIDAC)/qio$(TAG)
QLADIR = $(SCIDAC)/qla$(TAG)
QDPDIR = $(SCIDAC)/qdp$(TAG)
QOPDIR = $(SCIDAC)/qopqdp$(TAG)

############# end user settings #####################

ifneq ($(HAVENCN),)
  ifeq ($(COLORS),)
    COLORS = \'N\'
  endif
endif
ifneq ($(HAVENC3),)
  ifeq ($(COLORS),)
    COLORS = 3
  endif
  ifeq ($(DEFAULTNC),)
    DEFAULTNC = 3
  endif
endif
ifneq ($(HAVENC2),)
  ifeq ($(COLORS),)
    COLORS = 2
  endif
  ifeq ($(DEFAULTNC),)
    DEFAULTNC = 2
  endif
endif
ifneq ($(HAVENC1),)
  ifeq ($(COLORS),)
    COLORS = 1
  endif
  ifeq ($(DEFAULTNC),)
    DEFAULTNC = 1
  endif
endif
ifeq ($(COLORS),)
  error "must specify at least one HAVENC*"
endif
ifeq ($(DEFAULTNC),)
  DEFAULTNC = 3
endif
ifneq ($(HAVEMG),)
  HAVENCN = true
endif

INC = -DDEFAULTNC=$(DEFAULTNC)
ifneq ($(HAVENC1),)
  INC += -DHAVENC1
endif
ifneq ($(HAVENC2),)
  INC += -DHAVENC2
endif
ifneq ($(HAVENC3),)
  INC += -DHAVENC3
endif
ifneq ($(HAVENCN),)
  INC += -DHAVENCN
endif
ifneq ($(HAVEMG),)
  INC += -DHAVEMG
endif

IQMP = -I$(QMPDIR)/include
LQMP = -L$(QMPDIR)/lib -lqmp

IQIO = -I$(QIODIR)/include
LQIO = -L$(QIODIR)/lib -lqio -llime

#IQLA = -I$(QLADIR)/include -DQLA_Colors=$(COLORS) -DQLA_Nc=$(DEFAULTNC)
IQLA = -I$(QLADIR)/include -DQLA_Colors=$(COLORS)
#IQDP = -I$(QDPDIR)/include -DQDP_PROFILE -DQDP_Colors=$(COLORS) -DQDP_Nc=$(DEFAULTNC)
IQDP = -I$(QDPDIR)/include -DQDP_PROFILE -DQDP_Colors=$(COLORS)
#IQOP = -I$(QOPDIR)/include -DQOP_Colors=$(COLORS) -DQOP_Nc=$(DEFAULTNC)
IQOP = -I$(QOPDIR)/include -DQOP_Colors=$(COLORS)

LQLA = -L$(QLADIR)/lib
LQDP = -L$(QDPDIR)/lib
LQOP = -L$(QOPDIR)/lib -lqopqdp
ifneq "$(HAVENC1)" ""
  LQLA += -lqla_f1 -lqla_d1 -lqla_q1 -lqla_df1 -lqla_dq1
  LQDP += -lqdp_f1 -lqdp_d1 -lqdp_df1
endif
ifneq "$(HAVENC2)" ""
  LQLA += -lqla_f2 -lqla_d2 -lqla_q2 -lqla_df2 -lqla_dq2
  LQDP += -lqdp_f2 -lqdp_d2 -lqdp_df2
endif
ifneq "$(HAVENC3)" ""
  LQLA += -lqla_f3 -lqla_d3 -lqla_q3 -lqla_df3 -lqla_dq3
  LQDP += -lqdp_f3 -lqdp_d3 -lqdp_df3
endif
ifneq "$(HAVENCN)" ""
  LQLA += -lqla_fn -lqla_dn -lqla_qn -lqla_dfn -lqla_dqn
  LQDP += -lqdp_fn -lqdp_dn -lqdp_dfn
endif
LQLA += -lqla_f -lqla_d -lqla_q -lqla_df -lqla_dq -lqla_int -lqla_random -lqla_cmath
LQDP += -lqdp_f -lqdp_d -lqdp_df -lqdp_int -lqdp_common

QHMCQOPQDP = ../qopqdp/libqhmc_qopqdp.a

QOPQDPINC = $(INC) $(IQOP) $(IQDP) $(IQLA) $(IQIO) $(IQMP) $(DEFINES)
QOPQDPLIBS = $(QHMCQOPQDP) $(LQOP) $(LQDP) $(LQLA) $(LQIO) $(LQMP) -lm
QOPQDPDEPS = $(QHMCQOPQDP)
