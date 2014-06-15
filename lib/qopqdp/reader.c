#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.reader";

reader_t *
qopqdp_reader_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  reader_t *r = lua_touserdata(L, idx);
#if 1
  int hasmt = lua_getmetatable(L, idx);
  qassert(hasmt==1);
  luaL_getmetatable(L, mtname);
  int eq = lua_equal(L, -1, -2);
  qassert(eq==1);
  lua_pop(L, 2);
#endif
  return r;
}

static int
qopqdp_reader_close(lua_State *L)
{
  BEGIN_ARGS;
  GET_READER(r);
  END_ARGS;
  if(r->open) {
    QDP_close_read(r->qr);
    r->open = 0;
  }
  return 0;
}

/* get QIO record precision */
void
qopqdp_get_prec_type_nc(QDP_Reader *qr, int *prec, int *type, int *nc)
{
  QIO_RecordInfo *ri = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  QDP_String *md = QDP_string_create();
  QDP_read_qio_record_info(qr, ri, md);
  *prec = *QIO_get_precision(ri);
  *nc = QIO_get_colors(ri);
  char *typestr = QIO_get_datatype(ri);
  *type = -1;
  int k = 0;
  while(typestr[k++]!='_'); // go to first char after _
  switch(typestr[k]) {
  case 'R': *type = 'S'; break; // QDP_RandomState
  case 'I': *type = 'I'; break; // QDP_Int
  case 'F': // QDP_F...
  case 'D': // QDP_D...
    switch(typestr[k+2]) {
    case 'R': *type = 'R'; break; // QDP_P_Real
    case 'C': *type = 'C'; break; // QDP_P_Complex
    case '_': // QDP_PC_...
      switch(typestr[k+3]) {
      case 'C': // QDP_PC_Color...
	switch(typestr[k+8]) {
	case 'V': *type = 'V'; break; // QDP_PC_ColorVector
	case 'M': *type = 'M'; break; // QDP_PC_ColorMatrix
	}
	break;
      case 'H': *type = 'H'; break; // QDP_PC_HalfFermion
      case 'D': // QDP_PC_Dirac...
	switch(typestr[k+8]) {
	case 'F': *type = 'D'; break; // QDP_PC_DiracFermion
	case 'P': *type = 'P'; break; // QDP_PC_DiracPropagator
	}
	break;
      }
      break;
    }
    break;
  }
  if(*type==-1) {
    printerr0("unknown datatype: %s\n", typestr);
    ABORT(-1);
  }
  QIO_destroy_record_info(ri);
  QDP_string_destroy(md);
}

// 1: reader
// 2: field or table of fields of same type
// 3: (need to add) optional options table
// return: field metadata
static int
qopqdp_reader_read(lua_State *L)
{
#define NC nc
  BEGIN_ARGS;
  GET_READER(r);
  GET_AS_QOPQDP_DFERMION_ARRAY(nfields, wq);
  END_ARGS;
  QDP_Lattice *qlat = r->lat->qlat;
  //QDP_set_read_group_size(8);
  //printf0("saving lattice file %s\n", fn);
  double dt = -QDP_time();
  QDP_String *md = QDP_string_create();
  int prec, type, nc;
  qopqdp_get_prec_type_nc(r->qr, &prec, &type, &nc);
#if QDP_Precision == 'F'
#define QDPO(x) QDP_D_ ## x
#define QDPPO(x) QDP_FD_ ## x
#else
#define QDPO(x) QDP_F_ ## x
#define QDPPO(x) QDP_DF_ ## x
#endif
  if(type=='D') {
    if(prec==QDP_Precision) {
      QDP_DiracFermion *df[nfields];
      for(int i=0; i<nfields; i++) df[i] = wq[i]->field;
      QDP_vread_D(r->qr, md, df, nfields);
    } else {
      QDPO(DiracFermion) *df[nfields];
      for(int i=0; i<nfields; i++) df[i] = QDPO(create_D_L)(qlat);
      QDPO(vread_D)(r->qr, md, df, nfields);
      for(int i=0; i<nfields; i++) {
	QDPPO(D_eq_D)(wq[i]->field, df[i], QDP_all_L(qlat));
	QDPO(destroy_D)(df[i]);
      }
    }
  } else if(type=='P') {
    qassert(nfields==QLA_Nc*QLA_Ns);
    if(prec==QDP_Precision) {
      QDP_DiracPropagator *dp = QDP_create_P_L(qlat);
      QDP_read_P(r->qr, md, dp);
      for(int i=0; i<nfields; i++) {
	int color = i/QLA_Ns;
	int spin = i%QLA_Ns;
	QDP_D_eq_diracvec_P(wq[i]->field, dp, color, spin, QDP_all_L(qlat));
      }
      QDP_destroy_P(dp);
    } else {
      QDPO(DiracPropagator) *dp = QDPO(create_P_L)(qlat);
      QDPO(read_P)(r->qr, md, dp);
      QDPO(DiracFermion) *df = QDPO(create_D_L)(qlat);
      for(int i=0; i<nfields; i++) {
	int color = i/QLA_Ns;
	int spin = i%QLA_Ns;
	QDPO(D_eq_diracvec_P)(df, dp, color, spin, QDP_all_L(qlat));
	QDPPO(D_eq_D)(wq[i]->field, df, QDP_all_L(qlat));
      }
      QDPO(destroy_D)(df);
      QDPO(destroy_P)(dp);
    }
  }
#undef QDPO
#undef QDPPO
  lua_pushstring(L, QDP_string_ptr(md));
  QDP_string_destroy(md);
  dt += QDP_time();
  printf0(" read in %g seconds\n", dt);
  return 1;
#undef NC
}

static struct luaL_Reg reader_reg[] = {
  { "__gc",    qopqdp_reader_close },
  { "close",   qopqdp_reader_close },
  { "read",    qopqdp_reader_read },
  { NULL, NULL}
};

reader_t *
qopqdp_reader_create(lua_State *L, const char *fn, lattice_t *lat)
{
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  reader_t *r = lua_newuserdata(L, sizeof(reader_t));
  QDP_String *md = QDP_string_create();
  r->qr = QDP_open_read_L(lat->qlat, md, (char*)fn);
  r->lat = lat;
  r->open = 1;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, reader_reg);
  }
  lua_setmetatable(L, -2);
  lua_pushstring(L, QDP_string_ptr(md));
  QDP_string_destroy(md);
  return r;
}
