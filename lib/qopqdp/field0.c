#include <string.h>
#include <math.h>
#include <qdp_df.h>
#ifdef HAVE_NC1
#include <qdp_f1.h>
#include <qdp_d1.h>
#include <qdp_df1.h>
#endif
#ifdef HAVE_NC2
#include <qdp_f2.h>
#include <qdp_d2.h>
#include <qdp_df2.h>
#endif
#ifdef HAVE_NC3
#include <qdp_f3.h>
#include <qdp_d3.h>
#include <qdp_df3.h>
#endif

#define I(...) __VA_ARGS__
#define S2(a,b) T2(a,b)
#define T2(a,b) a ## b
#define S3(a,b,c) T3(a,b,c)
#define T3(a,b,c) a ## b ## c
#define S4(a,b,c,d) T4(a,b,c,d)
#define T4(a,b,c,d) a ## b ## c ## d
#define S5(a,b,c,d,e) T5(a,b,c,d,e)
#define T5(a,b,c,d,e) a ## b ## c ## d ## e
static char *mtname = "qopqdp." STR(FTYPE);

typedef S3(qopqdp_,FTYPE,_t) ftype_t;
#define ftype_opt S3(qopqdp_,FTYPE,_opt)
#define ftype_as_array_opt_len S3(qopqdp_,FTYPE,_as_array_opt_len)
#define ftype_as_array_opt S3(qopqdp_,FTYPE,_as_array_opt)
#define ftype_gc S3(qopqdp_,FTYPE,_gc)
#define ftype_zero S3(qopqdp_,FTYPE,_zero)
#define ftype_create S3(qopqdp_,FTYPE,_create)
#define ftype_create_unset S3(qopqdp_,FTYPE,_create_unset)
#define ftype_set S3(qopqdp_,FTYPE,_set)
#define ftype_clone S3(qopqdp_,FTYPE,_clone)
#define ftype_read S3(qopqdp_,FTYPE,_read)
#define ftype_write S3(qopqdp_,FTYPE,_write)
#define GET_FTYPE S2(GET_QOPQDP_,FTYPEC)
#define OPT_FTYPE S2(OPT_QOPQDP_,FTYPEC)
#define GET_AS_FTYPE_ARRAY S3(GET_AS_QOPQDP_,FTYPEC,_ARRAY)
#define OPT_AS_FTYPE_ARRAY S3(OPT_AS_QOPQDP_,FTYPEC,_ARRAY)

typedef S3(qopqdp_,FTYPE,O_t) ftypeO_t;
#define OPT_FTYPEO S3(OPT_QOPQDP_,FTYPEC,O)
#define ftypeO_create_unset S3(qopqdp_,FTYPE,O_create_unset)

#ifdef COLORED
#define NCPROTT int NC,
#define NCARGT NC,
#define NCV NC
#define SP(x,y,m,...)					\
  switch(NC) {						\
  case 1: I(S3(x,1_,y)m(qdptype1,__VA_ARGS__)); break;	\
  case 2: I(S3(x,2_,y)m(qdptype2,__VA_ARGS__)); break;	\
  case 3: I(S3(x,3_,y)m(qdptype3,__VA_ARGS__)); break;	\
  default: S2(x,y)(__VA_ARGS__);			\
  }
typedef S2(QDP_1_,QDPT) qdptype1;
typedef S2(QDP_2_,QDPT) qdptype2;
typedef S2(QDP_3_,QDPT) qdptype3;
#else
#define NCPROTT
#define NCARGT
#define NCV 1
#define SP(x,y,m,...) S2(x,y)(__VA_ARGS__)
#endif
#define SP2(a,b,m,...) SP(a,b,m,__VA_ARGS__)
#define SP3(a,b,c,m,...) SP(a,S2(b,c),m,__VA_ARGS__)
#define SP4(a,b,c,d,m,...) SP(a,S3(b,c,d),m,__VA_ARGS__)
#define SP5(a,b,c,d,e,m,...) SP(a,S4(b,c,d,e),m,__VA_ARGS__)
#define SP6(a,b,c,d,e,f,m,...) SP(a,S5(b,c,d,e,f),m,__VA_ARGS__)
#define C1x(t,a,b) ((t*)a,b)
#define Cx1x(t,a,b,c) (a,(t*)b,c)
#define C11x(t,a,b,c) ((t*)a,(t*)b,c)
#define Cxx2x(t,a,b,c,d) (a,b,(t**)c,d)
#define C1x1x(t,a,b,c,d) ((t*)a,b,(t*)c,d)
#define Cx11x(t,a,b,c,d) (a,(t*)b,(t*)c,d)
#define Cx1xx(t,a,b,c,d) (a,(t*)b,c,d)
#define Cx11xx(t,a,b,c,d,e) (a,(t*)b,(t*)c,d,e)
#define SP3x(a,b,c,m,...) S3(a,b,c)(__VA_ARGS__)
#define SP4x(a,b,c,d,m,...) S4(a,b,c,d)(__VA_ARGS__)
#define SP5x(a,b,c,d,e,m,...) S5(a,b,c,d,e)(__VA_ARGS__)

typedef S2(QDP_,QDPT) qdptype;
#ifdef PRECISE
typedef S2(QDP_F_,QDPT) qdptypeF;
typedef S2(QDP_D_,QDPT) qdptypeD;
#endif
#define qlatype S2(QLA_,QDPT)
#define qlazero(f) SP3x(QLA_,A,_eq_zero,Cl,f)
#define qlaeq(f1,f2) SP4x(QLA_,A,_eq_,A,Cll,f1,f2)
#define qdpcreate(l) S3(QDP_create_,A,_L)(l)
#define qdpdestroy(f) S2(QDP_destroy_,A)(f)
#define qdpzero(f,s) SP3(QDP_,A,_eq_zero,C1x,f,s)
#define qdpTeqt(f,c,s) SP4x(QDP_,A,_eq_,AL,C1lx,f,c,s)
#define qdpeq(f1,f2,s) SP4(QDP_,A,_eq_,A,C11x,f1,f2,s)
#define qdpeqs(f1,f2,sh,sd,s) S4(QDP_,A,_eq_s,A)(f1,f2,sh,sd,s)
#define qdpptrreadonly(f,s) S2(QDP_site_ptr_readonly_,A)(f,s)
#define qdpptrreadwrite(f,s) S2(QDP_site_ptr_readwrite_,A)(f,s)
#define qdpgaussian(f,r,s) SP3x(QDP_,A,_eq_gaussian_S,C1x,f,r,s)
#define qdptimesr(f1,r,f2,s) SP4(QDP_,A,_eq_r_times_,A,C1x1x,f1,r,f2,s)
#define qdppeqtimesr(f1,r,f2,s) SP4(QDP_,A,_peq_r_times_,A,C1x1x,f1,r,f2,s)
#define qdptimesc(f1,c,f2,s) SP4(QDP_,A,_eq_c_times_,A,C1x1x,f1,c,f2,s)
#define qdppeqtimesc(f1,c,f2,s) SP4(QDP_,A,_peq_c_times_,A,C1x1x,f1,c,f2,s)
#define qdpsum(r,f,s) SP4x(QDP_,AL,_eq_sum_,A,Cl1x,r,f,s)
#define qdpsummulti(r,f,s,n) SP5x(QDP_,AL,_eq_sum_,A,_multi,Cl1xx,r,f,s,n)
#define qdpnorm2(r,f,s) SP3(QDP_,r_eq_norm2_,A,Cx1x,r,f,s)
#define qdpnorm2multi(r,f,s,n) SP4(QDP_,r_eq_norm2_,A,_multi,Cx1xx,r,f,s,n)
#define qdplnorm2(r,f,s) SP3(QDP_,R_eq_norm2_,A,Cx1x,r,f,s)
#define qdpdot(c,f1,f2,s) SP5(QDP_,c_eq_,A,_dot_,A,Cx11x,c,f1,f2,s)
#define qdpdotmulti(c,f1,f2,s,n) SP6(QDP_,c_eq_,A,_dot_,A,_multi,Cx11xx,c,f1,f2,s,n)
#ifdef ISREAL
#define qdpredot(r,f1,f2,s) SP5(QDP_,r_eq_,A,_dot_,A,Cx11x,r,f1,f2,s)
#define qdpredotmulti(r,f1,f2,s,n) SP6(QDP_,r_eq_,A,_dot_,A,_multi,Cx11xx,r,f1,f2,s,n)
#else
#define qdpredot(r,f1,f2,s) SP5(QDP_,r_eq_re_,A,_dot_,A,Cx11x,r,f1,f2,s)
#define qdpredotmulti(r,f1,f2,s,n) SP6(QDP_,r_eq_re_,A,_dot_,A,_multi,Cx11xx,r,f1,f2,s,n)
#endif
#define qdpvread(r,m,f,n) SP3(QDP_,vread_,A,Cxx2x,r,m,f,n)
#define qdpvwrite(r,m,f,n) SP3(QDP_,vwrite_,A,Cxx2x,r,m,f,n)
#define qdpcreateF(l) S3(QDP_F_create_,A,_L)(l)
#define qdpdestroyF(f) S2(QDP_F_destroy_,A)(f)
#define qdpcreateD(l) S3(QDP_D_create_,A,_L)(l)
#define qdpdestroyD(f) S2(QDP_D_destroy_,A)(f)
#if QOP_Precision == 'F'
#define qdptypeO qdptypeD
#define qdpcreateO qdpcreateD
#define qdpdestroyO qdpdestroyD
#define qdpeqO(f1,f2,s) S4(QDP_FD_,A,_eq_,A)(f1,f2,s)
#define qdpeqsO(f1,f2,sh,sd,s) S4(QDP_D_,A,_eq_s,A)(f1,f2,sh,sd,s)
#define qdpOeq(f1,f2,s) S4(QDP_DF_,A,_eq_,A)(f1,f2,s)
#define qdpvreadO(r,m,f,n) S2(QDP_D_vread_,A)(r,m,f,n)
#define qdpvwriteO(r,m,f,n) S2(QDP_D_vwrite_,A)(r,m,f,n)
#else
#define qdptypeO qdptypeF
#define qdpcreateO qdpcreateF
#define qdpdestroyO qdpdestroyF
#define qdpeqO(f1,f2,s) S4(QDP_DF_,A,_eq_,A)(f1,f2,s)
#define qdpeqsO(f1,f2,sh,sd,s) S4(QDP_F_,A,_eq_s,A)(f1,f2,sh,sd,s)
#define qdpOeq(f1,f2,s) S4(QDP_FD_,A,_eq_,A)(f1,f2,s)
#define qdpvreadO(r,m,f,n) S2(QDP_F_vread_,A)(r,m,f,n)
#define qdpvwriteO(r,m,f,n) S2(QDP_F_vwrite_,A)(r,m,f,n)
#endif

ftype_t *
ftype_opt(lua_State *L, int *idx, int required, ftype_t *def)
{
  ftype_t *t;
  if(required) {
    t = luaL_checkudata(L, *idx, mtname);
    (*idx)++;
  } else {
    t = luaL_testudata(L, *idx, mtname);
    if(t==NULL) t = def;
    else (*idx)++;
  }
  return t;
}

int
ftype_as_array_opt_len(lua_State *L, int idx, int required, int def)
{
  int len = def, tlen = 1, i = idx;
  int type = lua_type(L, idx);
  if(type==LUA_TTABLE) {
    tlen = tableLength(L, idx);
    tableGetIndex(L, idx, 1);
    i = -1;
  }
  ftype_t *t = NULL;
  if(required) {
    t = luaL_checkudata(L, i, mtname);
  } else {
    t = luaL_testudata(L, i, mtname);
  }
  if(t!=NULL) len = tlen;
  if(type==LUA_TTABLE) {
    lua_pop(L, 1);
  }
  return len;
}

void
ftype_as_array_opt(lua_State *L, int *idx, int required,
		   int n, ftype_t **t, int dn, ftype_t **def)
{
  int type = lua_type(L, *idx);
  if(type==LUA_TTABLE) {
    tableGetIndex(L, *idx, 1);
    if(required) {
      t = luaL_checkudata(L, -1, mtname);
    } else {
      t = luaL_testudata(L, -1, mtname);
    }
    lua_pop(L, 1);
    if(t==NULL) {
      for(int i=0; i<n; i++) t[i] = def[i];
    } else {
      for(int i=0; i<n; i++) {
	tableGetIndex(L, *idx, i+1);
	t[i] = luaL_checkudata(L, -1, mtname);
	lua_pop(L, 1);
      }
      (*idx)++;
    }
  } else {
    t[0] = ftype_opt(L, idx, required, def[0]);
  }
}

static int
ftype_gc(lua_State *L)
{
  BEGIN_ARGS;
  GET_FTYPE(t);
  END_ARGS;
  qdpdestroy(t->field);
  return 0;
}

// 1: (optional) field or table of fields of same type (ignored)
// 2: field or table of fields of same type
// 3: reader
// return: field metadata
static int
ftype_read(lua_State *L)
{
#define NC nc
  BEGIN_ARGS;
  OPT_AS_FTYPE_ARRAY(n1,t1,1,(ftype_t*[]){NULL});
  OPT_AS_FTYPE_ARRAY(n2,t2,1,(ftype_t*[]){NULL});
  GET_READER(r);
  END_ARGS;
  ftype_t **t = t2;
  int n = n2;
  if(t2[0]==NULL) { t = t1; n = n1; }
  double dt = -QDP_time();
  QDP_String *md = QDP_string_create();
  int prec, type, nc;
  qopqdp_get_prec_type_nc(r->qr, &prec, &type, &nc);
#ifdef PRECISE
  if(prec==QDP_Precision) {
#endif
    qdptype *q[n];
    for(int i=0; i<n; i++) q[i] = t[i]->field;
    qdpvread(r->qr, md, q, n);
#ifdef PRECISE
  } else {
    qdptypeO *q[n];
    for(int i=0; i<n; i++) q[i] = qdpcreateO(t[0]->qlat);
    qdpvreadO(r->qr, md, q, n);
    for(int i=0; i<n; i++) {
      qdpeqO(t[i]->field, q[i], QDP_all_L(t[0]->qlat));
      qdpdestroyO(q[i]);
    }
  }
#endif
  lua_pushstring(L, QDP_string_ptr(md));
  QDP_string_destroy(md);
  dt += QDP_time();
  printf0("%s: %g seconds\n", __func__, dt);
  return 1;
#undef NC
}

// 1: (optional) field or table of fields of same type (ignored)
// 2: field or table of fields of same type
// 3: writer
// 4: metadata string
// 5: (optional) precision string
// return: field metadata
static int
ftype_write(lua_State *L)
{
#define NC QDP_get_nc(t[0]->field)
  BEGIN_ARGS;
  OPT_AS_FTYPE_ARRAY(n1,t1,1,(ftype_t*[]){NULL});
  OPT_AS_FTYPE_ARRAY(n2,t2,1,(ftype_t*[]){NULL});
  GET_WRITER(w);
  GET_STRING(mds);
#ifdef PRECISE
  OPT_STRING(precision, "F");
#endif
  END_ARGS;
  ftype_t **t = t2;
  int n = n2;
  if(t2[0]==NULL) { t = t1; n = n1; }
  double dt = -QDP_time();
  QDP_String *md = QDP_string_create();
  QDP_string_set(md, (char *)mds);
#ifdef PRECISE
  if(*precision==QDP_Precision) {
#endif
    qdptype *q[n];
    for(int i=0; i<n; i++) q[i] = t[i]->field;
    qdpvwrite(w->qw, md, q, n);
#ifdef PRECISE
  } else {
    qdptypeO *q[n];
    for(int i=0; i<n; i++) {
      q[i] = qdpcreateO(t[0]->qlat);
      qdpOeq(q[i], t[i]->field, QDP_all_L(t[0]->qlat));
    }
    qdpvwriteO(w->qw, md, q, n);
    for(int i=0; i<n; i++) qdpdestroyO(q[i]);
  }
#endif
  QDP_string_destroy(md);
  dt += QDP_time();
  printf0("%s: %g seconds\n", __func__, dt);
  return 0;
#undef NC
}

// create clone of field
static int
ftype_clone(lua_State *L)
{
#define NC QDP_get_nc(t2->field)
  BEGIN_ARGS;
  GET_FTYPE(t2);
#ifdef PRECISE
  OPT_STRING(precision, t2->lat->defaultPrecision);
#endif
  OPT_QSUBSET(sub, t2->lat, QDP_all_L(t2->qlat));
  END_ARGS;
#ifdef PRECISE
  if(*precision == QOP_Precision) {
#endif
    ftype_t *t1 = ftype_create_unset(L, NCARGT t2->lat);
    qdpeq(t1->field, t2->field, sub);
#ifdef PRECISE
  } else {
    ftypeO_t *t1 = ftypeO_create_unset(L, NCARGT t2->lat);
    qdpOeq(t1->field, t2->field, sub);
  }
#endif
  return 1;
#undef NC
}

// set from another field
static int
ftype_set(lua_State *L)
{
#define NC QDP_get_nc(t1->field)
  BEGIN_ARGS;
  GET_FTYPE(t1);
  OPT_FTYPE(t2, NULL);
#ifdef PRECISE
  OPT_FTYPEO(t2o, NULL);
#endif
  END_ARGS;
  lattice_t *lat2 = NULL;
  if(t2) lat2 = t2->lat;
#ifdef PRECISE
  else if(t2o) lat2 = t2o->lat;
#endif
  qassert(lat2!=NULL);
  if(t1->qlat==lat2->qlat) {
    if(t2) {
      qdpeq(t1->field, t2->field, QDP_all_L(t1->qlat));
#ifdef PRECISE
    } else {
      qdpeqO(t1->field, t2o->field, QDP_all_L(t1->qlat));
#endif
    }
  } else {
    // copy with truncation/replication
    QDP_Lattice *rlat = t1->qlat;
    QDP_Lattice *slat = lat2->qlat;
    int rnd = QDP_ndim_L(rlat);
    int snd = QDP_ndim_L(slat);
    //printf0("rls:");
    //for(int i=0; i<rnd; i++) printf0(" %i", QDP_coord_size_L(rlat,i));
    //printf0("\n");
    //printf0("sls:");
    //for(int i=0; i<snd; i++) printf0(" %i", QDP_coord_size_L(slat,i));
    //printf0("\n");
    int roff[rnd], rlen[rnd], sdir[rnd], soff[snd];
    for(int i=0; i<rnd; i++) {
      roff[i] = 0;
      rlen[i] = QDP_coord_size_L(rlat,i);
      sdir[i] = i;
    }
    for(int i=0; i<snd; i++) {
      soff[i] = 0;
    }
#ifdef PRECISE
    qdptypeF *qt2 = NULL;
    if(t2==NULL) {
      qt2 = qdpcreateF(lat2->qlat);
    }
#endif
    for(int s=0; ; s++) {
      printf0("s: %i\n", s);
      QDP_Subset *sub;
      QDP_Shift shift;
      qhmc_qopqdp_getCopyHyper(&shift, &sub, rlat, roff, rlen, sdir, slat, soff, s);
      if(sub==NULL) break;
      //printf0("created subset and map\n");
      if(t2) {
	qdpeqs(t1->field, t2->field, shift, QDP_forward, sub[0]);
#ifdef PRECISE
      } else {
#if QOP_Precision == 'F'
	qdpeqO(qt2, t2o->field, sub[0]);
	qdpeqs(t1->field, qt2, shift, QDP_forward, sub[0]);
#else
	qdpeqsO(qt2, t2o->field, shift, QDP_forward, sub[0]);
	qdpeqO(t1->field, qt2, sub[0]);
#endif
#endif
      }
      //printf0("finished shift\n");
      QDP_destroy_shift(shift);
      //printf0("destroyed map\n");
      QDP_destroy_subset(sub);
      //printf0("destroyed subset\n");
    }
#ifdef PRECISE
    if(t2==NULL) {
      qdpdestroyF(qt2);
    }
#endif
  }
  return 0;
#undef NC
}

#ifdef ISRSTATE
static void
lex_int(QLA_Int *li, int coords[], void *args)
{
  QDP_Lattice *lat = (QDP_Lattice *) args;
  int nd = QDP_ndim_L(lat);
  int t = coords[nd-1];
  for(int i=nd-2; i>=0; i--) {
    t = t*QDP_coord_size_L(lat,i) + coords[i];
  }
  *li = t;
}

void
qhmc_qopqdp_seed_func(QDP_RandomState *r, int seed, int uniform, QDP_Subset s)
{
  QDP_Lattice *qlat = QDP_get_lattice_S(r);
  QDP_Int *li = QDP_create_I_L(qlat);
  if(uniform!=-1) { // uniform seed
    QDP_I_eq_i(li, &seed, s);
  } else {
    QDP_I_eq_funca(li, lex_int, (void*)qlat, s);
  }
  QDP_S_eq_seed_i_I(r, seed, li, s);
  QDP_destroy_I(li);
}

static int
ftype_seed(lua_State *L)
{
  BEGIN_ARGS;
  GET_FTYPE(t);
  GET_INT(seed);
  OPT_INT(uniform, -1);
  OPT_QSUBSET(sub, t->lat, QDP_all_L(t->qlat));
  END_ARGS;
  qhmc_qopqdp_seed_func(t->field, seed, uniform, sub);
  return 0;
}
#endif

#ifdef ARITH
static int
ftype_zero(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_QSUBSET(sub, t->lat, QDP_all_L(t->qlat));
  END_ARGS;
  qdpzero(t->field, sub);
  return 0;
#undef NC
}

static int
ftype_unit(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_QSUBSET(sub, t->lat, QDP_all_L(t->qlat));
  END_ARGS;
  GET_QLA_UNIT(c);
  qdpTeqt(t->field, &c, sub);
  return 0;
#undef NC
}

// 1: field
// 2: coords table
// 3: color/spin values (if type has)
// 4: (optional) real or complex value
// FIXME: only returns value on node containing site
static int
ftype_point(lua_State *L)
{
  BEGIN_ARGS;
  GET_FTYPE(t);
  GET_TABLE_LEN_INDEX(nd,ip);
  GET_COLOR_SPIN;
  OPT_AS_COMPLEX_PTR(z,NULL);
  END_ARGS;
  int rv = 0;
  int site[nd]; qhmc_get_int_array(L, ip, nd, site);
  int node = QDP_node_number(site);
  if(z) {
    if(node==QDP_this_node) {
      int index = QDP_index(site);
      qlatype *q = qdpptrreadwrite(t->field, index);
#ifdef ISREAL
      *q = z->r;
#else
      QLA_c_eq_r_plus_ir(QLAELEM(*q), z->r, z->i);
#endif
    }
  } else {
    rv = 1;
#ifdef ISREAL
    QLA_Real q = 0;
#else
    QLA_Complex q;
    QLA_c_eq_r(q, 0);
#endif
    if(node==QDP_this_node) {
      int index = QDP_index(site);
      qlatype *p = qdpptrreadonly(t->field, index);
#ifdef ISREAL
      q = *p;
#else
      q = QLAELEM(*p);
#endif
    }
#ifdef ISREAL
    sum_real_array(&q, 1);
    lua_pushnumber(L, q);
#else
    sum_real_array((QLA_Real*)&q, 2);
    qhmc_complex_create(L, QLA_real(q), QLA_imag(q));
#endif
  }
  return rv;
}

static int
ftype_site(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  GET_TABLE_LEN_INDEX(nd,ip);
  //OPT_AS_COMPLEX_PTR(v,NULL);
  END_ARGS;
  int rv = 0;
  int site[nd]; qhmc_get_int_array(L, ip, nd, site);
  int node = QDP_node_number(site);
#if 0
  if(v) {
    if(node==QDP_this_node) {
      int index = QDP_index(site);
      qlatype *q = qdpptrreadwrite(t->field, index);
      qlaeq(*q, v);
    }
  } else {
#endif
    rv = 1;
    qlatype q;
    if(node==QDP_this_node) {
      int index = QDP_index(site);
      qlatype *p = qdpptrreadonly(t->field, index);
#ifdef ISREAL
      q = *p;
#else
      qlaeq(&q, p);
#endif
    } else {
#ifdef ISREAL
      q = 0;
#else
      qlazero(&q);
#endif
    }
    sum_real_array((QLA_Real *)&q, sizeof(q)/sizeof(QLA_Real));
    pushqlatype(L, NCARGT &q);
#if 0
  }
#endif
  return rv;
}

// normalized to sigma = 1
static int
ftype_random(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_QOPQDP_QRSTATE(rs, t->lat->rs);
  OPT_DOUBLE(s, 1);
  OPT_QSUBSET(sub, t->lat, QDP_all_L(t->qlat));
  END_ARGS;
  qdpgaussian(t->field, rs, sub);
  if(s!=1) {
    QLA_Real r = s;  
    qdptimesr(t->field, &r, t->field, sub);
  }
  return 0;
#undef NC
}

static int
ftype_randomU1(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_QOPQDP_QRSTATE(rs, t->lat->rs);
  OPT_QSUBSET(sub, t->lat, QDP_all_L(t->qlat));
  END_ARGS;
  qdpgaussian(t->field, rs, sub);
  int i;
  QDP_loop_sites(i, sub, {
      qlatype *x = qdpptrreadwrite(t->field, i);
#ifdef ISREAL
      *x = (*x>=0) ? 1 : -1;
#else
      LOOP_FTYPE_ELEM {
	QLA_Complex z = QLAELEM(*x);
	QLA_Real n = QLA_norm2_c(z);
	if(n==0) {
	  QLA_c_eq_r(QLAELEM(*x), 1);
	} else {
	  n = 1/sqrt(n);
	  QLA_c_eq_r_times_c(QLAELEM(*x), n, z);
	}
      } END_LOOP_FTYPE_ELEM;
#endif
    });
  return 0;
#undef NC
}

static int
ftype_normalize(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_DOUBLE(s, 1);
  OPT_QSUBSET(sub, t->lat, QDP_all_L(t->qlat));
  END_ARGS;
  QLA_Real r;
  qdpnorm2(&r, t->field, sub);
  if(r==0) {
    GET_QLA_UNIT(c);
    qdpTeqt(t->field, &c, sub);
    qdpnorm2(&r, t->field, sub);
  }
  r = s/sqrt(r);
  qdptimesr(t->field, &r, t->field, sub);
  return 0;
#undef NC
}

static int
ftype_lnormalize(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_DOUBLE(s, 1);
  OPT_QSUBSET(sub, t->lat, QDP_all_L(t->qlat));
  END_ARGS;
  int i;
  QDP_loop_sites(i, sub, {
      qlatype *x = qdpptrreadwrite(t->field, i);
#ifdef ISREAL
      *x = (*x>=0) ? s : -s;
#else
      LOOP_FTYPE_ELEM {
	QLA_Complex z = QLAELEM(*x);
	QLA_Real n = QLA_norm2_c(z);
	if(n==0) {
	  QLA_c_eq_r(QLAELEM(*x), s);
	} else {
	  n = s/sqrt(n);
	  QLA_c_eq_r_times_c(QLAELEM(*x), n, z);
	}
      } END_LOOP_FTYPE_ELEM;
#endif
    });
  return 0;
#undef NC
}

static int
ftype_makeGroup(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_INT(g, GROUP_GL);
  OPT_QSUBSET(sub, t->lat, QDP_all_L(t->qlat));
  END_ARGS;
  int i;
  QDP_loop_sites(i, sub, {
      qlatype *x = qdpptrreadwrite(t->field, i);
      qlamakegroup(NCARGT x, g);
    });
  return 0;
#undef NC
}

static int
ftype_combine(lua_State *L)
{
#define NC QDP_get_nc(td->field)
  BEGIN_ARGS;
  GET_FTYPE(td);
  GET_AS_FTYPE_ARRAY(ns,ts);
  GET_TABLE_LEN_INDEX(cl,ci);
  OPT_QSUBSET(sub, td->lat, QDP_all_L(td->qlat));
  END_ARGS;
  for(int i=0; i<ns; i++) {
    lua_pushinteger(L, i+1);
    lua_gettable(L, ci);
    if(lua_type(L,-1)==LUA_TNUMBER) {
      QLA_Real r = lua_tonumber(L, -1);
      if(i==0) { qdptimesr(td->field, &r, ts[0]->field, sub); }
      else { qdppeqtimesr(td->field, &r, ts[i]->field, sub); }
#ifndef ISREAL
    } else { // complex
      qhmc_complex_t *c = qhmc_complex_check(L, -1);
      QLA_Complex z;
      QLA_c_eq_r_plus_ir(z, c->r, c->i);
      if(i==0) { qdptimesc(td->field, &z, ts[0]->field, sub); }
      else { qdppeqtimesc(td->field, &z, ts[i]->field, sub); }
#endif
    }
  }
  return 0;
#undef NC
}

static int
ftype_sum(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_AS_QSUBSET_ARRAY(ns, subs, t->lat, 1, QDP_all_and_empty_L(t->qlat));
  END_ARGS;
  if(ns==1) {
    qlatype r;
    qdpsum(&r, t->field, subs[0]);
    pushqlatype(L, NCARGT &r);
  } else {
    qlatype r[ns];
    qdpsummulti(r, t->field, subs, ns);
    lua_createtable(L, ns, 0);
    for(int i=0; i<ns; i++) {
      pushqlatype(L, NCARGT &r[i]);
      lua_rawseti(L, -2, i+1);
    }
  }
  return 1;
#undef NC
}

static int
ftype_norm2(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_AS_QSUBSET_ARRAY(ns, subs, t->lat, 1, QDP_all_and_empty_L(t->qlat));
  END_ARGS;
  if(ns==1) {
    QLA_Real r;
    qdpnorm2(&r, t->field, subs[0]);
    lua_pushnumber(L, r);
  } else {
    QLA_Real r[ns];
    qdpnorm2multi(r, t->field, subs, ns);
    qhmc_push_real_array(L, ns, r);
  }
  return 1;
#undef NC
}

static int
ftype_lnorm2(lua_State *L)
{
#define NC QDP_get_nc(t->field)
  BEGIN_ARGS;
  GET_FTYPE(t);
  OPT_QOPQDP_REAL(r, NULL);
  OPT_QSUBSET(sub, t->lat, QDP_all_L(t->qlat));
  END_ARGS;
  int rv = 0;
  if(r==NULL) {
    r = qopqdp_real_create_unset(L, t->lat);
    rv = 1;
  }
#ifdef ISREAL
  QDP_R_eq_fabs_R(r->field, t->field, sub);
#else
  qdplnorm2(r->field, t->field, sub);
#endif
  return rv;
#undef NC
}

static int
ftype_redot(lua_State *L)
{
#define NC QDP_get_nc(t1->field)
  BEGIN_ARGS;
  GET_FTYPE(t1);
  GET_FTYPE(t2);
  OPT_AS_QSUBSET_ARRAY(ns, subs, t1->lat, 1, QDP_all_and_empty_L(t1->qlat));
  END_ARGS;
  if(ns==1) {
    QLA_Real r;
    qdpredot(&r, t1->field, t2->field, subs[0]);
    lua_pushnumber(L, r);
  } else {
    QLA_Real r[ns];
    qdpredotmulti(r, t1->field, t2->field, subs, ns);
    qhmc_push_real_array(L, ns, r);
  }
  return 1;
#undef NC
}

#ifndef ISREAL
static int
ftype_dot(lua_State *L)
{
#define NC QDP_get_nc(t1->field)
  BEGIN_ARGS;
  GET_FTYPE(t1);
  GET_FTYPE(t2);
  OPT_AS_QSUBSET_ARRAY(ns, subs, t1->lat, 1, QDP_all_and_empty_L(t1->qlat));
  END_ARGS;
  if(ns==1) {
    QLA_Complex c;
    qdpdot(&c, t1->field, t2->field, subs[0]);
    qhmc_complex_create(L, QLA_real(c), QLA_imag(c));
  } else {
    QLA_Complex c[ns];
    qdpdotmulti(c, t1->field, t2->field, subs, ns);
    qhmc_complex_t cc[ns];
    for(int i=0; i<ns; i++) {
      cc[i].r = QLA_real(c[i]);
      cc[i].i = QLA_imag(c[i]);
    }
    qhmc_push_complex_array(L, ns, cc);
  }
  return 1;
#undef NC
}
#endif

#endif


static struct luaL_Reg ftype_reg[] = {
  { "__gc",       ftype_gc },
  { "read",       ftype_read },
  { "write",      ftype_write },
  { "clone",      ftype_clone },
  { "set",        ftype_set },
#ifdef ISRSTATE
  { "seed",       ftype_seed },
#endif
#ifdef ARITH
  { "zero",       ftype_zero },
  { "unit",       ftype_unit },
  { "point",      ftype_point },
  { "site",       ftype_site },
  { "random",     ftype_random },
  { "randomU1",   ftype_randomU1 },
  { "normalize",  ftype_normalize },
  { "lnormalize", ftype_lnormalize },
  //{ "checkU",   qopqdp_gauge_checkU },
  //{ "checkSU",  qopqdp_gauge_checkSU },
  { "makeGroup",  ftype_makeGroup },
  { "combine",   ftype_combine },
  { "sum",       ftype_sum },
  { "norm2",       ftype_norm2 },
  { "lnorm2",      ftype_lnorm2 },
  //{ "infnorm",    qopqdp_force_infnorm },
  { "reDot",       ftype_redot },
#ifdef ISREAL
  { "dot",         ftype_redot },
#else
  { "dot",         ftype_dot },
#endif
  //{ "gamma",      qopqdp_wquark_gamma },
  //{ "smearGauss", qopqdp_wquark_smearGauss },
#endif
  { NULL, NULL}
};
// det
// shift
// transport
// multiply (R*T,C*T,M*M,M1*M)
// split/combine spin
// nc*CV <-> CM
// check/make: GL,U,H,AH,SL,SU,TGL,TH,TAH

ftype_t *
ftype_create_unset(lua_State *L, NCPROTT lattice_t *lat)
{
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  ftype_t *t = lua_newuserdata(L, sizeof(ftype_t));
  t->lat = lat;
  t->qlat = lat->qlat;
#ifdef COLORED
  if(NC<0) NC = QOPQDP_DEFAULTNC;
#endif
  t->nc = NCV;
  t->field = qdpcreate(lat->qlat);
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, ftype_reg);
  }
  lua_setmetatable(L, -2);
  return t;
}

ftype_t *
ftype_create(lua_State *L, NCPROTT lattice_t *lat)
{
  ftype_t *t = ftype_create_unset(L, NCARGT lat);
#ifdef ARITH
  qdpzero(t->field, QDP_all_L(t->qlat));
#endif
  return t;
}
