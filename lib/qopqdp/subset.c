#include "qhmc_qopqdp_common.h"
#include <math.h>

static char *mtname = "qopqdp.subset";

subsetGroup_t *
qhmc_qopqdp_subsetGroup_create(lua_State *L, QDP_Subset *qgroup)
{
  subsetGroup_t *sg = (subsetGroup_t *) malloc(sizeof(subsetGroup_t));
  sg->qgroup = qgroup;
  sg->refcount = 0;
  return sg;
}

void
qhmc_qopqdp_subsetGroup_free(lua_State *L, subsetGroup_t *sg)
{
  if(sg) {
    QDP_destroy_subset(sg->qgroup);
    free(sg);
  }
}

subset_t *
qhmc_qopqdp_opt_subset(lua_State *L, int *idx, int required, subset_t *def)
{
  subset_t *lat;
  if(required) {
    lat = luaL_checkudata(L, *idx, mtname);
    (*idx)++;
  } else {
    lat = luaL_testudata(L, *idx, mtname);
    if(lat==NULL) lat = def;
    else (*idx)++;
  }
  return lat;
}

int
qhmc_qopqdp_opt_as_subset_array_len(lua_State *L, int idx, int required,
				    int def)
{
  int len = def, tlen = -1, i = idx;
  int type = lua_type(L, idx);
  if(type==LUA_TTABLE) {
    tlen = tableLength(L, idx);
    tableGetIndex(L, idx, 1);
    i = -1;
  }
  subset_t *t = NULL;
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
qhmc_qopqdp_opt_as_subset_array(lua_State *L, int *idx, int required,
				int n, subset_t **t, int dn, subset_t **def)
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
      for(int i=0; i<abs(n); i++) t[i] = def[i];
    } else {
      for(int i=0; i<n; i++) {
        tableGetIndex(L, *idx, i+1);
        t[i] = luaL_checkudata(L, -1, mtname);
        lua_pop(L, 1);
      }
      (*idx)++;
    }
  } else {
    t[0] = qhmc_qopqdp_opt_subset(L, idx, required, def[0]);
  }
}

static int
qopqdp_subset_gc(lua_State *L)
{
  BEGIN_ARGS;
  GET_LSUBSET(s);
  END_ARGS;
  if(s->group) {
    if(--s->group->refcount==0) {
      qhmc_qopqdp_subsetGroup_free(L, s->group);
    }
  }
  return 0;
}

static int
slice_func(QDP_Lattice *lat, int x[], void *args)
{
  int dir = *(int *)args;
  return x[dir];
}

QDP_Subset *
qhmcqdp_get_timeslices(lattice_t *lat)
{
  if(lat->timeslices==NULL) {
    int dir = lat->nd - 1;
    int n = QDP_coord_size_L(lat->qlat, dir);
    lat->timeslices =
      QDP_create_subset_L(lat->qlat, slice_func, (void*)&dir, sizeof(dir), n);
  }
  return lat->timeslices;
}

static int
staggered_func(QDP_Lattice *lat, int x[], void *args)
{
  int nd = *(int *)args;
  int s = 0;
  for(int i=nd-1; i>=0; i--) {
    s *= 2;
    s += x[i]&1;
  }
  return s;
}

QDP_Subset *
qhmcqdp_get_staggered(lattice_t *lat)
{
  if(lat->staggered==NULL) {
    int nd = QDP_ndim_L(lat->qlat);
    int n = 1 << nd;
    //printf("creating staggered: %i %i\n", nd, n);
    lat->staggered =
      QDP_create_subset_L(lat->qlat, staggered_func, (void*)&nd,sizeof(nd),n);
    //printf("... done.\n");
  }
  return lat->staggered;
}

// This is a private function which defines even/odd subsections
// in different directions. Return 0 if an even coordinate, 1 if odd.
// nd is a special case: we want coordinates all even or all odd.
static int
eodir_func(QDP_Lattice *lat, int x[], void *args)
{
  int dir = ((int*)args)[0];
  int nd = ((int*)args)[1];
  if (dir != nd) {
    return x[dir]&1;
  } else {
    int k=0;
    for(int i=0; i<nd; i++) k += x[i]&1;
    if(k==0) return 0;
    else if(k==nd) return 1;
    return 2;
  }
}

QDP_Subset *
qhmcqdp_get_eodir(lattice_t *lat, int dir)
{
  if(lat->eodir[dir]==NULL) {
    int nd = lat->nd;
    int args[2] = {dir, nd};
    // We get only even or odd for 0->nd-1.
    // We get even, odd, or neither for nd.
    lat->eodir[dir] = QDP_create_subset_L(lat->qlat, eodir_func, (void*)args,
					  sizeof(args), (dir==nd?3:2));
  }
  return lat->eodir[dir];
}

QDP_Subset *
qhmc_qopqdp_qsubset_from_string(lua_State *L, lattice_t *lat,
				const char *s, int *n)
{
  *n = 0;
  if(s==NULL) return NULL;
  QDP_Subset *subs = NULL;
  QDP_Lattice *qlat = lat->qlat;
  switch(s[0]) {
  case 'a':
    subs = QDP_all_and_empty_L(qlat);
    if(strcmp(s,"all")==0) *n = 1;
    else *n = 2;
    break;
  case 'e':
    *n = 1;
    if(strncmp(s,"evenodd",7)==0 || strncmp(s,"evenandodd",10)==0) *n = 2;
    subs = QDP_even_and_odd_L(qlat);
    {
      int t=0;
      int nn = sscanf(s,"%*[^0-9]%i",&t);
      if(nn && t>=1 && t<=lat->nd) {
	subs = qhmcqdp_get_eodir(lat, t-1);
      }
    }
    break;
  case 'o':
    *n = 1;
    subs = 1 + QDP_even_and_odd_L(qlat);
    {
      int t=0;
      int nn = sscanf(s,"%*[^0-9]%i",&t);
      if(nn && t>=1 && t<=lat->nd) {
	subs = 1 + qhmcqdp_get_eodir(lat, t-1);
      }
    }
    break;
  case 's':
    if(strncmp(s,"staggered",9)==0) {
      int ns = 1 << QDP_ndim_L(qlat);
      if(strcmp(s+9,"")==0) {
	subs = qhmcqdp_get_staggered(lat);
	*n = ns;
      } else {
	int t=0;
	int nn = sscanf(s+9,"%i",&t);
	if(nn && t>=0 && t<ns) {
	  subs = &qhmcqdp_get_staggered(lat)[t];
	  *n = 1;
	}
      }
    }
    break;
  case 't':
    if(strncmp(s,"timeslice",9)==0) {
      int nt = QDP_coord_size_L(qlat,QDP_ndim_L(qlat)-1);
      if(strcmp(s+9,"s")==0) {
	subs = qhmcqdp_get_timeslices(lat);
	*n = nt;
      } else {
	int t=0;
	int nn = sscanf(s+9,"%i",&t);
	if(nn && t>=0 && t<nt) {
	  subs = &qhmcqdp_get_timeslices(lat)[t];
	  *n = 1;
	}
      }
    }
    break;
  }
  return subs;
}

static int
opt_qsubset_real(lua_State *L, int idx, lattice_t *lat, QDP_Subset *t)
{
  int len = 0;
  if(lua_isstring(L, idx)) {
    lua_pushvalue(L, idx);
    const char *s = lua_tostring(L, -1);
    lua_pop(L, 1);
    int n;
    QDP_Subset *qs = qhmc_qopqdp_qsubset_from_string(L, lat, s, &n);
    if(n) {
      len += n;
      if(t) for(int i=0; i<n; i++) t[i] = qs[i];
    }
  } else {
    int tidx = idx;
    subset_t *ls = qhmc_qopqdp_opt_subset(L, &tidx, 0, NULL);
    if(ls) {
      len++;
      if(t) t[0] = ls->qsub;
    }
  }
  return len;
}

static int
opt_as_qsubset_array_real(lua_State *L, int idx, lattice_t *lat, QDP_Subset *t)
{
  int len = 0;
  if(lua_istable(L, idx)) {
    int n = tableLength(L, idx);
    for(int i=1; i<=n; i++) {
      tableGetIndex(L, idx, i);
      int nn = opt_qsubset_real(L, -1, lat, t?t+len:t);
      lua_pop(L, 1);
      if(nn==0) {
	if(i==1) break;
	qlerror0(L, 1, "unknown subset at array index %i\n", i);
      }
      len += nn;
    }
  } else {
    len = opt_qsubset_real(L, idx, lat, t);
    if(len==1) len = -1; // use -1 for non-array return
  }
  return len;
}

// we allow one of:
// string
// subset
// array of mixed strings and subsets
int
qopqdp_opt_as_qsubset_array_len(lua_State *L, int idx, int required,
				lattice_t *lat, int def)
{
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  int len = opt_as_qsubset_array_real(L, idx, lat, NULL);
  if(len==0) {
    if(required) {
      qlerror0(L, 1, "invalid subset array at argument %i\n", idx);
    } else {
      len = def;
    }
  }
  return len;
}

void
qopqdp_opt_as_qsubset_array(lua_State *L, int *idx, int required,
			    lattice_t *lat, int n, QDP_Subset *t,
			    int dn, QDP_Subset *def)
{
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  int len = opt_as_qsubset_array_real(L, *idx, lat, t);
  if(len==0) {
    if(required) {
      qlerror0(L, 1, "invalid subset array at argument %i\n", *idx);
    } else {
      for(int i=0; i<abs(dn); i++) t[i] = def[i];
    }
  } else {
    (*idx)++;
  }
}

QDP_Subset
qopqdp_opt_qsubset(lua_State *L, int *idx, int required,
		   lattice_t *lat, QDP_Subset def)
{
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  int len = opt_qsubset_real(L, *idx, lat, NULL);
  QDP_Subset qsub = def;
  if(len==0) {
    if(required) {
      qlerror0(L, 1, "invalid subset array at argument %i\n", *idx);
    }
  } else {
    opt_qsubset_real(L, *idx, lat, &qsub);
    (*idx)++;
  }
  return qsub;
}

static struct luaL_Reg subset_reg[] = {
  { "__gc",           qopqdp_subset_gc },
  { NULL, NULL}
};

subset_t *
qopqdp_subset_create(lua_State *L, QDP_Subset qsub, subsetGroup_t *group)
{
  subset_t *s = lua_newuserdata(L, sizeof(subset_t));
  s->qsub = qsub;
  s->group = group;
  if(group) group->refcount++;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, subset_reg);
  }
  lua_setmetatable(L, -2);
  return s;
}
