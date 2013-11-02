#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.lattice";

lattice_t *
qopqdp_opt_lattice(lua_State *L, int *idx, int required, lattice_t *def)
{
  lattice_t *lat;
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

static void
qopqdp_lattice_free(lua_State *L, int idx)
{
  lattice_t *lat = qopqdp_check_lattice(L, idx);
  if(lat->timeslices) {
    QDP_destroy_subset(lat->timeslices);
    lat->timeslices = NULL;
  }
  if(lat->qlat!=QDP_get_default_lattice()) {
    QDP_destroy_lattice(lat->qlat);
    lat->qlat = NULL;
  }
}

static int
qopqdp_lattice_gc(lua_State *L)
{
  qopqdp_lattice_free(L, -1);
  return 0;
}

static struct luaL_Reg lattice_reg[] = {
  { "__gc",    qopqdp_lattice_gc },
  { NULL, NULL}
};

lattice_t *
qopqdp_lattice_create(lua_State *L, int nd, int size[])
{
  lattice_t *lat = lua_newuserdata(L, sizeof(lattice_t));
  lat->qlat = QDP_create_lattice(NULL, NULL, nd, size);
  lat->nd = nd;
  lat->timeslices = NULL;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, lattice_reg);
  }
  lua_setmetatable(L, -2);
  return lat;
}
