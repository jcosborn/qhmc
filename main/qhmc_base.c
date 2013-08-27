static struct luaL_Reg qhmc_base[] = {
  { "complex",    qhmc_complex },
  { NULL, NULL}
};

void
open_qopqdp(lua_State* L)
{
  luaL_register(L, "qopqdp", qopqdp_reg);
  int jobnum = QMP_get_job_number();
  int numjobs = QMP_get_number_of_jobs();
  lua_pushinteger(L, jobnum);
  lua_setglobal(L, "jobnum");
  lua_pushinteger(L, numjobs);
  lua_setglobal(L, "numjobs");
  lua_getglobal(L, "qopqdp");
  lua_pushinteger(L, QLA_Nc);
  lua_setfield(L, -2, "Nc");
  lua_pop(L, 1);
  open_qopqdp_smear(L);
}
