-- A sample lua script to measure symmE and symmQ
--  using O(a^2) or O(a^4) improved definitions.

require 'common'
require 'gaugeact'
require 'topo'

--fn="test_config_scidac"
--latsize = { 24, 24, 24, 48 }
latsize = { 4, 4, 4, 4 }
--latsize = { 4, 4, 4, 8 }
--latsize = { 8, 8, 8, 8 }

if fn then 
  latsize = qopqdp.getFileLattice(fn)
end

L = qopqdp.lattice(latsize)
qopqdp.profile(profile or 0)
qopqdp.verbosity(0)
Nc = qopqdp.Nc
Ns = 4

seed = seed or os.time()
qopqdp.seed(seed)

function getplaq(g)
  local ps,pt = g:action{plaq=1}
  local lat = qopqdp.lattice()
  local nd,vol = #lat,1
  for i=1,nd do vol=vol*lat[i] end
  local s = 0.25*nd*(nd-1)*vol
  printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", ps/s, pt/s, 0.5*(ps+pt)/s)
end

g = qopqdp.gauge()
if fn then 
  printf("CFG= %s\n",fn)
  g:load(fn)
else
  --g:unit()
  g:random()
  local nrep = 10
  local nhb = 1
  local nor = 1
  local beta = 10
  local coeffs = {plaq=1}
  g:heatbath(nrep, nhb, nor, beta, coeffs)
end

getplaq(g)

pe = plaqE(g)
printf("#pe: %.10g\n", pe)

t0 = clock()
ses,set,sq = symmEQ(g)
t0 = clock() - t0
printf("#time: %.10g\n", t0)
printf("WFLOW 0 %g %g %g\n", set, ses, sq)

--[[
t0 = clock()
ses,set,sq = symmEQ(g, 1)
t0 = clock() - t0
printf("time: %.10g\n", t0)
printf("set1: %.10g\n", set)
printf("ses1: %.10g\n", ses)
printf("sq1: %.10g\n", sq)
--]]

--[[
t0 = clock()
ses,set,sq = symmEQ(g, 1, "timeslices")
t0 = clock() - t0
printf("time: %.10g\n", t0)
for i=1,#ses do
  printf("%i  set: %.10g  ses: %.10g  sq: %.10g\n", i, set[i], ses[i], sq[i])
end
--]]

eps = eps or 0.02
tmax = tmax or 0
t2e = 0
i = 1
while true do
  wflow(g, {plaq=1}, eps, 1)
  ses,set,sq = symmEQ(g)
  t = eps * i
  printf("WFLOW %g %g %g %g\n", t, set, ses, sq)
  if tmax>0 then
    if t>(tmax-0.5*eps) then break end
  else
    t2eo = t2e
    t2e = t*t*(ses+set)
    t2ed = (t2e-t2eo)/eps
    --printf("t2e: %g  t2ed: %g\n", t2e, t2ed)
    if t>1 and i%20==0 then
      if t2e>0.45 and t2ed>0.35 then break end
    end
  end
  i = i + 1
end
