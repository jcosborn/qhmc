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


t0 = clock()
se,sq = symmEQ(g)
t0 = clock() - t0
printf("time: %.10g\n", t0)
printf("se: %.10g\n", se)
printf("sq: %.10g\n", sq)
pe = plaqE(g)
printf("pe: %.10g\n", pe)

t0 = clock()
se,sq = symmEQ(g, 1)
t0 = clock() - t0
printf("time: %.10g\n", t0)
printf("se: %.10g\n", se)
printf("sq: %.10g\n", sq)

t0 = clock()
se,sq = symmEQ(g, 1, "timeslices")
t0 = clock() - t0
printf("time: %.10g\n", t0)
for i=1,#se do
  printf("%i  se: %.10g  sq: %.10g\n", i, se[i], sq[i])
end

eps = 0.01
for i=1,100 do
  wflow(g, {plaq=1}, eps, 1)
  se,sq = symmEQ(g, 1)
  t = eps * i
  printf("t: %-5g  se: %-15.10g  se*t^2: %-15.10g  sq: %-15.10g\n",
	 t, se, t*t*se, sq)
end
