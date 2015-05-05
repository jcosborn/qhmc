-- This example is part of the regression test suite.  Functions starting
-- with TEST are part of the test framework and can be ignored.

require 'topo'

nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)

TESTON()

seed = 987654321
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
g:unit()
getplaq(g)

g:random()
getplaq(g)

local nrep = 10
local nhb = 1
local nor = 1
local beta = 10
local coeffs = {plaq=1}
g:heatbath(nrep, nhb, nor, beta, coeffs)
getplaq(g)

t0 = clock()
ses,set,sq = symmEQ(g)
t0 = clock() - t0
printf("#time: %.8g\n", t0)
printf("se: %.8g\n", ses+set)
printf("sq: %.8g\n", sq)

t0 = clock()
ses,set,sq = symmEQ(g,1)
t0 = clock() - t0
printf("#time: %.8g\n", t0)
printf("se: %.8g\n", ses+set)
printf("sq: %.8g\n", sq)

t0 = clock()
ses,set,sq = symmEQ(g,1,"timeslices")
t0 = clock() - t0
printf("#time: %.8g\n", t0)
for i=1,#ses do
  printf("%i  se: %.7g  sq: %.7g\n", i, ses[i]+set[i], sq[i])
end

TESTOFF()
