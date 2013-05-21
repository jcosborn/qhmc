package.path = arg[0]:gsub("[^/]*.lua","?.lua") .. ";./hmc/?.lua;" .. package.path
require 'common'
require 'run'

trace(doTrace)

local nx = nx or 4
local nt = nt or 8
local seed = seed or os.time()
latsize = { nx, nx, nx, nt }
local nd = #latsize
qopqdp.lattice(latsize)
qopqdp.profile(0)
qopqdp.verbosity(0)
qopqdp.seed(seed)

vol = 1
printf("latsize =")
for k,v in ipairs(latsize) do vol=vol*v; printf(" %i",v); end
printf("\n")
printf("seed = %i\n", seed)

function plaq(g)
  local ss,st = g:action({plaq=1})
  local s = vol*qopqdp.Nc
  return ss/s, st/s
end

U = qopqdp.gauge()
if inlat then
  U:load(inlat)
else
  U:unit()
end

function rtloop(U, r, t)
  local lr,li = 0,0
  for mu=1,nd-1 do
    local p = {}
    for i=1,t do p[#p+1] = nd end
    for i=1,r do p[#p+1] = mu end
    for i=1,t do p[#p+1] = -nd end
    for i=1,r do p[#p+1] = -mu end
    local tr,ti = U:loop(p)
    lr = lr + tr
    li = li + ti
  end
  return lr/(nd-1),li/(nd-1)
end

function stats(U)
  local ps,pt = plaq(U)
  printf(" plaq ss: %g  st: %g  tot: %g\n", ps, pt, 0.5*(ps+pt))

  for i=1,nd do
    for j=i+1,nd do
      local lr,li = U:loop({i,j,-i,-j})
      printf(" plaq%i%i:  %g\t%g\n", i, j, lr, li)
    end
  end

  local plpath = rep(-nd, latsize[nd])
  plpr,plpi = U:loop(plpath)
  printf(" ploop:  %g\t%g\n", plpr, plpi)

  --[[
  for r=1,nx/2 do
    local wr0,wi0 = rtloop(U,r,0)
    local wr1,wi1 = rtloop(U,r,1)
    local wl = wr1/wr0
    printf(" wl%i%i: %12g   %g\n", r, 1, wl, -math.log(wl))
  end
  --]]
end


printf("plain\n")
stats(U)

function ape4d(alpha)
  local c = {}
  for mu = 1,4 do
    c[mu] = { [0] = 1-alpha }
    for nu = 1,4 do
      if nu ~= mu then
	c[mu][nu] = alpha/6
	c[mu][-nu] = alpha/6
      end
    end
  end
  return c
end
function ape3d(alpha)
  local c = {}
  for mu = 1,3 do
    c[mu] = { [0] = 1-alpha }
    for nu = 1,3 do
      if nu ~= mu then
	c[mu][nu] = alpha/4
	c[mu][-nu] = alpha/4
      end
    end
  end
  c[4] = { [0] = 1 }
  return c
end

--smear = { type="fat7", coeffs={one_link=1} }
--smear = { type="fat7", coeffs={one_link=0.5} }
--smear[#smear+1] = { type="fat7", coeffs={one_link=2} }
--smear[#smear+1] = { type="fat7", coeffs={three_staple=0.01} }
--smear = { type="fat7", coeffs={one_link=0.9,three_staple=0.1/6} }
--smear = { type="staples", coeffs=ape4d(0.1) }
--smear[#smear+1] = { type="stout", rho=0.01 }
smear = { type="stout", rho=0.14 }
--smear[#smear+1] = { type="stout", rho=0.14 }
--smear = { type="hyp", alpha={0.7,0.5,0.3} }

g = {U}
s = U
for i=1,5 do
  s = smearGauge({g=s}, {smear})
  smear.sg = nil
  printf("smeared %i\n", i)
  stats(s)
end
