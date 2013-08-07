package.path = arg[0]:gsub("[^/]*.lua","?.lua") .. ";./hmc/?.lua;" .. package.path
require 'common'
require 'gaugeact'

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

U = qopqdp.gauge()
if inlat then
   U:load(inlat)
else
   U:unit()
end

function plaq(g)
  local ss,st = g:action({plaq=1})
  local s = vol*qopqdp.Nc
  return ss/s, st/s
end

function wflow(u, coeffs, eps, nsteps)
  -- u1 = exp((eps/4)f0) u0
  -- u2 = exp((eps*8/9)f1-(eps*17/36)f0) u1
  -- u3 = exp((eps*3/4)f2-(eps*8/9)f1+(eps*17/36)f0) u2
  local f = qopqdp.force()
  local ft = qopqdp.force()
  for i=1,nsteps do
    --[[
    u:force(f, coeffs)
    u:update(f, eps)
    --]]
    --[[
    local alpha = 0.5
    u:force(ft, coeffs)
    u:update(ft, eps*alpha)
    u:force(f, coeffs)
    f:update(ft, 2*alpha-1-2*alpha*alpha)
    u:update(f, eps*0.5/alpha)
    --]]
    u:force(f, coeffs)
    u:update(f, eps/4)
    u:force(ft, coeffs)
    f:update(ft, -32/17)
    u:update(f, eps*-17/36)
    u:force(ft, coeffs)
    f:update(ft, 27/17)
    u:update(f, eps*17/36)
  end
end

coeffs = { plaq=1 }
--coeffs = { plaq=1, rect=-0.05 }
eps = 0.025
tmax = nx*nx/16
tmax = 1
nsteps = tmax/eps

function stats(u, i)
  local ss,st = plaq(u)
  local t = eps*i
  local e = 6*(2-ss-st)
  printf("%g\t%g\n", t, e*t*t)
end

t0 = clock()

stats(U, 0)
for i=1,nsteps do
   wflow(U, coeffs, eps, 1)
  stats(U, i)
  if(i%10==0) then collectgarbage() end
end

t1 = clock()

printf("# time: %g seconds\n", t1-t0)
