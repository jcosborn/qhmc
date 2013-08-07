package.path = arg[0]:gsub("[^/]*.lua","?.lua") .. ";./hmc/?.lua;" .. package.path
require 'common'
require 'gaugeact'

latsize = { 4, 4, 4, 8 }
prec = 1
restart = 500
resid = 1e-6
mass = 0

qopqdp.lattice(latsize)
qopqdp.profile(profile or 0)
qopqdp.verbosity(0)

wr = qopqdp.writer("test.out", "<metadata/>")

function getplaq(g)
  local ps,pt = g:plaq()
  printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", ps, pt, 0.5*(ps+pt))
end

g = qopqdp.gauge()
if fn then g:load(fn)
else
  seed = seed or os.time()
  qopqdp.seed(seed)
  --g:unit()
  g:random()
end

w = qopqdp.wilson()
w:printcoeffs()
w:set(g, prec)

opts = { prec=prec, restart=restart }
src = w:quark()
dest = {}
local Nc = qopqdp.Nc
local Ns = 4

for color = 0,Nc-1 do
  for spin = 0,Ns-1 do
    src:zero()
    -- point({coord},color,spin,re,im)
    src:point({0,0,0,0},color,spin,1,0)
    printf("src norm2: %g\n", src:norm2())
    src:smearGauss(g, 0.9, 10);
    printf("src norm2: %g\n", src:norm2())
    dest[#dest+1] = w:quark()

    t0 = qopqdp.dtime()
    w:solve(dest[#dest], src, mass, resid, "all", opts)
    dt = qopqdp.dtime() - t0
    mf = 1e-6 * w:flops() / dt
    printf("its: %g  secs: %g  Mflops: %g\n", w:its(), dt, mf)
    --wr:write(dest[#dest], "<field metadata/>")
  end
end
--wr:write(dest, "<field metadata/>")
wr:prop(dest, "<field metadata/>")
for i=1,#dest do
  printf("%i norm2: %g\n", i, dest[i]:norm2())
  dest[i]:zero()
end

rd,md = qopqdp.reader("test.out")
printf("%s\n", md)
md = rd:read(dest)
printf("%s\n", md)

for i=1,#dest do
  printf("%i norm2: %g\n", i, dest[i]:norm2())
end

local gammaelem = {}
gammaelem[0] = { {1,0},{2,0},{3,0},{4,0} }
gammaelem[1] = { {4,1},{3,1},{2,3},{1,3} }
gammaelem[2] = { {4,2},{3,0},{2,0},{1,2} }
gammaelem[4] = { {3,1},{4,3},{1,3},{2,1} }
gammaelem[8] = { {3,0},{4,0},{1,0},{2,0} }
function gtimes(a, b)
  local c = {}
  for i=1,4 do
    c[i] = { b[a[i][1]][1], (a[i][2]+b[a[i][1]][2])%4 }
  end
  return c
end
for g=0,15 do
  if not gammaelem[g] then
    local m = gammaelem[0]
    local h,b = g,1
    repeat
      if h%2==1 then
	--printf("%i\t%i\t%i\n", g, h, b)
	m = gtimes(gammaelem[b], m)
      end
      h,b = math.floor(h/2),2*b
    until h == 0
    gammaelem[g] = m
  end
  --myprint("g[",tostring(g),"]=",gammaelem[g],"\n")
end

qt = w:quark()
function mydot(v, g, color, spin)
  local spin2,phase = table.unpack(gammaelem[g][spin])
  local k = (color-1)*Ns + spin
  local k2 = (color-1)*Ns + spin2
  local vk = v[k]
  if g ~= 0 then
    qt:gamma(g, vk)
    vk = qt
  end
  local t
  if phase==0 or phase==2 then
    if k==k2 and g==0 then
      t = v[k2]:norm2("timeslices")
    else
      t = v[k2]:Re_dot(vk,"timeslices")
    end
  else
    t = v[k2]:Im_dot(vk,"timeslices")
  end
  --printf("%i\t%i\t%g\n", k2, k, t[1])
  if phase==1 or phase==2 then
    for i=1,#t do t[i] = -t[i]; end
  end
  return t
end

local pions = {}
for g = 0,#gammaelem do
  pions[g] = {}
  for color = 1,Nc do
    for spin = 1,Ns do
      local t = mydot(dest, g, color, spin)
      --printf("%g\n", t[1])
      for j=1,#t do
	pions[g][j] = t[j] + (pions[g][j] or 0)
      end
    end
  end
end

for g = 0,#gammaelem do
  myprint("pion["..tostring(g).."] = ",pions[g],"\n")
end
