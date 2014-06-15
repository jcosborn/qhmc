require 'common'
require 'gaugeact'
require 'wilsonObservables'

--trace(true)
latsize = { 4, 4, 4, 8 }
aniso = 1/2.38
--aniso = 1
mass = -0.4125
--rsmear = 0.6
--nsmear = 30
rsmear = 0
nsmear = 0
prec = 1
restart = 500
resid = 1e-12

L = qopqdp.lattice(latsize)
qopqdp.profile(profile or 0)
qopqdp.verbosity(0)

--wr = qopqdp.writer("test.out", "<metadata/>")

function getplaq(g)
  local ps,pt = g:action{plaq=1}
  local lat = qopqdp.lattice()
  local nd,vol = #lat,1
  for i=1,nd do vol=vol*lat[i] end
  local s = 0.25*nd*(nd-1)*vol
  printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", ps/s, pt/s, 0.5*(ps+pt)/s)
end

g = qopqdp.gauge()
if fn then g:load(fn)
else
  seed = seed or os.time()
  qopqdp.seed(seed)
  --g:unit()
  g:random()
end

getplaq(g)
-- coulomb(j_decay, error, max iterations, overrelaxation param)
-- note that here 0->x, ..., 3->t
--g:coulomb(3, 1e-7, 1000, 1.2)
--getplaq(g)

w = qopqdp.wilson()
w:printcoeffs()
w:set(g, {aniso=aniso}, prec)

opts = { prec=prec, restart=restart }
src = w:quark()
dest = {}
dest2 = {}
local Nc = qopqdp.Nc
local Ns = 4

local cvCS = {}
local cvSC = {}
for spin = 1,Ns do
  cvSC[spin] = {}
  dest2[spin] = {}
  for color = 1,Nc do
    if spin==1 then cvCS[color] = {} end
    cvCS[color][spin] = L:colorVector()
    cvSC[spin][color] = cvCS[color][spin]
  end
end

for spin = 0,Ns-1 do
  for color = 0,Nc-1 do
    local k = Ns*color + spin + 1
    src:zero()
    -- point({coord},color,spin,re,im)
    src:point({0,0,0,0},color,spin,1)
    printf("src norm2: %g\n", src:norm2())
    src:smearGauss(g, 4, rsmear, nsmear);
    printf("src norm2: %g\n", src:norm2())
    dest[k] = w:quark()

    t0 = qopqdp.dtime()
    w:solve(dest[k], src, mass, resid, "all", opts)
    dt = qopqdp.dtime() - t0
    mf = 1e-6 * w:flops() / dt
    printf("its: %g  secs: %g  Mflops: %g\n", w:its(), dt, mf)
    dest[k]:smearGauss(g, 4, rsmear, nsmear);
    printf("dest norm2: %.16g\n", dest[k]:norm2())
    --wr:write(dest[#dest], "<field metadata/>")
  end
  for color = 0,Nc-1 do
    local k = Ns*color + spin + 1
    dest[k]:splitSpin(cvCS[color+1])
  end
  for spin2 = 0,Ns-1 do
    local cm = L:colorMatrix()
    dest2[spin2+1][spin+1] = cm
    cm:combineColor(cvSC[spin2+1])
  end
end

--[[
vk2 = w:quark()
for spin = 0,Ns-1 do
  for spin2 = 0,Ns-1 do
    local s,s2 = 0,0
    for color = 0,Nc-1 do
      local k = Ns*color + spin + 1
      local k2 = Ns*color + spin2 + 1
      vk2:gamma(15, dest[k2])
      s = s + dest[k]:dot(vk2,"timeslice0");
    end
    for spin3 = 0,Ns-1 do
      local z = 1
      if spin3>=2 then z = -1 end
      s2 = s2 + z*dest2[spin+1][spin3+1]:dot(dest2[spin2+1][spin3+1],"timeslice0")
    end
    printf("%i\t%i\t%s\t%s\n", spin, spin2, s, s2)
  end
end
--]]

--wr:write(dest, "<field metadata/>")
--wr:prop(dest, "<field metadata/>")
--for i=1,#dest do
  --printf("%i norm2: %g\n", i, dest[i]:norm2())
  --dest[i]:zero()
--end

--rd,md = qopqdp.reader("test.out")
--printf("%s\n", md)
--md = rd:read(dest)
--printf("%s\n", md)
--for i=1,#dest do
  --printf("%i norm2: %g\n", i, dest[i]:norm2())
--end

--mesons = wilsonMesons(dest)
--mesons = wilsonMesons2(dest2)
mesons = wilsonMesons3(dest2)

--[[
for g = 0,#mesons do
  myprint("meson["..tostring(g).."] = ",mesons[g],"\n")
end
--]]

printf("-= mesons =-\n")
for t=1,latsize[4] do
  for g = 0,#mesons do
    printf("%i\t%i\t%i\t%i\t%g\n", 0, t-1, g, g, mesons[g][t])
  end
end

baryons = wilsonBaryons3(dest2)

printf("-= baryons =-\n")
for t=1,latsize[4] do
  printf("%i", t-1)
  for g = 0,#baryons do
    printf("\t%g\t%g", baryons[g][t].r, baryons[g][t].i)
  end
  printf("\n")
end
