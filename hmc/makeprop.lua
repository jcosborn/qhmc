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
