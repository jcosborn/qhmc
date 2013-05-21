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

g = qopqdp.gauge()
if fn then g:load(fn)
else g:unit() end

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
    dest[#dest+1] = w:quark()

    t0 = qopqdp.dtime()
    w:solve(dest[#dest], src, mass, resid, "all", opts)
    dt = qopqdp.dtime() - t0
    mf = 1e-6 * w:flops() / dt
    printf("its: %g  secs: %g  Mflops: %g\n", w:its(), dt, mf)
    --wr:write(dest[#dest], "<field metadata/>")
  end
end
wr:prop(dest, "<field metadata/>")
