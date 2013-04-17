function clock()
  --return os.clock()
  --return os.time()
  return qopqdp.dtime()
end

function globalRand()
  local r = qopqdp.random()
  return r
end

local oldprintf = printf
function printf(...)
  if(qopqdp.master()) then
    oldprintf(...)
  end
end

function exit(n)
  os.exit(n)
end

local actmt = {}
actmt.__index = actmt

local gaugemt = {}
gaugemt.__index = gaugemt

local forcemt = {}
forcemt.__index = forcemt

-- gauge action

local gaugecoeffs={}
function gaugecoeffs.plaquette(p)
  return { plaq=1, rect=0, pgm=0, adjplaq=0 }
end
function gaugecoeffs.plaquette_adjoint(p)
  return { plaq=1, rect=0, pgm=0, adjplaq=p.adjFac }
end
function gaugecoeffs.symanzik_tree(p)
  local u0 = p.u0 or 1
  local u2 = u0*u0
  local c = { plaq=1, pgm=0, adjplaq=0 }
  c.rect = -1/(20*u2)
  return c
end
function gaugecoeffs.iwasaki(p)
  local u0 = p.u0 or 1
  local u2 = u0*u0
  local c1 = -0.331
  local c = { pgm=0, adjplaq=0 }
  c.plaq = 1 - 8*c1
  c.rect = c1/u2
  return c
end
function gaugecoeffs.dbw2(p)
  local u0 = p.u0 or 1
  local u2 = u0*u0
  local c1 = -1.4067
  local c = { pgm=0, adjplaq=0 }
  c.plaq = 1 - 8*c1
  c.rect = c1/u2
  return c
end
function gaugecoeffs.symanzik_1loop(p)
  local u0 = p.u0 or 1
  local nf = p.nf or 0
  local u2 = u0*u0
  local lu0 = math.log(u0)
  local c = { plaq=1, adjplaq=0 }
  c.rect = -(1 - (0.6264-0.4742*nf)*lu0 ) / (20*u2)
  c.pgm = (0.0433-0.012*nf)*lu0 / u2
  return c
end
function gaugecoeffs.symanzik_1loop_hisq(p)
  local u0 = p.u0 or 1
  local nf = p.nf or 0
  local u2 = u0*u0
  local lu0 = math.log(u0)
  local c = { plaq=1, adjplaq=0 }
  c.rect = -(1 - (0.6264-1.1746*nf)*lu0 ) / (20*u2)
  c.pgm = (0.0433-0.0156*nf)*lu0 / u2
  return c
end

function gaugeact(p)
  local a = {}
  a.latsize = p.latsize
  a.vol = p.vol
  a.beta = p.beta
  a.params = p.gaugeact or { type="plaquette" }
  local gcfunc = gaugecoeffs[a.params.type]
  if not gcfunc then
    printf("unknown gauge action type %s\n", a.params.type)
    exit(1)
  end
  a.coeffs = gcfunc(a.params)
  qopqdp.lattice(a.latsize)
  qopqdp.profile(profile or 0)
  qopqdp.verbosity(0)
  qopqdp.seed(p.seed)
  --printf("gauge coeffs: ")
  --for k,v in pairs(a.coeffs) do printf(" %-5s = % g\n", k, v) end
  myprint("gauge coeffs: ", a.coeffs, "\n")
  a.act0 = a.vol*(6*a.coeffs.plaq + 12*a.coeffs.rect + 16*a.coeffs.pgm + 6*a.coeffs.adjplaq)
  a.gf = actmt.forceNew(a)
  actmt.clearStats(a)
  return setmetatable(a, actmt)
end

function actmt.clearStats(a)
  a.GFtime = 0
  a.GFflops = 0
  a.GFn = 0
  a.GFnorm2 = 0
  a.GFmax = 0
  a.GUtime = 0
  a.GUn = 0
end

function actmt.updateStats(a)
  a.GFmflops = 1e-6 * a.GFflops / a.GFtime
  a.GFrms = math.sqrt(a.GFnorm2/(4*a.GFn*a.vol))
end

function actmt.gaugeNew(a)
  local g = {}
  g.a = a
  g.g = qopqdp.gauge()
  g.nupdate = 0
  return setmetatable(g, gaugemt)
end

function actmt.forceNew(a)
  local f = {}
  f.a = a
  f.f = qopqdp.force()
  return setmetatable(f, forcemt)
end

function actmt.action(a, g)
  local ss,st,sm = g.g:action(a.coeffs)
  return a.beta*(a.act0-ss-st) + sm
end

function actmt.force(a, f, g)
  local t0 = clock()
  g.g:force(f.f, a.coeffs, a.beta)
  a.GFtime = a.GFtime + clock() - t0
  --a.GFtime = a.GFtime + f.f:time()
  a.GFflops = a.GFflops + f.f:flops()
  a.GFn = a.GFn + 1
end

function actmt.updateMomentum(a, f, g, eps)
  local gf = a.gf
  a:force(gf, g)
  local s = eps
  f.f:update(gf.f, s)
  local gf2 = s*s*gf.f:norm2()
  local gfi = s*gf.f:infnorm()
  if a.printforce then
    printf("gauge force: norm %g  inf %g\n", gf2, gfi)
  end
  a.GFnorm2 = a.GFnorm2 + gf2
  if a.GFmax < gfi then a.GFmax = gfi end
end

-- gauge methods

function gaugemt.nupdates(g)
  return g.nupdate
end

function gaugemt.unit(g)
  g.g:unit()
  g.nupdate = g.nupdate + 1
end

function gaugemt.set(g, g2)
  g.g:set(g2.g)
  g.nupdate = g.nupdate + 1
end

function gaugemt.load(g, fn)
  g.g:load(fn)
  g.nupdate = g.nupdate + 1
end

function gaugemt.save(g, fn)
  g.g:save(fn, "")
end

function gaugemt.checkSU(g)
  return g.g:checkSU()
end

function gaugemt.makeSU(g)
  g.g:makeSU()
  g.nupdate = g.nupdate + 1
end

function gaugemt.plaq(g)
  local ss,st = g.g:action({plaq=1})
  local nd = #qopqdp.lattice()
  local s = 0.25*nd*(nd-1)*g.a.vol
  return ss/s, st/s
end

function gaugemt.ploop(g)
  return 0
end

function gaugemt.update(g, f, eps)
  local t0 = clock()
  g.g:update(f.f, eps)
  g.a.GUtime = g.a.GUtime + clock() - t0
  g.a.GUn = g.a.GUn + 1
  g.nupdate = g.nupdate + 1
end

-- force methods

function forcemt.random(f)
  f.f:random()
end

function forcemt.norm2(f)
  return f.f:norm2()
end
