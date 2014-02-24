require 'common'
require 'gaugeact'
require 'hmc'

nx = nx or 4
nt = nt or 8
beta = beta or 6
--seed = seed or 54321
seed = seed or os.time()
tau = tau or 1
nsteps = nsteps or 40
ntraj = ntraj or 1
first = first or -1

latsize = { nx, nx, nx, nt }
--latsize = { nx, nt }
vol = 1; for k,v in ipairs(latsize) do vol=vol*v end
p = {}
p.latsize = latsize
p.vol = vol
p.beta = beta
p.seed = seed
if gact then
  p.gaugeact = { type = gact }
end
act = {}
act.g = gaugeact(p)

myprint("latsize = ", latsize, "\n")
printf("volume = %i\n", vol)
printf("seed = %i\n", seed)
printf("beta = %g\n", beta)
printf("tau = %g\n", tau)
printf("nsteps = %g\n", nsteps)
printf("ntraj = %g\n", ntraj)

basefn = string.format("l4x%i%ib%03.0f.", nx, nt, 100*beta)
last = ntraj * tau
if first >= 0 then
  infn = basefn .. tostring(first)
  U = loadGauge(infn)
  last = last + first
else
  G = act.g:gaugeNew()
  G:unit()
end
outfn = basefn .. tostring(last)

local fields = {}
function fields.save(f)
  f.GSave:set(f.G)
  f.F:random()
end
function fields.action(f)
  local Sgq = act.g:action(f.G)
  local Sgp = 0.5 * f.F:norm2()
  local S = Sgq + Sgp
  printf("Sgq: %-8.6g  Sgp: %-8.6g\n", Sgq, Sgp)
  return S
end
function fields.accept(f)
  local devavg,devmax = f.G:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  f.G:makeSU()
  devavg,devmax = f.G:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end
function fields.reject(f)
  f.G:set(f.GSave)
end
function fields.nfields(f)
  return 1
end
function fields.nforces(f)
  return { 1 }
end
function fields.updateField(f, i, eps)
  --printf("F1: fnorm: %g\tact: %g\n", 0.5*f.F:norm2(), act.g:action(f.G))
  f.G:update(f.F, eps)
  --printf("F2: fnorm: %g\tact: %g\n", 0.5*f.F:norm2(), act.g:action(f.G))
end
function fields.updateMomentum(f, i, j, eps)
  --printf("P1: fnorm: %g\tact: %g\n", 0.5*f.F:norm2(), act.g:action(f.G))
  act.g:updateMomentum(f.F, f.G, eps[1])
  --printf("P2: fnorm: %g\tact: %g\n", 0.5*f.F:norm2(), act.g:action(f.G))
end
--act.g.printforce = true

fields.G = G
fields.GSave = act.g:gaugeNew()
fields.F = act.g:forceNew()

local fp = {}
fp.nsteps = nsteps
fp.intalg = { type = "leapfrog" }
--fp.intalg = { type = "omelyan", lambda = 0.21 }
--hmcparams.traceS = true
--qcd.defaults{qdpProfcontrol=1}
hmcparams = {}
hmcparams.tau = tau
hmcparams.forceparams = {{ fp }}

function measure(G)
  local t0 = clock()
  local ps,pt = G:plaq()
  printf("plaq ss: %g  st: %g  tot: %g\n", ps, pt, 0.5*(ps+pt))

  local nd = #G.a.latsize
  for i=1,nd do
    for j=i+1,nd do
      local l = G.g:loop({i,j,-i,-j})
      printf("plaq%i%i:  %g\t%g\n", i, j, l.r, l.i)
    end
  end

  local plpath = rep(-nd, G.a.latsize[nd])
  local plp = G.g:loop(plpath)
  printf("ploop:  %g\t%g\n", plp.r, plp.i)

  t0 = clock() - t0
  printf("meas time: %g\n", t0)

  io.stdout:flush()
end

totaltime = clock()
measure(G)
for traj=1,ntraj do
  t0 = clock()
  hmcstep(fields, hmcparams)
  --act.g:checkSU(fields.G)
  t0 = clock() - t0

  printf("traj %i secs: %g\n", traj, t0)
  act.g:updateStats()
  printf("GF secs: %8.3f  mflops: %5.0f\n", act.g.GFtime, act.g.GFmflops)
  printf("GU secs: %8.3f\n", act.g.GUtime)
  ot = t0 - act.g.GFtime - act.g.GUtime
  printf("?? secs: %8.3f\n", ot)
  act.g:clearStats()

  measure(fields.G)
end

--if outfn then saveGauge(G, outfn) end

totaltime = clock() - totaltime
printf("total time: %g seconds\n", totaltime)
