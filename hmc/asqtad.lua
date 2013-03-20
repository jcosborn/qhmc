package.path = arg[0]:gsub("[^/]*.lua","?.lua") .. ";./hmc/?.lua;" .. package.path
require 'common'
require 'gaugeact'
require 'asqtadact'
require 'hmc'

trace(doTrace)

local nx = nx or 4
local nt = nt or 4
local beta = beta or 6
local u0 = u0 or 1
local nf = nf or 4
local mass = mass or 0.1
local ntraj = ntraj or 10
local tau = tau or 1
local nsteps = nsteps or 120
seed = seed or os.time()
local first = first or -1
local ffreq = ffreq or 3
local prec = prec or 1

local latsize = { nx, nx, nx, nt }
local vol = 1
printf("latsize =")
for k,v in ipairs(latsize) do vol=vol*v; printf(" %i",v); end
printf("\nvolume = %i\n", vol)
printf("beta = %g\n", beta)
printf("u0 = %g\n", u0)
printf("nf = %g\n", nf)
printf("mass = %g\n", mass)
printf("ntraj = %g\n", ntraj)
printf("tau = %g\n", tau)
printf("nsteps = %g\n", nsteps)
printf("seed = %i\n", seed)

p = {}
p.latsize = latsize
p.vol = vol
p.beta = beta
p.seed = seed
p.gaugeact = {type="plaquette", u0=u0}
--p.gaugeact = {type="symanzik_1loop", u0=u0}

function setpseudo(rhmc, mf, mb)
  local s1 = 4*mf*mf
  if mb then -- term (A+4mb^2)/(A+4mf^2)
    local s2 = 4*mb*mb
    local sr = math.sqrt(s1*s2)
    local c1 = s1 + s2 - 2*sr
    rhmc[#rhmc+1] = {
      --[[GR = {
	{ {1}, {sr-s2, s2} },
	{ {math.sqrt(c1), s2}, allfaceven=0, allfacodd=1 }
      },--]]
      GR = {
	{ {2, s2}, allfaceven=2*mf, allfacodd=1, allmass2=mb }
      },
      FA = {
	{ {1} },
	{ {math.sqrt(s2-s1), s1}, allfaceven=2*mf, allfacodd=1 }
      },
      MD = { {1}, {s2-s1, s1} }
    }
  else -- term 1/(A+4mf^2)
    rhmc[#rhmc+1] = {
      GR = {{ {1}, allfaceven=2*mf, allfacodd=1 }},
      FA = {{ {1,s1}, allfaceven=2*mf, allfacodd=1 }},
      MD = { {1,s1} }
    }
  end
end

local rhmc = {}
setpseudo(rhmc, mass)
setpseudo(rhmc, mass2, mass)

local smear = {}
--smear[#smear+1] = { type="fat7", coeffs={one_link=1} }
--smear[#smear+1] = { type="fat7", coeffs={one_link=0.5} }
--smear[#smear+1] = { type="fat7", coeffs={one_link=2} }
--smear[#smear+1] = { type="fat7", coeffs={three_staple=0.2} }
--smear[#smear+1] = { type="fat7", coeffs={one_link=0.4,three_staple=0.1} }
--smear[#smear+1] = { type="fat7", coeffs={one_link=0.4,three_staple=0.1,five_staple=0.2} }
--smear[#smear+1] = { type="fat7", coeffs={one_link=0.4,three_staple=0.1,five_staple=0.1,seven_staple=0.1} }
--smear[#smear+1] = { type="stout", rho=0.01 }
--smear[#smear+1] = { type="stout", rho=0.14 }

function projRat(s)
  local o = s.coeffs.one_link
  s.coeffs.one_link = nil
  return { type="projectRat", rho=1/o, smear=s }
  --return { type="projectRat", rho=4, smear=s }
end

scoeffs = asqtad_coeffs(u0)
scoeffs.three_staple = -scoeffs.three_staple
scoeffs.seven_staple = -scoeffs.seven_staple
scoeffs.one_link = scoeffs.one_link + 3*scoeffs.naik
scoeffs.naik = nil
scoeffs.one_link = scoeffs.one_link + 6*scoeffs.lepage
scoeffs.lepage = nil
--for k,v in pairs(scoeffs) do scoeffs[k] = -v end
--scoeffs={one_link=0.9,three_staple=0.1}
--scoeffs={one_link=0.9,five_staple=0.1}
--smear[#smear+1] = { type="fat7", coeffs=scoeffs }
--smear[#smear+1] = { type="projectU" }
--smear[#smear] = projRat(smear[#smear])
smear[#smear+1] = { type="hyp", alpha={0.5,0.5,0.4} }
myprint("smear = ", smear, "\n")

--coeffs = asqtad_coeffs(u0)
--coeffs.one_link = coeffs.one_link + 3*coeffs.naik
--coeffs.naik = nil
--coeffs.one_link = coeffs.one_link - 6*coeffs.lepage
--coeffs.lepage = 2*coeffs.lepage
--coeffs.naik = 0
coeffs = { one_link=1 }
--coeffs={one_link=0.4,three_staple=0.1,five_staple=0.1,seven_staple=0.1}
--coeffs={one_link=0.9,three_staple=0.1}
--coeffs={one_link=0.9,five_staple=0.1}

local act = {}
act.g = gaugeact(p)
act.f = asqtadact(act.g, {smear=smear,coeffs=coeffs,rhmc=rhmc})

do
  local mstring = string.format("%g", mass)
  mstring = mstring:sub(3)
  basefn = string.format("l%i%ib%03.0ff%im%s.", nx, nt, 100*beta, nf, mstring)
end
G = act.g:gaugeNew()
last = ntraj * tau
if first >= 0 then
  infn = basefn .. tostring(first)
  G:load(infn)
  last = last + first
else
  G:unit()
end
outfn = basefn .. tostring(last)

do
  local devavg,devmax = G.g:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  G.g:makeSU()
  devavg,devmax = G.g:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end

local fields = {}
function fields.save(f)
  f.GSave:set(f.G)
  f.F:random()
  act.f:refresh(f.G)
end
function fields.action(f)
  local Sgq = act.g:action(f.G)
  local Sgp = 0.5 * f.F:norm2() - 16*vol
  local Sfq = act.f:action(f.G)
  local S = Sgq + Sgp + Sfq
  printf("Sgq: %-8.6g  Sgp: %-8.6g  Sfq: %-8.6g\n", Sgq, Sgp, Sfq)
  return S
end
function fields.accept(f)
  local devavg,devmax = f.G.g:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  f.G.g:makeSU()
  devavg,devmax = f.G.g:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end
function fields.reject(f)
  f.G:set(f.GSave)
end
function fields.nfields(f)
  return 1
end
function fields.nforces(f)
  return { 2 }
end
--local gt = 0
function fields.updateField(f, i, eps)
  --gt = gt + eps
  --printf("Gupdate %g %g %g\n", eps, gt, 0.5*f.F:norm2()-16*vol)
  f.G:update(f.F, eps)
  --local ss,st = f.G:plaq()
  --printf("plaq %g %g\n", 3*ss, 3*st);
end
function fields.updateMomentum(f, i, tj, teps)
  for k,j in ipairs(tj) do
    if(j==1) then
      --printf("begin GFupdate %g %g\n", eps, 0.5*f.F:norm2()-16*f.vol)
      act.g:updateMomentum(f.F, f.G, teps[k])
      --printf("end   GFupdate %g %g\n", eps, 0.5*f.F:norm2()-16*f.vol)
      table.remove(tj, k)
      table.remove(teps, k)
    end
  end
  if #tj > 0 then
    for k=1,#tj do tj[k] = tj[k] - 1 end
    --myprint("forces ", tj, "\n")
    --printf("begin FFupdate[%i] %g %g\n", j-1, eps, 0.5*f.F:norm2()-16*f.vol)
    act.f:updateMomentum(f.F, f.G, teps, tj)
    --printf("end   FFupdate[%i] %g %g\n", j-1, eps, 0.5*f.F:norm2()-16*f.vol)
  end
end

fields.G = G
fields.GSave = act.g:gaugeNew()
fields.F = act.g:forceNew()

hmcparams = {}
hmcparams.tau = tau
local fp = {}
hmcparams.forceparams = fp
fp[1] = {}
fp[1][1] = {nsteps=ffreq*nsteps, intalg={type="omelyan", lambda=0.2}}
local rhmc1 = {}
for j=1,#rhmc do
  --printf("j = %i\n", j)
  --printf("ns = %i\n", nfsteps[j])
  fp[1][j+1] = {nsteps=nsteps, intalg={type="omelyan", lambda=0.2}}
  rhmc1[j] = {GR={},FA={},MD={}}
  rhmc1[j].GR[1] = {}
  rhmc1[j].GR[1].resid = 1e-6
  rhmc1[j].FA.resid = 1e-6
  rhmc1[j].MD.resid = 1e-5
  rhmc1[j].GR[1].solveopts = {
    prec = prec,
    restart = 500
  }
  rhmc1[j].FA.solveopts = {
    prec = prec,
    restart = 500
  }
  rhmc1[j].MD.solveopts = {
    prec = prec,
    restart = 500
  }
  rhmc1[j].MD.ffprec = prec
end
--myprint("rhmc1 = ", rhmc1, "\n")
copyto(rhmc, rhmc1)
--hmcparams.traceS = true
--qcd.defaults{qdpProfcontrol=1}

function measure(G)
  local t0 = clock()
  local ps,pt = G:plaq()
  printf("plaq ss: %-8g  st: %-8g  tot: %g\n", ps, pt, 0.5*(ps+pt))

  local nd = #G.a.latsize
  for i=1,nd do
    for j=i+1,nd do
      local lr,li = G.g:loop({i,j,-i,-j})
      printf("plaq%i%i:  %g\t%g\n", i, j, lr, li)
    end
  end

  local plpath = rep(-nd, G.a.latsize[nd])
  plpr,plpi = G.g:loop(plpath)
  printf("ploop:  %g\t%g\n", plpr, plpi)

  t0 = clock() - t0
  printf("meas time: %g\n", t0)

  io.stdout:flush()
end

totaltime = clock()
measure(G)
for traj=1,ntraj do
  local t0 = clock()
  hmcstep(fields, hmcparams)
  --act.g:checkSU(fields.G)
  t0 = clock() - t0

  printf("traj %i secs: %g\n", traj, t0)
  act.g:updateStats()
  act.f:updateStats()
  printf("GF     secs: %8.3f %3.0f%%  calls: %4.0f  mflops: %5.0f  rms: %6.4f  max: %6.4f\n", act.g.GFtime, 100*act.g.GFtime/t0, act.g.GFn, act.g.GFmflops, act.g.GFrms, act.g.GFmax)
  printf("GU     secs: %8.3f %3.0f%%  calls: %4.0f\n", act.g.GUtime, 100*act.g.GUtime/t0, act.g.GUn)
  printf("LL     secs: %8.3f %3.0f%%  calls: %4.0f  mflops: %5.0f\n", act.f.LLtime, 100*act.f.LLtime/t0, act.f.LLn, act.f.LLmflops)

  local fft,ffn,fff = 0,0,0
  for i=1,act.f.nff do
    printf("FF[%02i] secs: %8.3f %3.0f%%  calls: %4.0f  mflops: %5.0f  rms: %6.4f  max: %6.4f\n", i, act.f.FFtime[i], 100*act.f.FFtime[i]/t0, act.f.FFn[i], act.f.FFmflops[i], act.f.FFrms[i], act.f.FFmax[i])
    fft = fft + act.f.FFtime[i]
    ffn = ffn + act.f.FFn[i]
    fff = fff + act.f.FFflops[i]
  end
  printf("FFtot  secs: %8.3f %3.0f%%  calls: %4.0f  mflops: %5.0f\n", fft, 100*fft/t0, ffn, 1e-6*fff/fft)

  local cgt,cgn,cgf,cgi,cgm = 0,0,0,0,0
  for i=1,act.f.ncg do
    local ai = 0
    if act.f.CGn[i] > 0 then ai = act.f.CGits[i]/act.f.CGn[i] end
    printf("CG[%02i] secs: %8.3f %3.0f%%  calls: %4.0f  mflops: %5.0f  avgits: %5.0f  max: %5.0f\n", i, act.f.CGtime[i], 100*act.f.CGtime[i]/t0, act.f.CGn[i], act.f.CGmflops[i], ai, act.f.CGmaxits[i])
    cgt = cgt + act.f.CGtime[i]
    cgn = cgn + act.f.CGn[i]
    cgf = cgf + act.f.CGflops[i]
    cgi = cgi + act.f.CGits[i]
    cgm = math.max(cgm, act.f.CGmaxits[i])
  end
  printf("CGtot  secs: %8.3f %3.0f%%  calls: %4.0f  mflops: %5.0f  avgits: %5.0f  max: %5.0f\n", cgt, 100*cgt/t0, cgn, 1e-6*cgf/cgt, cgi/cgn, cgm)
  local ot = t0 - act.g.GFtime - act.g.GUtime - act.f.LLtime - fft - cgt
  printf("other  secs: %8.3f %3.0f%%\n", ot, 100*ot/t0)

  act.g:clearStats()
  act.f:clearStats()

  measure(fields.G)
end

if outfn then G:save(outfn) end

totaltime = clock() - totaltime
printf("total time: %g seconds\n", totaltime)
