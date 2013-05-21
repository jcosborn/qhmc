package.path = arg[0]:gsub("[^/]*.lua","?.lua") .. ";./hmc/?.lua;" .. package.path
require 'common'
require 'run'
require 'mg'

trace(doTrace)
noMG = true

local nx = nx or 4
local nt = nt or 8
local beta = beta or 4
local beta_a = beta_a or 0
local u0 = u0 or 1
local nf = nf or 2
local mass = mass or 0.0
local rho = rho or 0.02
_G.mass = mass

local rhmc = {}
local hmcmasses = { mass }
--local hmcmasses = { mass, 20*mass }
local seed = 1316844761

--local inlat = inlat or nil
local inlat = inlat or nil
--local inlat = "l84f8b40m04a.2700.scidac"
local outlat = outlat or nil
--local outlat = "f8x88b40m01.100"
local warmup = warmup or nil
local ntraj = ntraj or 10
local tau = tau or 0.1
local ngsteps = ngsteps or 240
local nfsteps = nfsteps or { 80 }
--local nfsteps = { 80, 80 }

nfsteps = repelem(nfsteps, nf/2)
local grcg = { prec=2, resid=1e-10, restart=500 }
local facg = { prec=2, resid=1e-10, restart=500 }
local mdcg = { prec=2, resid=1e-8, restart=500 }
local ffprec = 2
--local gintalg = {type="leapfrog"}
--local gintalg = {type="omelyan"}
local gintalg = {type="omelyan", lambda=0.22}
--local gintalg = {type="omelyan", lambda=0.33}
--local fintalg = {type="leapfrog"}
local fintalg = {type="omelyan", lambda=0.22}

local pbp = {}
pbp[1] = { reps=1 }
pbp[1].mass = mass
pbp[1].resid = 1e-6
pbp[1].opts = { restart=500, max_restarts=5, max_iter=2000 }

local smear = {}
--smear[#smear+1] = { type="fat7", coeffs={one_link=1} }
--smear[#smear+1] = { type="fat7", coeffs={one_link=0.5} }
--smear[#smear+1] = { type="fat7", coeffs={one_link=2} }
--smear[#smear+1] = { type="fat7", coeffs={three_staple=0.2} }
--smear[#smear+1] = { type="fat7", coeffs={one_link=0.4,three_staple=0.1} }
--smear[#smear+1] = { type="stout", rho=0.01 }
--smear[#smear+1] = { type="stout", rho=rho}
smear[#smear+1] = { type="stout", rho=rho}

--- end of parameters

-- A = 
function setpseudo(rhmc, mf, mb)
  --kf = 0.5/(4+mf)
  if mb then -- term (A+4mb^2)/(A+4mf^2)
    local s2 = 4*mb*mb
    local sr = math.sqrt(s1*s2)
    local c1 = s1 + s2 - 2*sr
    rhmc[#rhmc+1] = {
      GR = {
	{ {1}, {sr-s2, s2} },
	{ {math.sqrt(c1), s2}, allfaceven=0, allfacodd=1 }
      },
      FA = {
	{ {1} },
	{ {math.sqrt(s2-s1), s1}, allfaceven=2*mf, allfacodd=1 }
      },
      MD = { {1}, {s2-s1, s1} }
    }
  else -- term 1/[(1-A)(1-A^+)]
    rhmc[#rhmc+1] = {
      GR = { mf },
      FA = { mf },
      MD = { mf }
    }
  end
end
for i=1,#hmcmasses do
  for j=1,nf/2 do
    if i<#hmcmasses then
      setpseudo(rhmc, hmcmasses[i], hmcmasses[i+1])
    else
      setpseudo(rhmc, hmcmasses[i])
    end
  end
end
local npseudo = #rhmc

local p = {}
p.latsize = { nx, nx, nx, nt }
latsize = p.latsize
p.seed = seed or os.time()
p.beta = beta
p.nf = nf
p.u0 = u0
--p.gaugeact = {type="symanzik_1loop_hisq", u0=p.u0, nf=p.nf}
p.gaugeact = {type="plaquette_adjoint",adjFac=beta_a}
p.npseudo = npseudo
p.fermact = {type="wilson", rhmc=rhmc}
p.fermact.smear = smear

local rhmc0 = copy(rhmc)
local acts = setupacts(p)
--myprint("rhmc0 = ", rhmc0, "\n")

local r = {}
r.tau = tau

local fp = {}
r.forceparams = fp
fp[1] = {}
fp[1][1] = {nsteps=ngsteps, intalg=gintalg}
local rhmc1 = {}
for j=1,npseudo do
  --printf("j = %i\n", j)
  --printf("ns = %i\n", nfsteps[j])
  fp[1][j+1] = {nsteps=nfsteps[j], intalg=fintalg}
  rhmc1[j] = {GR={},FA={},MD={}}
  rhmc1[j].GR.resid = grcg.resid
  rhmc1[j].FA.resid = facg.resid
  rhmc1[j].MD.resid = mdcg.resid
  rhmc1[j].GR.solveopts = {
    prec = grcg.prec,
    restart = grcg.restart
  }
  rhmc1[j].FA.solveopts = {
    prec = facg.prec,
    restart = facg.restart
  }
  rhmc1[j].MD.solveopts = {
    prec = mdcg.prec,
    restart = mdcg.restart
  }
  rhmc1[j].MD.ffprec = ffprec
end
--myprint("rhmc1 = ", rhmc1, "\n")
copyto(rhmc, rhmc1)

r.pbp = pbp
--myprint("runparams = ", r, "\n")

if inlat then
  acts:load(inlat)
  local ps,pt = acts.fields.G:plaq()
  printf("plaq ss: %g  st: %g  tot: %g\n", ps, pt, 0.5*(ps+pt))
else
  acts:unit()
end

if warmup then
  r.ntraj = warmup
  r.md = true
  acts:run(r)
  r.md = false
end

r.ntraj = ntraj
acts:run(r)

if outlat then
  acts:save(outlat)
end
