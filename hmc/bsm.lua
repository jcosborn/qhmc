package.path = "./hmc/?.lua;" .. package.path
require 'common'
require 'run'

local nx = 8
local nt = 4
local beta = 4
local u0 = 0.8163
local nf = 8
local mass = 0.04
local rhmc = {}
--local hmcmasses = { mass }
local hmcmasses = { mass, 20*mass }
--local hmcmasses = { mass, 2*mass, 3*mass, 4*mass, 5*mass, 6*mass, 7*mass, 8*mass }
local seed = 1316844761

--local inlat = nil
--local inlat = "f8x88b40m01.100"
local inlat = "l84f8b40m04a.2700.scidac"
local outlat = nil
--local outlat = "f8x88b40m01.100"

local ntraj = 10
local tau = 2
local ngsteps = 240
--local nfsteps = { 80 }
local nfsteps = { 80, 80 }
--local nfsteps = { 100, 100, 100, 100, 100, 100, 100, 200 }
nfsteps = repelem(nfsteps, nf/4)
local grcg = { prec=1, resid=1e-5, restart=500 }
local facg = { prec=1, resid=1e-5, restart=500 }
local mdcg = { prec=1, resid=1e-5, restart=500 }
local ffprec = 1
--local gintalg = {type="leapfrog"}
--local gintalg = {type="omelyan"}
--local gintalg = {type="omelyan", lambda=0.2}
local gintalg = {type="omelyan", lambda=0.33}
--local fintalg = {type="omelyan", lambda=0.2}
local fintalg = {type="omelyan", lambda=0.33}

local pbp = {}
pbp[1] = { reps=1 }
pbp[1].mass = mass
pbp[1].resid = 1e-6
pbp[1].opts = { restart=500, max_restarts=5, max_iter=2000 }

--- end of parameters

function setpseudo(rhmc, mf, mb)
  local s1 = 4*mf*mf
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
  else -- term 1/(A+4mf^2)
    rhmc[#rhmc+1] = {
      GR = {{ {1}, allfaceven=2*mf, allfacodd=1 }},
      FA = {{ {1,s1}, allfaceven=2*mf, allfacodd=1 }},
      MD = { {1,s1} }
    }
  end
end
for i=1,#hmcmasses do
  for j=1,nf/4 do
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
p.seed = seed or os.time()
p.beta = beta
p.nf = nf
p.u0 = u0
p.gaugeact = {type="symanzik_1loop_hisq", u0=p.u0, nf=p.nf}
p.npseudo = npseudo
p.fermact = {type="hisq", rhmc=rhmc}

local rhmc0 = copy(rhmc)
local acts = setupacts(p)
--myprint("rhmc0 = ", rhmc0, "\n")

local r = {}
r.ntraj = ntraj
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
  rhmc1[j].GR[1] = {}
  rhmc1[j].GR[1].resid = grcg.resid
  rhmc1[j].FA[1] = {}
  rhmc1[j].FA[1].resid = facg.resid
  rhmc1[j].MD.resid = mdcg.resid
  rhmc1[j].GR[1].solveopts = {
    prec = grcg.prec,
    restart = grcg.restart
  }
  rhmc1[j].FA[1].solveopts = {
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

acts:run(r)

if outlat then
  acts:save(outlat)
end
