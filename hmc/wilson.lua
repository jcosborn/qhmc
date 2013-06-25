package.path = arg[0]:gsub("[^/]*.lua","?.lua") .. ";./hmc/?.lua;" .. package.path
require 'common'
require 'run'
require 'mg'

trace(doTrace)
--noMG = true

local nx = nx or 4
local nt = nt or 8
local beta = beta or 4
local beta_a = beta_a or 0
local u0 = u0 or 1
local aniso = aniso or {}
local nf = nf or 2
local mass = mass or 0.0
local rho = rho or 0.0
local clov = clov or 0
local clov_s = clov_s or clov
local clov_t = clov_t or clov
_G.mass = mass

aniso.xi0 = aniso.xi0 or 1
aniso.nu = aniso.nu or 1
aniso.gmom = aniso.gmom or 1

local rhmc = {}
local hmcmasses = masses or {-0.2}
--local hmcmasses = { mass, 1.4*mass }
local seed = 1316844761

--local inlat = inlat or nil
local inlat = inlat or nil
--local inlat = "l84f8b40m04a.2700.scidac"
local outlat = outlat or nil
--local outlat = "f8x88b40m01.100"
local warmup = warmup or 0
local ntraj = ntraj or 10
local tau = tau or 1
local ngsteps = ngsteps or {24}
local nfsteps = nfsteps or {4, 8}
--local nfsteps = { 80, 80 }

nfsteps = repelem(nfsteps, nf/2)
local grcg = { prec=2, resid=1e-8, restart=500 }
local facg = { prec=2, resid=1e-8, restart=500 }
local mdcg = { prec=2, resid=1e-8, restart=500 }
local ffprec = 2
--local gintalg = {type="leapfrog"}
--local gintalg = {type="omelyan"}
--local gintalg = {type="omelyan", lambda=0.22}
local gintalg = {type="2MNV", lambda=0.1932}
--local gintalg = {type="omelyan", lambda=0.33}

--local fintalg = {type="leapfrog"}
--local fintalg = {type="omelyan", lambda=0.22}
local fintalg = {type="2MNV", lambda=0.1932}

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
--smear[#smear+1] = { type="stout", rho=rho}

--- end of parameters

--[[
--mf: fermionic mass
--mb: bosonic mass (Hasenbusch mass preconditioning, or Pauli-Villars field in DWF)
--]]--
function setpseudo(rhmc, mf, mb)
  if mb then
    rhmc[#rhmc+1] = {
      GR = {mf, mb},
      FA = {mf, mb},
      MD = {mf, mb}
    }
  else
    rhmc[#rhmc+1] = {
      GR = { mf },
      FA = { mf },
      MD = { mf }
    }
  end
end

for i=1,#hmcmasses do
  for j=1,nf/2 do
    if i < #hmcmasses then
      setpseudo(rhmc, hmcmasses[i], hmcmasses[i+1])
    else
      setpseudo(rhmc, hmcmasses[i])
    end
  end
end
local npseudo = #rhmc
-- p: lattice parameters, including fermion and gauge actions, smearing, anisotropy, and lattice dimensions, etc. 
local p = {}
p.latsize = { nx, nx, nx, nt }
latsize = p.latsize
p.seed = seed or os.time()
p.beta = beta
p.nf = nf
p.u0 = u0
p.xi0 = aniso.xi0
p.gmom_var = { 1, 1, 1, aniso.gmom }
--p.gaugeact = {type="symanzik_1loop_hisq", u0=p.u0, nf=p.nf}
p.gaugeact = {type="plaquette"}
p.npseudo = npseudo
p.fermact = {type="wilson", rhmc=rhmc}
p.fermact.smear = smear
p.fermact.coeffs = {clov_s=clov_s, clov_t=clov_t, aniso=aniso.nu/aniso.xi0}

local rhmc0 = copy(rhmc)
local acts = setupacts(p)
--myprint("rhmc0 = ", rhmc0, "\n")

-- r: run parameters
local r = {}
r.tau = tau

-- fp: force parameters
local fp = {}
r.forceparams = fp

-- Lua index starts from 1
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
--  acts:random()
end

if warmup then
  r.ntraj = warmup
  r.md = true
  acts:run(r)
  r.md = false
end

r.ntraj = ntraj - warmup -- warmup trajectories are counted towards the total.
acts:run(r)

if outlat then
  acts:save(outlat)
end
