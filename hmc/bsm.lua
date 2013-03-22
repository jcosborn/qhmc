package.path = arg[0]:gsub("[^/]*.lua","?.lua") .. ";./hmc/?.lua;" .. package.path
require 'common'
require 'run'

--profile = 1

local nx = nx or 8
local nt = nt or 4
local beta = beta or 6
local beta_a = beta_a or 0
local u0 = u0 or 1
local nf = nf or 4
local mass = mass or 0.01
local mass2 = mass2 or 0.1
local prec = prec or 1
local faresid = faresid or 1e-5
local grresid = faresid
local mdresid = mdresid or 1e-5
local restart = 2000
local use_prev_soln = use_prev_soln or 0
--printf("faresid: %g\n", faresid)

local rhmc = {}
--local hmcmasses = { mass }
local hmcmasses = { mass, mass2 }
--local hmcmasses = { mass, 2*mass, 3*mass, 4*mass, 5*mass, 6*mass, 7*mass, 8*mass }
--local seed = 1316844761
local seed = seed or os.time()

local inlat = inlat or nil
--local inlat = "f8x88b40m01.100"
--local inlat = "l84f8b40m04a.2700.scidac"
local outlat = outlat or nil
--local outlat = "f8x88b40m01.100"

local ntraj = ntraj or 10
local tau = tau or 1
ngpfs = ngpfs or 5
--lambdaG = lambdaG or 0.2
--lambdaF = lambdaF or 0.2
local nsteps = nsteps or 80
local nsteps2 = nsteps2 or 80

local ngsteps = ngsteps or ngpfs*nsteps
--local nfsteps = nfsteps or { nsteps }
local nfsteps = nfsteps or { nsteps, nsteps2 }
--local nfsteps = { 80, 80 }
--local nfsteps = { 100, 100, 100, 100, 100, 100, 100, 200 }
nfsteps = repelem(nfsteps, nf/4)
local grcg = { prec=prec, resid=grresid, restart=restart }
local facg = { prec=prec, resid=faresid, restart=restart }
local mdcg = { prec=prec, resid=mdresid, restart=restart }
local ffprec = prec
--local gintalg = {type="leapfrog"}
--local gintalg = {type="omelyan"}
--local gintalg = {type="omelyan", lambda=0.2}
local gintalg = {type="omelyan", lambda=lambdaG}
--local gintalg = {type="2MNV", lambda=lambdaG}
--local fintalg = {type="omelyan", lambda=0.2}
local fintalg = {type="omelyan", lambda=lambdaF}
--local fintalg = {type="2MNV", lambda=lambdaF}

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
      --[[GR = {
	{ {1}, {sr-s2, s2} },
	{ {math.sqrt(c1), s2}, allfaceven=0, allfacodd=1 }
      },--]]
      GR = {
	{ {2, s2}, allfaceven=-2*mf, allfacodd=1, allmass2=-mb }
      },
      FA = {
        --{ {1,s2} },
	{ {1} },
	{ {math.sqrt(s2-s1), s1}, allfaceven=2*mf, allfacodd=1 }
      },
      --MD = { {1}, {s2-s1, s1} }
      MD = { {s2-s1, s1} }
    }
  else -- term 1/(A+4mf^2)
    rhmc[#rhmc+1] = {
      GR = {{ {1}, allfaceven=2*mf, allfacodd=-1 }},
      FA = {{ {1,s1}, allfaceven=2*mf, allfacodd=1 }},
      --FA = {{ {1,s1} }},
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
--p.gaugeact = {type="plaquette"}
p.gaugeact = {type="plaquette_adjoint", adjFac=beta_a}
--p.gaugeact = {type="symanzik_1loop_hisq", u0=p.u0, nf=p.nf}
p.npseudo = npseudo

function ape4d(alpha)
  local c = {}
  for mu = 1,4 do
    c[mu] = { [0] = 1-alpha }
    for nu = 1,4 do
      if nu ~= mu then
        c[mu][nu] = alpha/6
        c[mu][-nu] = alpha/6
      end
    end
  end
  return c
end
function ape3d(alpha)
  local c = {}
  for mu = 1,3 do
    c[mu] = { [0] = 1-alpha }
    for nu = 1,3 do
      if nu ~= mu then
        c[mu][nu] = alpha/4
        c[mu][-nu] = alpha/4
      end
    end
  end
  c[4] = { [0] = 1 }
  return c
end

local smear = {}
--smear[#smear+1] = { type="staples", coeffs={{[0]=1},{[0]=1},{[0]=1},{[0]=1}} }
--smear[#smear+1] = { type="staples", coeffs=ape4d(0.5) }
--smear[#smear+1] = { type="hyp", alpha={0.,0.,0.} }
--smear[#smear+1] = { type="hyp", alpha={0.0,0.5,0.5} }
smear[#smear+1] = { type="hyp", alpha={0.4,0.5,0.5} }
--smear[#smear+1] = { type="hyp", alpha={0.5,0.6,0.6} }
myprint("smear = ", smear, "\n")
coeffs = { one_link=1 }
p.fermact = {type="asqtad", smear=smear, coeffs=coeffs, rhmc=rhmc}

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
  rhmc1[j].FA.resid = facg.resid
  rhmc1[j].MD.resid = mdcg.resid
  rhmc1[j].GR[1].solveopts = {
    prec = grcg.prec,
    restart = grcg.restart
  }
  rhmc1[j].FA.solveopts = {
    prec = facg.prec,
    restart = facg.restart
  }
  rhmc1[j].MD.solveopts = {
    prec = mdcg.prec,
    restart = mdcg.restart,
    use_prev_soln = use_prev_soln
  }
  rhmc1[j].MD.ffprec = ffprec
end
myprint("rhmc1 = ", rhmc1, "\n")
copyto(rhmc, rhmc1)
myprint("rhmc = ", rhmc, "\n")

r.pbp = pbp
--myprint("runparams = ", r, "\n")

if inlat then
  acts:load(inlat)
else
  acts:unit()
end
local ps,pt = acts.fields.G:plaq()
printf("plaq ss: %g  st: %g  tot: %g\n", ps, pt, 0.5*(ps+pt))

acts:run(r)

if outlat then
  acts:save(outlat)
end
