require 'common'
require 'run'

function milcin(fn)
  local t = readfile(fn, {["warms"]={l=1,t="runs"}})
  return t
end

function milcrat(fn)
  local t = readfile(fn, {["naik_term_epsilon"]={l=1,t="pseudo"}})
  return t
end

local infile = infile or arg[1]
if not infile then
  printf("no input file: %s %s\n", arg[-1], arg[0])
  exit(1)
end
printf("reading input file %s\n", infile)
local mi = milcin(infile)
--myprint(mi)

local mr = nil
if mi["load_rhmc_params"] then
  local ratfile = mi["load_rhmc_params"][1]
  mr = milcrat(ratfile)
  --myprint(mr)
end

ngpfs = ngpfs or 3
lambdaG = lambdaG or 0.2
lambdaF = lambdaF or 0.2

local p = {}
p.latsize = { mi.nx[1], mi.ny[1], mi.nz[1], mi.nt[1] }
p.seed = mi.iseed[1]
p.beta = mi.beta[1]
p.nf = sum(mi.dyn_flavors)
p.u0 = mi.u0[1]
--p.gaugeact = {type="symanzik_1loop", u0=p.u0, nf=p.nf}
p.gaugeact = {type="symanzik_1loop_hisq", u0=p.u0, nf=p.nf}
p.npseudo = #mr.pseudo
local rhmc = {}
--p.fermact = {type="asqtad", u0=p.u0, rhmc=rhmc}
p.fermact = {type="hisq", rhmc=rhmc}
for i,v in ipairs(mr.pseudo) do
  v.pole_GR[1] = nil
  v.pole_FA[1] = nil
  v.pole_MD[1] = nil
  rhmc[i] = {}
  rhmc[i].GR = {}
  rhmc[i].FA = {}
  rhmc[i].GR[1] = paste(v.res_GR, v.pole_GR)
  rhmc[i].FA[1] = paste(v.res_FA, v.pole_FA)
  rhmc[i].MD = paste(v.res_MD, v.pole_MD)
end
local rhmc0 = copy(rhmc)

local acts = setupacts(p)

myprint("rhmc0 = ", rhmc0, "\n")

for i,v in ipairs(mi.runs) do
  local r = {}
  r.ntraj = v.trajecs[1]
  local nsteps = v.steps_per_trajectory[1]
  r.tau = nsteps*v.microcanonical_time_step[1]
  local intalgG = {type="omelyan", lambda=lambdaG}
  local intalgF = {type="omelyan", lambda=lambdaF}
  local fp = {}
  r.forceparams = fp
  fp[1] = {}
  fp[1][1] = {nsteps=ngpfs*nsteps, intalg=intalgG}
  local rhmc1 = {}
  for j=1,p.npseudo do
    fp[1][j+1] = {nsteps=nsteps, intalg=intalgF}
    rhmc1[j] = {GR={},FA={},MD={}}
    rhmc1[j].GR[1] = {}
    rhmc1[j].GR[1].resid = tonumber(v.cgresid_md_fa_gr[3*(j-1)+3])
    rhmc1[j].FA.resid = tonumber(v.cgresid_md_fa_gr[3*(j-1)+2])
    rhmc1[j].MD.resid = tonumber(v.cgresid_md_fa_gr[3*(j-1)+1])
    rhmc1[j].GR[1].solveopts = {
      prec = v.cgprec_md_fa_gr[3*(j-1)+3],
      restart = v.max_multicg_md_fa_gr[3*(j-1)+3]
    }
    rhmc1[j].FA.solveopts = {
      prec = v.cgprec_md_fa_gr[3*(j-1)+2],
      restart = v.max_multicg_md_fa_gr[3*(j-1)+2]
    }
    rhmc1[j].MD.solveopts = {
      prec = v.cgprec_md_fa_gr[3*(j-1)+1],
      restart = v.max_multicg_md_fa_gr[3*(j-1)+1],
      use_prev_soln = use_prev_soln
    }
    rhmc1[j].MD.ffprec = v.prec_ff[1]
  end
  myprint("rhmc1 = ", rhmc1, "\n")
  copyto(rhmc, rhmc1)

  local pbpreps = (v.npbp_reps or {1})[1]
  local pbpmasses = v.mass or mi.dyn_mass
  local pbpresid = tonumber(v.error_for_propagator) or { 1e-6 }
  if #pbpresid ~= #pbpmasses then
    pbpresid = rep(pbpresid[1], #pbpmasses)
  end
  local cgmax = v.max_cg_prop[1]
  local cgrest = v.max_cg_prop_restarts[1]
  local pbpopts = { restart=cgmax, max_restarts=cgrest, max_iter=cgmax*cgrest }
  r.pbp = {}
  for i=1,#pbpmasses do
    local t = {}
    t.reps = pbpreps
    t.mass = pbpmasses[i]
    t.resid = pbpresid[i]
    t.opts = pbpopts
    r.pbp[i] = t;
  end
  if checkReverse then r.checkReverse = true end
  myprint("runparams = ", r, "\n")

  if(v.fresh) then
    acts:unit()
  end
  for k,f in pairs(v) do
    if(string.match(k,"^reload")) then
      acts:load(f[1])
      local ps,pt = acts.fields.G:plaq()
      printf("plaq ss: %g  st: %g  tot: %g\n", ps, pt, 0.5*(ps+pt))
    end
  end
  acts:run(r)
  for k,f in pairs(v) do
    if(string.match(k,"^save")) then
      acts:save(f[1])
    end
  end
end
