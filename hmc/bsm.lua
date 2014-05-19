require 'common'
require 'run'

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
local grresid = grresid or faresid
local mdresid = mdresid or 1e-5
local restart = restart or 2000
local use_prev_soln = use_prev_soln or 0
local mixedRsq = mixedRsq or 0

local doGfix = doGfix or false
local doSpectrum = doSpectrum or false
local doS4obs = doS4obs or false

local rhmc = {}
local hmcmasses = hmcmasses or { mass, mass2 }
local seed = seed or os.time()

local inlat = inlat or nil
local outlat = outlat or nil

local ntraj = ntraj or 10
local tau = tau or 1
local ngpfs = ngpfs or 5
--lambdaG = lambdaG or 0.2
--lambdaF = lambdaF or 0.2
local nsteps = nsteps or 80
local nsteps2 = nsteps2 or 80

local ngsteps = ngsteps or ngpfs*nsteps
local nfsteps = nfsteps or { nsteps, nsteps2 }
nfsteps = repelem(nfsteps, nf/4)
local grcg = { prec=prec, resid=grresid, restart=restart }
local facg = { prec=prec, resid=faresid, restart=restart }
local mdcg = { prec=prec, resid=mdresid, restart=restart }
local ffprec = prec
--local gintalg = {type="leapfrog"}
--local gintalg = {type="omelyan"}
--local gintalg = {type="omelyan", lambda=0.2}
local gintalg = gintalg or {type="omelyan", lambda=lambdaG}
--local gintalg = {type="2MNV", lambda=lambdaG}
--local fintalg = {type="omelyan", lambda=0.2}
local fintalg = fintalg or {type="omelyan", lambda=lambdaF}
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
      printf("Setting pseudo: i=mass j=species -> %i\t%i\t%g\t%g\n",i,j,hmcmasses[i],hmcmasses[i+1])
      setpseudo(rhmc, hmcmasses[i], hmcmasses[i+1])
    else
      printf("Setting pseudo: i=mass j=species -> %i\t%i\t%g\n",i,j,hmcmasses[i])
      setpseudo(rhmc, hmcmasses[i])
    end
  end
end
local npseudo = #rhmc
myprint("Number of pseudo: ",npseudo,"\n")

local mdcgresid = {}
for i=1,#hmcmasses do
  mdcgresid[i] = (type(mdresid)=="table") and mdresid[i] or mdresid
end

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
  rhmc1[j].MD.resid = mdcgresid[1+math.floor(((j-1)*4+0.5)/nf)]
  rhmc1[j].GR[1].solveopts = {
    prec = grcg.prec,
    restart = grcg.restart,
    mixed_rsq = mixedRsq
  }
  rhmc1[j].FA.solveopts = {
    prec = facg.prec,
    restart = facg.restart,
    mixed_rsq = mixedRsq
  }
  rhmc1[j].MD.solveopts = {
    prec = mdcg.prec,
    restart = mdcg.restart,
    mixed_rsq = mixedRsq,
    use_prev_soln = use_prev_soln
  }
  rhmc1[j].MD.ffprec = ffprec
end
myprint("rhmc1 = ", rhmc1, "\n")
copyto(rhmc, rhmc1)
myprint("rhmc = ", rhmc, "\n")

if checkReverse then r.checkReverse = true end
r.pbp = pbp

function r.meas(a, r)
  local G = a.fields.G
  if doGfix then
    G.g:coulomb(3, 1e-6, 1500, 1.8);
    G.nupdate = G.nupdate + 1
  end

  if doSpectrum then
    -- Wall pions!
    for j=2,nt-1,8 do
      pions_wall = a.f:pions_wall(G, {j}, mass, 1e-6, {prec = prec, restart = 500})
      --pions_wall = a.f:pions_wall(G,{2,10,18,26},mass,1e-6, {prec = prec, restart = 5000})

      printf("Source %i\n", j);
      printf("Local Pions\n");
      printf("\tPion5\t\tPion5_4\n");
      for i = 1,#(pions_wall.pion5) do
	printf("%i\t%.6e\t%.6e\n", i-1, pions_wall.pion5[i], pions_wall.pion5_gamma4[i]);
      end

      printf("Local Baryon\n")
      printf("\tNucleon\n")
      for i = 1,#(pions_wall.pion5) do
	printf("%i\t%.6e\n", i-1, pions_wall.nucleon[i].r)
	--printf("%i\t%s\n", i-1, tostring(pions_wall.nucleon[i]))
      end

      printf("\nNonlocal Pions + Check\n")
      printf("\tPion5\tPion5_4\tPion_i5\tPion_ij\n");
      for i = 1,#(pions_wall.pion5) do
	printf("%i\t%.6e\t%.6e\t%.6e\t%.6e\n", i-1, pions_wall.pion5_ck[i],
	       pions_wall.pion5_gamma4_ck[i], pions_wall.pion_i5[i], pions_wall.pion_ij[i]);
      end

      printf("Nonlocal Baryons\n");
      printf("\tNucleon\tDelta\n");
      for i = 1,#(pions_wall.pion5) do
	printf("%i\t%.6e\t%.6e\n", i-1, pions_wall.nucleon_ck[i].r, pions_wall.delta[i].r);
      end
    end
  end

  if doS4obs then
    local masses = {mass, 2*mass};
    local s4_check = a.f:s4_broken_observe(G, masses, 1e-6, {prec = prec, restart = 5000}, 1)

    printf("MEASplaq_ss %.12e\n", s4_check.s4_g_plaq);
    printf("MEASplaq_t even %.12e odd %.12e\n", s4_check.s4_g_even[1], s4_check.s4_g_odd[1]);
    printf("MEASplaq_x even %.12e odd %.12e\n", s4_check.s4_g_even[2], s4_check.s4_g_odd[2]);
    printf("MEASplaq_y even %.12e odd %.12e\n", s4_check.s4_g_even[3], s4_check.s4_g_odd[3]);
    printf("MEASplaq_z even %.12e odd %.12e\n", s4_check.s4_g_even[4], s4_check.s4_g_odd[4]);
    printf("MEASplaq_a even %.12e odd %.12e\n", s4_check.s4_g_even[5], s4_check.s4_g_odd[5]);

    for i,v in ipairs(masses) do
      printf("MEASpbp_all m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][1].r, s4_check.s4_f_even[i][1].i, s4_check.s4_f_odd[i][1].r, s4_check.s4_f_odd[i][1].i);
      printf("MEASpbp_1 m %.4e even %.12e %.12e odd %.12e %.12e\n",  v,s4_check.s4_f_even[i][2].r, s4_check.s4_f_even[i][2].i, s4_check.s4_f_odd[i][2].r, s4_check.s4_f_odd[i][2].i);
      printf("MEASpbp_t m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][3].r, s4_check.s4_f_even[i][3].i, s4_check.s4_f_odd[i][3].r, s4_check.s4_f_odd[i][3].i);
      printf("MEASpbp_x m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][4].r, s4_check.s4_f_even[i][4].i, s4_check.s4_f_odd[i][4].r, s4_check.s4_f_odd[i][4].i);
      printf("MEASpbp_y m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][5].r, s4_check.s4_f_even[i][5].i, s4_check.s4_f_odd[i][5].r, s4_check.s4_f_odd[i][5].i);
      printf("MEASpbp_z m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][6].r, s4_check.s4_f_even[i][6].i, s4_check.s4_f_odd[i][6].r, s4_check.s4_f_odd[i][6].i);
    end
  end
end

--myprint("runparams = ", r, "\n")

local traj = traj or 0

if inlat then
  acts:load(inlat)
else
  if latpat then
    if traj > 0 then
      inlat = string.format(latpat, traj)
      acts:load(inlat)
    else
      acts:unit()
    end
  else
    acts:unit()
  end
end
if qopverb then qopqdp.verbosity(qopverb) end

local ps,pt = acts.fields.G:plaq()
printf("plaq ss: %g  st: %g  tot: %g\n", ps, pt, 0.5*(ps+pt))

local nlats = nlats or 1

--[[
L = qopqdp.defaultLattice()
rdr = L:reader("test.rs")
rs = L:getRstate()
rs:read(rdr)
rdr:close()
--]]

for nl=1,nlats do
  if warmup and traj<warmup then
    r.md = true
  else
    r.md = false
  end

  r.ntraj = ntraj
  acts:run(r)
  traj = traj + ntraj
  if latpat then
    outlat = string.format(latpat, traj)
  end
  if outlat then
    acts:save(outlat)
  end
end

--[[
L = qopqdp.defaultLattice()
rs = L:getRstate()
w = L:writer("test.rs","metadata")
rs:write(w,"metadata")
w:close()
--]]
