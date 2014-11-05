package.path = (userpath or "") .. "/projectnb/qcd/Staggered8f/qhmc/?.lua;/projectnb/qcd/Staggered8f/qhmc/hmc/?.lua;" .. arg[0]:gsub("[^/]*.lua","?.lua;") .. package.path
require 'common'
require 'gaugeact'
require 'asqtadact'
require 'hmc'

--trace(doTrace)

-- General lattice parameters.
local nx = 16
local nt = 32
local beta = 5.0
local beta_a = -0.25 -- for adjoint gauge action
local u0 = 1.0
local nf = 8
local mass = 0.05
local prec = 1 -- Hm? Oh, single precision.

-- So I feel like these aren't needed anymore.
local ntraj = 10
local tau = 1
local nsteps = 120
seed = seed or os.time()
local first = -1
local ffreq = 3

-- I added this!
local inlat = "cfg/hyp1632f8b50m05a_coulomb.1134.lime"

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
p.gaugeact = {type="plaquette_adjoint", adjFac=beta_a}
--p.gaugeact = {type="plaquette", u0=u0}
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
--setpseudo(rhmc, mass2, mass)

local smear = {}

function projRat(s)
  local o = s.coeffs.one_link
  s.coeffs.one_link = nil
  return { type="projectRat", rho=1/o, smear=s }
  --return { type="projectRat", rho=4, smear=s }
end

-- I'm just going to trust these.
scoeffs = asqtad_coeffs(u0)
scoeffs.three_staple = -scoeffs.three_staple
scoeffs.seven_staple = -scoeffs.seven_staple
scoeffs.one_link = scoeffs.one_link + 3*scoeffs.naik
scoeffs.naik = nil
scoeffs.one_link = scoeffs.one_link + 6*scoeffs.lepage
scoeffs.lepage = nil


smear[#smear+1] = { type="hyp", alpha={0.4,0.5,0.5} }
myprint("smear = ", smear, "\n")

coeffs = { one_link=1 } -- presumably no smearing.

local act = {}
act.g = gaugeact(p)
act.f = asqtadact(act.g, {smear=smear,coeffs=coeffs,rhmc=rhmc})

G = act.g:gaugeNew(act)
last = ntraj * tau
G:load(inlat)

--G:save("better_lat.lime");


do
  local devavg,devmax = G.g:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  G.g:makeSU()
  devavg,devmax = G.g:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end


function measure(G)
  local t0 = clock()
  local ps,pt = G:plaq()
  printf("plaq ss: %-8g  st: %-8g  tot: %g\n", ps, pt, 0.5*(ps+pt))

  local nd = #G.a.latsize
  for i=1,nd do
    for j=i+1,nd do
      --local lp = G.g:loop({i,j,-i,-j}) -- changed convention
      local lp = G.g:loop({-i,-j,i,j})
      printf("plaq%i%i:  %g\t%g\n", i, j, lp.r, lp.i)
    end
  end

  --local plpath = rep(-nd, G.a.latsize[nd]) -- changed convention
  local plpath = rep(nd, G.a.latsize[nd])
  plp = G.g:loop(plpath)
  printf("ploop:  %g\t%g\n", plp.r, plp.i)

  t0 = clock() - t0
  printf("meas time: %g\n", t0)

  io.stdout:flush()
end



totaltime = clock()
measure(G)

-- Let's try doing something fermion related.
-- So... action, gauge field, mass, residual, {single/double precision, max iters} ?
--local pbp1, pbp2 = act.f:pbp(G, mass, 1e-4, {prec = prec, restart = 1000})

--printf("PBP results: %g %g\n", pbp1, pbp2);

-- Okay, cool! So we have access to a :dot(other), similar to :norm2().

-- Dot versus norm2.
--[[squark_field = act.f.h:quark();
squark_field:zero();
squark_field:point({0,0,0,0}, 1, 0, 1); -- Set to i.
local from_norm2 = squark_field:norm2();
local from_dot = squark_field:Re_dot(squark_field);
printf("From norm: %g, From dot: %g\n", from_norm2, from_dot);
io.stdout:flush();]]


--local t0 = clock();

--G.g:gaugefix(3, 1e-4, 1000, 1, 1.2);

--t0 = clock() - t0
--printf("Coulgauge meas time: %g\n", t0)

io.stdout:flush()

--measure(G);

-- Pions!
--local pions = act.f:pions(G,{2,10,18,26},mass,1e-7, {prec = prec, restart = 5000})
--myprint("Pion5", pions.pion5, "\n");
--myprint("Rhoi", pions.rhoi, "\n");
--myprint("PionS", pions.pionS, "\n");

-- S4 broken phase check.

-- S4 broken phase check.

-- Arguments:
-- 1. gauge action, NOT FIELD. (should be called G in your code, too)
-- 2. The mass you want to use to measure pbp. Does not support more than one mass at once, so it must be a single value, i.e., 0.005.
-- 3. Stopping condition of CG.
-- 4. Other parameters of CG. prec = 1 for single precision, 2 for double, and restart is the restart value.
-- 5. The number of stoichastic pbp measurements to perform. I just did 1, but its your call. It automatically averages over all of them.
-- Note: the gauge observables use the nHYP smeared gauge links, not the bare ones. If they don't agree with what other plaquette measurements, that's why!
--[[local masses = {mass, 2*mass};
local s4_check = act.f:s4_broken_observe(G, masses, 1e-6, {prec = prec, restart = 5000}, 1)

printf("MEASplaq_ss %.12e\n", s4_check.s4_g_plaq);
printf("MEASplaq_t e %.12e o %.12e\n", s4_check.s4_g_even[1], s4_check.s4_g_odd[1]);
printf("MEASplaq_x e %.12e o %.12e\n", s4_check.s4_g_even[2], s4_check.s4_g_odd[2]);
printf("MEASplaq_y e %.12e o %.12e\n", s4_check.s4_g_even[3], s4_check.s4_g_odd[3]);
printf("MEASplaq_z e %.12e o %.12e\n", s4_check.s4_g_even[4], s4_check.s4_g_odd[4]);
printf("MEASplaq_a e %.12e o %.12e\n", s4_check.s4_g_even[5], s4_check.s4_g_odd[5]);

for i,v in ipairs(masses) do
	printf("MEASpbp_all m %.4e e %.12e %.12e o %.12e %.12e\n", v, s4_check.s4_f_even[i][1].r, s4_check.s4_f_even[i][1].i, s4_check.s4_f_odd[i][1].r, s4_check.s4_f_odd[i][1].i);
	printf("MEASpbp_1 m %.4e e %.12e %.12e o %.12e %.12e\n",  v,s4_check.s4_f_even[i][2].r, s4_check.s4_f_even[i][2].i, s4_check.s4_f_odd[i][2].r, s4_check.s4_f_odd[i][2].i);
	printf("MEASpbp_t m %.4e e %.12e %.12e o %.12e %.12e\n", v, s4_check.s4_f_even[i][3].r, s4_check.s4_f_even[i][3].i, s4_check.s4_f_odd[i][3].r, s4_check.s4_f_odd[i][3].i);
	printf("MEASpbp_x m %.4e e %.12e %.12e o %.12e %.12e\n", v, s4_check.s4_f_even[i][4].r, s4_check.s4_f_even[i][4].i, s4_check.s4_f_odd[i][4].r, s4_check.s4_f_odd[i][4].i);
	printf("MEASpbp_y m %.4e e %.12e %.12e o %.12e %.12e\n", v, s4_check.s4_f_even[i][5].r, s4_check.s4_f_even[i][5].i, s4_check.s4_f_odd[i][5].r, s4_check.s4_f_odd[i][5].i);
	printf("MEASpbp_z m %.4e e %.12e %.12e o %.12e %.12e\n", v, s4_check.s4_f_even[i][6].r, s4_check.s4_f_even[i][6].i, s4_check.s4_f_odd[i][6].r, s4_check.s4_f_odd[i][6].i);
end
]]


-- Wall pions!
--for j=2,26,8 do
	pions_wall = act.f:pions_wall(G,{2,10,18,26},mass,1e-6, {prec = prec, restart = 5000})
	--pions_wall = act.f:pions_wall(G,{2,10,18,26},mass,1e-6, {prec = prec, restart = 5000})

	--printf("Source %i\n", j);
-- This is all to reproduce what MILC does.
	printf("STARTPROP\n");
	printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
	printf("SOURCE: EVEN_WALL\n");
	printf("SINKS: PION_PS PION_SC\n");
	--printf("Local Pions\n");
	--printf("\tPion5\t\tPion5_4\n");

	for i = 1,#(pions_wall.pion5) do
	   printf("%i %.6e %f %.6e %f\n", i-1, pions_wall.pion5[i], 0.0, pions_wall.pion5_gamma4[i], 0.0);
	end

	printf("ENDPROP\n");
	
	printf("STARTPROP\n");
        printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
        printf("SOURCE: EVEN_WALL\n");
        printf("SINKS: NUCLEON\n");


	--printf("Local Baryon\n");
	--printf("\tNucleon\n");

	for i = 1,#(pions_wall.pion5) do
	   printf("%i %.6e %f\n", i-1, pions_wall.nucleon[i].r, 0.0);
	end

	printf("ENDPROP\n");
        
        printf("STARTPROP\n");
        printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
        printf("SOURCE: EVENANDODD_WALL\n");
        printf("SINKS: PION_PS PION_SC PION_i5 PION_ij\n");


	--printf("\nNonlocal Pions + Check\n")
	--printf("\tPion5\tPion5_4\tPion_i5\tPion_ij\n");

	for i = 1,#(pions_wall.pion5) do
	   printf("%i %.6e %f %.6e %f %.6e %f %.6e %f\n", i-1, pions_wall.pion5_ck[i], 0.0, pions_wall.pion5_gamma4_ck[i], 0.0, pions_wall.pion_i5[i], 0.0, pions_wall.pion_ij[i], 0.0);
	end
	
		printf("STARTPROP\n");
        printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
        printf("SOURCE: EVENANDODD_WALL\n");
        printf("SINKS: RHO_0 RHO_is RHO_ij RHO_i5\n");

	for i = 1,#(pions_wall.pion5) do
	   printf("%i %.6e %f %.6e %f %.6e %f %.6e %f\n", i-1, pions_wall.rho_0[i], 0.0, pions_wall.rho_is[i], 0.0, pions_wall.rho_ij[i], 0.0, pions_wall.rho_i5[i], 0.0);
	end

	printf("ENDPROP\n");

        printf("STARTPROP\n");
        printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
        printf("SOURCE: EVENANDODD_WALL\n");
        printf("SINKS: NUCLEON DELTA\n");

	
	--printf("Nonlocal Baryons\n");
	--printf("\tNucleon\tDelta\n");

	for i = 1,#(pions_wall.pion5) do
	   printf("%i %.6e %f %.6e %f\n", i-1, pions_wall.nucleon_ck[i].r, 0.0, pions_wall.delta[i].r, 0.0);
	end

	printf("ENDPROP\n");

--end

--myprint("Pion5", pions_wall.pion5, "\n");
--myprint("Pion5_Gamma4", pions_wall.pion5_gamma4, "\n");
--myprint("Pion5_i5", pions_wall.pion_i5, "\n");

--[[
squark_field = act.f.h:quark();
squark_field:zero();
squark_field:point({0,0,0,0}, 1, 0.5, -0.3);
squark_field:point({1,0,0,0},1,0.5,-0.3);
squark_field:printpoint({0,0,0,0});
squark_field:printpoint({1,0,0,0});
printf("Rephase.\n");
io.stdout:flush();
squark_field:rephase(1, {1,0,0,0});
squark_field:printpoint({0,0,0,0});
squark_field:printpoint({1,0,0,0}); ]]

io.stdout:flush()


-- Measure CG info. 

--[[
do

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

	act.f:clearStats()

	io.stdout:flush();

end

totaltime = clock() - totaltime
printf("total time: %g seconds\n", totaltime)]]

