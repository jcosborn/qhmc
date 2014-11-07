-- ESW 09-24-2014
-- This standalone file gauge fixes a configuration to coulomb gauge, 
-- then measure various spectral quantities.

require 'common'
require 'smear'
require 'run'

-- Parameters we need to set.
--[[
local nx = 24 -- Lattice spatial dimension
local nt = 48 -- Lattice temporal dimension
local mass = 0.005 -- Light quark mass.
local inlat = "configurations/f4plus8l24t48b40m005m080_:CONFIGNUM:_scidac" -- Input configuration
local gfix_prec = 1e-7 -- Max residual of gauge fixing.
local gfix_max = 4000 -- Maximum number of gauge fixing steps.
local cg_prec = 1e-6 -- Max residual of CG.
local cg_max = 4000 -- Maximum number of CG steps.
local src_start = 2 -- Where to put the first source.
local src_num = math.floor(nt/8) -- How many sources to place. (From this, a uniform spacing is derived)
local gfix_or = 1.75 -- The overrelaxation parameter for gauge fixing.
--]]

local nx = nx or 4
local nt = nt or 8
local mass = mass or 0.05
--local inlat = "hyp1632f8b50m05a_coulomb.1002.lime"
local gfix_prec = gfix_prec or 1e-6
local gfix_max = gfix_max or 200
local gfix_or = gfix_or or 1.75
local cg_prec = cg_prec or 1e-6
local cg_max = cg_max or 2000
local src_start = src_start or 2
local src_num = src_num or math.floor(nt/8)
local prec = prec or 1 -- Use a multi-precision inverter. Default.

-- Start preparing to load the gauge field, spit out basic info.
local latsize = { nx, nx, nx, nt }
local vol = 1
local spatvol = nx*nx*nx;
local seed = seed or os.time()
printf("latsize =")
for k,v in ipairs(latsize) do vol=vol*v; printf(" %i",v); end
printf("\nvolume = %i\n", vol)
printf("mass = %g\n", mass)
printf("seed = %i\n", seed)
printf("gfix_prec = %g\n", gfix_prec)
printf("cg_prec = %g\n", cg_prec)
printf("src_num = %i\n", src_num)

-- Set up qopqdp.
qopqdp.lattice(latsize);
qopqdp.profile(profile or 0);
qopqdp.verbosity(0);
qopqdp.seed(seed);

-- Start a timer.
totaltime = qopqdp.dtime()

-- Load the gauge field.
g = qopqdp.gauge();
if inlat then g:load(inlat)
else
  g:random()
end

-- Reunitarize, just to be safe.
do
  local devavg,devmax = g:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  g:makeSU()
  devavg,devmax = g:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end

-- Print some basic information about the configuration.
function getplaq(g)
  local ps,pt = g:action{plaq=1}
  local lat = qopqdp.lattice()
  local nd,vol = #lat,1
  for i=1,nd do vol=vol*lat[i] end
  local s = 0.25*nd*(nd-1)*vol
  printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", ps/s, pt/s, 0.5*(ps+pt)/s)
end
getplaq(g);


-- First, we need to gauge fix!
local t0 = qopqdp.dtime();

-- coulomb(j_decay, error, max iterations, overrelaxation param)
-- note that here 0->x, ..., 3->t
g:coulomb(3, gfix_prec, gfix_max, gfix_or);

t0 = qopqdp.dtime() - t0
printf("Coulgauge meas time: %g\n", t0)

-- Prepare a smearing setting. This has been verified to be consistent
-- with MILC. Then smear the gauge field!
local smear = {}
smear[#smear+1] = { type="hyp", alpha={0.4,0.5,0.5} }
myprint("smear = ", smear, "\n")

printf("Start smearing.\n");

-- We need to do this because 'smearGauge' expects an action object.
local sg = smearGauge({g = g}, smear);
printf("Smearing done.\n");

-- Set the ASQTAD coefficients. This corresponds to just simple
-- Staggered fermions (supplemented with nHYP smearing)
coeffs = { one_link=1 }

-- Create an asqtad object, set coefficients, etc.
w = qopqdp.asqtad();
w:coeffs(coeffs);
w:set(sg, prec); 

-- By the way...
local Nc = qopqdp.defaultNc();

-- Now we can prepare to measure things! Measurements are
-- set up to output in a form similar to MILC's measurements.

-- Set up a function which does an inversion and prints
-- out timing information. Based on actmt.set in asqtadact.lua
function solve_printinfo(w, dest, src, m, res, sub, opts)
	local t0 = qopqdp.dtime();
	w:solve({dest}, src, {m}, res, sub, opts);
	local cgtime = qopqdp.dtime() - t0;
	local flops = w:flops();
	local its = w:its();
	
	printf("inversion its: %g  secs: %g  Mflops: %g\n", its, cgtime, flops*1e-6/cgtime);
end

-- Set up a function which adds a correlator to an existing
-- correlator array, properly shifting, normalizing, and 
-- compensating for mesonic/baryonic source.
-- corr: What we're adding a new piece of data to.
-- new_data: the new correlator we're adding in.
-- t_shift: the 't' value of the source.
-- norm: the normalization to multiply new_data by before adding it.
function add_correlator_meson(corr, new_data, t_shift, norm)

	for j=1,#new_data do
		-- Compensate for shifted wall source.
		j_real = (j+t_shift-1)%(#new_data)+1 
		corr[j] = new_data[j_real]*norm + (corr[j] or 0)
	end

end

-- Same deal, but for baryons we have to be careful about
-- going around the T end of the lattice.

function add_correlator_baryon(corr, new_data, t_shift, norm)

	for j=1,#new_data do
		-- Compensate for shifted wall source.
		j_real = (j+t_shift-1)%(#new_data)+1 
		if (((math.floor((j+t_shift-1)/(#new_data))-math.floor((t_shift-1)/(#new_data)))%2) == 0) then -- count the number of times we wrap.
			corr[j] = -new_data[j_real]*norm + (corr[j] or 0)
		else
			corr[j] = new_data[j_real]*norm + (corr[j] or 0)
		end
	end

end

-- Set up a function to build staggered phases.
-- Use 0 for no phase, 1 for phase.
function make_phase_term(xsign,ysign,zsign,tsign)
  return xsign + 2*ysign + 4*zsign + 8*tsign
end

-- First, set up where we put sources. Evenly space src_num sources
-- starting at t=src_start.
local time_sources = {};
for i=1,src_num do
	time_sources[i] = (src_start+(i-1)*math.floor(nt/src_num))%nt;
end

-- Next, measure the pion using a random wall source (equiv to
-- a point source in the infinite stochastic source limit). 
-- This is equivalent to the asqtadact function, fpi_pion, that I wrote.
-- Well... at least in function.

-- We want 3 random sources.
num_rand = 3;

-- We also need two quark objects---one to put the source in,
-- the other the sink.
src = w:quark();
src:zero();
dest = w:quark();
dest:zero();

-- A place for correlator output to temporarily go.

t = {}; 

pion_fpi_ptp = {};  -- Where point to point fpi goes.
pion_fpi_ptw = {};  -- Where point to wall fpi goes.

do
	-- Loop over all time sources and all stochastic sources.
	for srcnum=1,#time_sources do
		for randnum=1,num_rand do
			printf("Start random wall source %i at t=%i.\n", randnum, time_sources[srcnum]);
			io.stdout:flush();
			
			-- Create a random wall source (gaussian)
			src:wall_gaussian(time_sources[srcnum]);
			printf("gauss_src norm2: %g\n", src:norm2());
			
			-- And invert! This is a special function which also prints
			-- out CG info.
			solve_printinfo(w, dest, src, mass, cg_prec, "all", {prec = prec, restart = cg_max});
			
			-- Contract, point source. 
			t = dest:norm2("timeslices");
			add_correlator_meson(pion_fpi_ptp, t, time_sources[srcnum], 1.0/(#time_sources*num_rand*spatvol*spatvol));
			
			-- Contract, wall source.
			t = dest:norm2_wallsink("timeslices");
			add_correlator_meson(pion_fpi_ptw, t, time_sources[srcnum], 1.0/(#time_sources*num_rand*spatvol*spatvol));
		end
	end
end	

-- Next, prepare wall source solves.
-- We need a lot more quarks. Let's create these here!

-- Reuse some old things and allocate new ones.
even_src = src; even_src:zero();
odd_src = dest; even_src:zero();
odd_soln = w:quark(); odd_soln:zero();
o_gupta = w:quark(); o_gupta:zero();
Do_gupta = w:quark(); Do_gupta:zero();
Dq_gupta = w:quark(); Dq_gupta:zero();
temp1 = w:quark(); temp1:zero();
temp2 = w:quark(); temp2:zero();
-- We need to save even results for baryon contractions.
even_soln, q_gupta = {}, {};
Dq_gupta_all = {{}, {}, {}};
for i=1,3 do
	even_soln[i] = w:quark(); even_soln[i]:zero();
	q_gupta[i] = w:quark(); q_gupta[i]:zero();
	for j=1,3 do
		Dq_gupta_all[i][j] = w:quark(); Dq_gupta_all[i][j]:zero();
	end
end

-- Since we have guage fixed, we don't need to use the
-- gauge field in the parallel transporter.
local unitg = qopqdp.gauge();
unitg:unit();

-- Stick the results somewhere!
p5,p5_g4,pion_ps_ck, pion_4_ck, pion_i5,pion_ij = {},{},{},{},{},{}
rho_0, rho_is, rho_ij, rho_i5 = {}, {}, {}, {}
nucleon, nucleon_ck, delta = {}, {}, {} 
		
-- Let's go go go!
do
	-- Loop over all time sources.
	for srcnum=1,#time_sources do
		printf("Start wall source %i at t=%i.\n", srcnum, time_sources[srcnum]);
		
		-- Loop over all colors. This is for meson measurements,
		-- and preparing for nucleon measurements.
		
		for i=1,Nc do
			printf("Start color %i.\n", i);
			io.stdout:flush();
			
			-- Prepare the even source. Even wall, Norm matches milc.
			even_src:zero();
			even_src:wall(time_sources[srcnum], 0, i, -0.125); 
			printf("even_src norm2: %g\n", even_src:norm2());
			
			-- Prepare the odd source. Odd wall, Norm matches milc.
			odd_src:zero();
			odd_src:wall(time_sources[srcnum], 1, i, -0.125);
			printf("odd_src norm2: %g\n", odd_src:norm2());
			
			-- Invert on the even source. Note that we save each
			-- even solution for computing the nucleon.
			solve_printinfo(w, even_soln[i], even_src, mass, cg_prec, "all", {prec = prec, restart = cg_max});
			
			-- Invert on the odd source. 
			solve_printinfo(w, odd_soln, odd_src, mass, cg_prec, "all", {prec = prec, restart = cg_max});
			
			-- First off, let's follow MILC in computing the local
			-- pions. These are the Goldstone pion, and the extra
			-- gamma_4 pion. These are only printed for a
			-- consistency check.
			
			do -- Pion 0_A^(-+)
				t = even_soln[i]:norm2("timeslices");
				add_correlator_meson(pion_ps_ck, t, time_sources[srcnum], 1.0/(#time_sources));
			end

			do -- Pion 0_A^(-+) w/ extra gamma_4
				temp1:set(even_soln[i]); -- Make a copy. This is for rephasing.
				temp1:rephase(make_phase_term(1,1,1,0), {0,0,0,time_sources[srcnum]}) -- (-1)^(x+y+z)
				t = even_soln[i]:Re_dot(temp1, "timeslices");
				add_correlator_meson(pion_4_ck, t, time_sources[srcnum], 1.0/(#time_sources));
			end
			
			-- In the language of Gupta et al, the even_soln is "o+q".
			-- The odd_soln is "q-o". We can thus reconstruct Gupta's
			-- o, q.
			q_gupta[i]:zero();
			o_gupta:zero();
			q_gupta[i]:combine({even_soln[i], odd_soln},{1.0,1.0});
			o_gupta:combine({even_soln[i], odd_soln},{1.0,-1.0});
			
			-- We next prepare symmetric shifted values in the z direction.
			-- We then zero out odd z coordinates. These are in the wrong
			-- place with respect to the 2^4 hypercube, and will contribute
			-- to a wrong answer.
			
			-- First, symm shift the odd solution in the z-direction.
			-- Then remove odd z values.
			Do_gupta:symshift(o_gupta, unitg, 3);
			temp1:set(Do_gupta);
			temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,time_sources[srcnum]})
			temp2:set(Do_gupta);
			Do_gupta:combine({temp1, temp2}, {0.5,0.5});
			
			-- Do the same to the even solution.
			Dq_gupta:symshift(q_gupta[i], unitg, 3);
			temp1:set(Dq_gupta);
			temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,time_sources[srcnum]})
			temp2:set(Dq_gupta);
			Dq_gupta:combine({temp1, temp2}, {0.5,0.5});
			
			-- Using these states, we can compute the one-link separated 
			-- states denoted pion_i5 and pion_ij.

			-- First, pion_i5. This requires no phasing, just q Dq - o Do.

			do 
				-- q Dq
				t = q_gupta[i]:Re_dot(Dq_gupta, "timeslices");
				add_correlator_meson(pion_i5, t, time_sources[srcnum], 1.0/(#time_sources));
				
				-- minus o Do

				t = o_gupta:Re_dot(Do_gupta, "timeslices");
				add_correlator_meson(pion_i5, t, time_sources[srcnum], -1.0/(#time_sources));
			end

			-- Next, pion_ij. This requires rephasing by a factor of (-1)^(x+y+z),
			-- then we take o Dq - q Do.
			
			Dq_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]});
			Do_gupta:rephase(make_phase_term(1,1,1,0), {0,0,0,time_sources[srcnum]});
			do 
				-- q Dq
				t = o_gupta:Re_dot(Dq_gupta, "timeslices");
				add_correlator_meson(pion_ij, t, time_sources[srcnum], 1.0/(#time_sources));
				
				-- minus o Do

				t = q_gupta[i]:Re_dot(Do_gupta, "timeslices");
				add_correlator_meson(pion_ij, t, time_sources[srcnum], -1.0/(#time_sources));
			end


			-- Better values for the pion and the g4_pion. Taking
			-- combinations of q, o reduces errors.

			-- First, the Goldstone boson pion. This is q q + o o.
			do
				-- q q
				t = q_gupta[i]:Re_dot(q_gupta[i], "timeslices");
				add_correlator_meson(p5, t, time_sources[srcnum], 1.0/(#time_sources));
				
				-- plus o o

				t = o_gupta:Re_dot(o_gupta, "timeslices");
				add_correlator_meson(p5, t, time_sources[srcnum], 1.0/(#time_sources));

			end

			-- Next, the gamma4 pion. This is q o with a (-1)^(x+y+z) phase factor.
			temp1:set(o_gupta);
			temp1:rephase(make_phase_term(1,1,1,0), {0,0,0,time_sources[srcnum]})
			do
				t = q_gupta[i]:Re_dot(temp1, "timeslices");
				add_correlator_meson(p5_g4, t, time_sources[srcnum], 1.0/(#time_sources));
			end

			-- There was nothing special about shifting in the z direction.
			-- We can also shift in the x, y direction. This gets us
			-- some rho mesons!
			
			-- First, x shifted. This lets us grab:
			-- 1. rho_0: gamma_1 gamma_4 x taste_4
			-- 2. rho_is: gamma_1 x 1
			
			-- Odd soln.
			Do_gupta:symshift(o_gupta, unitg, 1);
			temp1:set(Do_gupta);
			temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,time_sources[srcnum]})
			temp2:set(Do_gupta);
			Do_gupta:combine({temp1, temp2}, {0.5,0.5});
			
			-- Do the same to the even solution.
			Dq_gupta:symshift(q_gupta[i], unitg, 1);
			temp1:set(Dq_gupta);
			temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,time_sources[srcnum]})
			temp2:set(Dq_gupta);
			Dq_gupta:combine({temp1, temp2}, {0.5,0.5});
			
			-- The combinations look familiar!
			-- First, rho_0. This requires no phasing, just q Dq - o Do.

			do 
				-- q Dq
				t = q_gupta[i]:Re_dot(Dq_gupta, "timeslices");
				add_correlator_meson(rho_0, t, time_sources[srcnum], 1.0/(#time_sources));
				
				-- minus o Do

				t = o_gupta:Re_dot(Do_gupta, "timeslices");
				add_correlator_meson(rho_0, t, time_sources[srcnum], -1.0/(#time_sources));
			end

			-- Next, rho_is. This requires rephasing by a factor of (-1)^(x+y+z),
			-- then we take o Dq - q Do.
			
			Dq_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]});
			Do_gupta:rephase(make_phase_term(1,1,1,0), {0,0,0,time_sources[srcnum]});
			do 
				-- q Dq
				t = o_gupta:Re_dot(Dq_gupta, "timeslices");
				add_correlator_meson(rho_is, t, time_sources[srcnum], 1.0/(#time_sources));
				
				-- minus o Do

				t = q_gupta[i]:Re_dot(Do_gupta, "timeslices");
				add_correlator_meson(rho_is, t, time_sources[srcnum], -1.0/(#time_sources));
			end
			
			-- Next, y shifted. This lets us grab:
			-- 1. rho_ij: gamma_3 x taste_1 taste_4 taste_5
			-- 2. rho_i5: gamma_3 gamma_4 x taste_1 taste_5
			
			-- Odd soln.
			Do_gupta:symshift(o_gupta, unitg, 2);
			temp1:set(Do_gupta);
			temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,time_sources[srcnum]})
			temp2:set(Do_gupta);
			Do_gupta:combine({temp1, temp2}, {0.5,0.5});
			
			-- Do the same to the even solution.
			Dq_gupta:symshift(q_gupta[i], unitg, 2);
			temp1:set(Dq_gupta);
			temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,time_sources[srcnum]})
			temp2:set(Dq_gupta);
			Dq_gupta:combine({temp1, temp2}, {0.5,0.5});
			
			-- The combinations look familiar!
			-- First, rho_ij. This requires no phasing, just q Dq - o Do.

			do 
				-- q Dq
				t = q_gupta[i]:Re_dot(Dq_gupta, "timeslices");
				add_correlator_meson(rho_ij, t, time_sources[srcnum], 1.0/(#time_sources));
				
				-- minus o Do

				t = o_gupta:Re_dot(Do_gupta, "timeslices");
				add_correlator_meson(rho_ij, t, time_sources[srcnum], -1.0/(#time_sources));
			end

			-- Next, rho_is. This requires rephasing by a factor of (-1)^(x+y+z),
			-- then we take o Dq - q Do.
			
			Dq_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]});
			Do_gupta:rephase(make_phase_term(1,1,1,0), {0,0,0,time_sources[srcnum]});
			do 
				-- q Dq
				t = o_gupta:Re_dot(Dq_gupta, "timeslices");
				add_correlator_meson(rho_i5, t, time_sources[srcnum], 1.0/(#time_sources));
				
				-- minus o Do

				t = q_gupta[i]:Re_dot(Do_gupta, "timeslices");
				add_correlator_meson(rho_i5, t, time_sources[srcnum], -1.0/(#time_sources));
			end

			printf("End color.\n", i);
		end -- color
	
	-- Now that we're done with all the mesons, we can use the even and q
    -- solutions we saved to build Nucleons. For nucleons, one just sums over
    -- the primary site in the hypercube. We're going to reuse the trick we
    -- used to zero out sites with z%2==1 to zero out all sites with x%2 or
    -- y%2 == 1 as well!
	-- Note that we do NOT do a spin projection. 
	
		do -- Consistency check on nucleon
		
			-- Perform the zeroing.
			for i=1,3 do
				-- odd x
				temp1:set(even_soln[i]);
				temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,time_sources[srcnum]})
				temp2:set(even_soln[i]); -- reuse variable
				even_soln[i]:combine({temp2, temp1}, {0.5,0.5}) 

				-- odd y
				temp1:set(even_soln[i]);
				temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,time_sources[srcnum]})
				temp2:set(even_soln[i]); -- reuse variable
				even_soln[i]:combine({temp2, temp1}, {0.5,0.5}) 

				-- odd z
				temp1:set(even_soln[i]);
				temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,time_sources[srcnum]})
				temp2:set(even_soln[i]); -- reuse variable
				even_soln[i]:combine({temp2, temp1}, {0.5,0.5}) 
			end
		end
		

		-- Now that we've zeroed out non-primary sites, we perform the epsilon
		-- contraction. This takes the determinant of a matrix where each column
		-- corresponds to one of the three color vectors. 
		t = even_soln[1]:epsContract({even_soln[2],even_soln[3]}, "timeslices");
		add_correlator_baryon(nucleon_ck, t, time_sources[srcnum], 1.0/(#time_sources));
		
		
		-- Now do it for real!
		
		
		do -- Delta!
			-- An operator which measures the Delta baryon can be found in the end of
			-- Golterman's Lattice baryon paper. The method here is the operator defined
			-- by equation 6.3. 

			-- We need symmetric shifts in all directions for all colors.
			-- We also need to zero out all non-primary lattice sites. 
			for i=1,3 do -- For all colors
				for j=1,3 do -- For all directions

					-- Perform the symmetric shift.
					Dq_gupta_all[i][j]:symshift(q_gupta[i], unitg, j);

					-- Zero out non-primary sites.
					temp1:set(Dq_gupta_all[i][j]);
					temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,time_sources[srcnum]})
					temp2:set(Dq_gupta_all[i][j]); 
					Dq_gupta_all[i][j]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd x.

					temp1:set(Dq_gupta_all[i][j]);
					temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,time_sources[srcnum]})
					temp2:set(Dq_gupta_all[i][j]);
					Dq_gupta_all[i][j]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd y.

					temp1:set(Dq_gupta_all[i][j]);
					temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,time_sources[srcnum]})
					temp2:set(Dq_gupta_all[i][j]); 
					Dq_gupta_all[i][j]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd z.
				end
			end

			-- While there's slicker ways to do this, I explicitly construct
			-- an antisymmetric tensor. 
			eps_symbol = {};
			for i=1,3 do
				eps_symbol[i] = {};
				for j=1,3 do
					eps_symbol[i][j] = {};
					for k=1,3 do
						eps_symbol[i][j][k] = 0;
					end
				end
			end
			eps_symbol[1][2][3] = 1;
			eps_symbol[2][3][1] = 1;
			eps_symbol[3][1][2] = 1;
			eps_symbol[2][1][3] = -1;
			eps_symbol[3][2][1] = -1;
			eps_symbol[1][3][2] = -1;
			t = {}; 
			for k=1,3 do -- All elements
				for l = 1,3 do -- of the
					for m = 1,3 do -- epsilon symbol (color)
						-- Don't waste time with zero elements!
						if not (eps_symbol[k][l][m] == 0) then 
							-- Construct the proper epsilon contraction.
							t = Dq_gupta_all[k][1]:epsContract({Dq_gupta_all[l][2],Dq_gupta_all[m][3]}, "timeslices");
							
							add_correlator_baryon(delta, t, time_sources[srcnum], eps_symbol[k][l][m]/(#time_sources));
						end
					end
				end
			end
		end
		
		
		do -- Nucleon check! This is similar to how local meson states have a check.
			-- Zeroing.
			for i=1,3 do
				temp1:set(q_gupta[i]);
				temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,time_sources[srcnum]})
				temp2:set(q_gupta[i]); 
				q_gupta[i]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd x.

				temp1:set(q_gupta[i]);
				temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,time_sources[srcnum]})
				temp2:set(q_gupta[i]); 
				q_gupta[i]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd y.

				temp1:set(q_gupta[i]);
				temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,time_sources[srcnum]})
				temp2:set(q_gupta[i]); 
				q_gupta[i]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd z.
			end
		end

		-- Contract.
		t = q_gupta[1]:epsContract({q_gupta[2], q_gupta[3]}, "timeslices");
		add_correlator_baryon(nucleon, t, time_sources[srcnum], 1.0/(#time_sources));
		
		
		printf("End wall source.\n");
	end -- number of sources

end
		
		
-- Print it out!

--This is all to reproduce what MILC does.
printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
printf("SOURCE: RANDOM_WALL\n");
printf("SINKS: POINT_KAON_5 WALL_KAON_5\n");

for i = 1,#(pion_fpi_ptp) do
	printf("%i %.6e %f %.6e %f\n", i-1, pion_fpi_ptp[i], 0.0, pion_fpi_ptw[i], 0.0);
end

-- This is all to reproduce what MILC does.
printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
printf("SOURCE: EVEN_WALL\n");
printf("SINKS: PION_PS PION_SC\n");
--printf("Local Pions\n");
--printf("\tPion5\t\tPion5_4\n");

for i = 1,#(pion_ps_ck) do
	printf("%i %.6e %f %.6e %f\n", i-1, pion_ps_ck[i], 0.0, pion_4_ck[i], 0.0);
end


printf("ENDPROP\n");

printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
printf("SOURCE: EVEN_WALL\n");
printf("SINKS: NUCLEON\n");


--printf("Local Baryon\n");
--printf("\tNucleon\n");

for i = 1,#(nucleon_ck) do
printf("%i %.6e %f\n", i-1, nucleon_ck[i].r, 0.0);
end

printf("ENDPROP\n");

printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
printf("SOURCE: EVENANDODD_WALL\n");
printf("SINKS: PION_PS PION_SC PION_i5 PION_ij\n");


--printf("\nNonlocal Pions + Check\n")
--printf("\tPion5\tPion5_4\tPion_i5\tPion_ij\n");

for i = 1,#(p5) do
	printf("%i %.6e %f %.6e %f %.6e %f %.6e %f\n", i-1, p5[i], 0.0, p5_g4[i], 0.0, pion_i5[i], 0.0, pion_ij[i], 0.0);
end

printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
printf("SOURCE: EVENANDODD_WALL\n");
printf("SINKS: RHO_0 RHO_is RHO_ij RHO_i5\n");

for i = 1,#(rho_0) do
	printf("%i %.6e %f %.6e %f %.6e %f %.6e %f\n", i-1, rho_0[i], 0.0, rho_is[i], 0.0, rho_ij[i], 0.0, rho_i5[i], 0.0);
end

printf("ENDPROP\n");

printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n",mass, mass);
printf("SOURCE: EVENANDODD_WALL\n");
printf("SINKS: NUCLEON DELTA\n");


--printf("Nonlocal Baryons\n");
--printf("\tNucleon\tDelta\n");

for i = 1,#(nucleon) do
	printf("%i %.6e %f %.6e %f\n", i-1, nucleon[i].r, 0.0, delta[i].r, 0.0);
end

printf("ENDPROP\n");

totaltime = qopqdp.dtime() - totaltime;
printf("Total time: %f seconds.\n", totaltime);
