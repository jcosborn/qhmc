-- ESW 10-10-2014
-- This standalone file measures the connected and
-- disconnected scalar using time, color, and space even-odd 
-- dilution. 

require 'common'
require 'smear'
require 'run'

-- This is a modified version of one of James' scripts for doing hmc generations.
-- It has been changed to the point that I'm not quite sure which one it's based off,
-- but in any case, it is designed in a way to only perform measurements, despite the
-- fact it looks like it's set up to do HMC generation.


-- Parameters we need to set.
--[[
local nx = 24 -- Lattice spatial dimension
local nt = 48 -- Lattice temporal dimension
local mass = 0.005 -- Light quark mass.
local inlat = "configurations/f4plus8l24t48b40m005m080_:CONFIGNUM:_scidac" -- Input configuration
local cg_prec = 1e-6 -- Max residual of CG.
local cg_max = 4000 -- Maximum number of CG steps.
--]]

local nx = nx or 4
local nt = nt or 8
local mass = mass or 0.1
--local inlat = "hyp1632f8b50m05a_coulomb.1002.lime"
local cg_prec = cg_prec or 1e-6
local cg_max = cg_max or 2000

-- Stuff for dilution!
local num_src = num_src or 2 -- Use 6 sources.
local color_dil = color_dil or qopqdp.defaultNc(); -- Set to 1 to not color dilute!
local source_type = source_type or "U1" -- other option is "Z2".

-- Must be set to not zero to use lowmemory dilution.
local time_dil = time_dil or nt -- set to 1 to not time dilute. 
local space_dil = space_dil or 2 -- set to 1 to not space even/odd dilute


-- These parameters matter!
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
printf("cg_prec = %g\n", cg_prec)
printf("num_src_src = %i\n", num_src)
printf("source_type = %s\n", source_type)

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

-- No need to gauge fix!

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

-- Now we can prepare to measure things!

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

-- A utility function. We have an index that goes from 1 to the
-- total number of inversions we need, and this function takes
-- that and spits out a noise source, color, timeslice, and even/odd
-- index.
function get_dilution_pattern_colorin(index, numstoch, numtimeslice, numcolor, numeo)
	local tmpnum = index-1 -- switch to zero indexing.
	local eo_tmp; -- 0 for even, 1 for odd
	local color_tmp; -- zero indexed color.
	local timeslice_tmp; -- zero indexed timeslice.
	local stoch_tmp; -- zero indexed source number.
	
	if numcolor == 1 then
	   color_tmp = -1;
	else
	   color_tmp = tmpnum % numcolor;
	   tmpnum = (tmpnum - color_tmp)/numcolor;
	end
	
	if numeo == 1 then
	   eo_tmp = -1;
	else
	   eo_tmp = tmpnum % numeo;
	   tmpnum = (tmpnum-eo_tmp)/numeo;
	end
	
	if numtimeslice == 1 then
	   timeslice_tmp = -1;
	else
	   timeslice_tmp = tmpnum % numtimeslice;
	   tmpnum = (tmpnum - timeslice_tmp)/numtimeslice;
	end
	
	stoch_tmp = tmpnum % numstoch;
	tmpnum = (tmpnum - stoch_tmp)/numstoch;
	
	if tmpnum ~= 0 then
	   printf("Inconsistent tmpnum? %i\n", tmpnum);
	end
	
	return {stoch = stoch_tmp, timeslice = timeslice_tmp, color = color_tmp, eo = eo_tmp};

end

-- We're now ready to measure! First, get 'num_src' sources.
-- eta_full = act.f:disc_get_sources(num_src, 0)

printf("Preparing to get %d %s noise sources.\n", num_src, source_type);

eta = {};

if (source_type == "U1") then
	for i=1,num_src do
		eta[i] = w:quark();
		eta[i]:randomU1();
	end
elseif (source_type == "Z2") then
	for i=1,num_src do
		eta[i] = w:quark();
		eta[i]:randomZ2();
	end
else -- Not a valid type.
	printf("Error: %s is not a valid noise source. U1 and Z2 are available.\n");
	os.exit();
end

printf("Obtained %d %s noise sources.\n", num_src, source_type);

-- Next, we actually perform the dilutions, inversions, and contractions.
--local retvals = act.f:disc_get_all_lowmem_fast(eta_full, time_dil, color_dil, space_dil, G, mass, cg_prec, {prec = prec, restart = cg_max});

-- We use a low memory trick to perform the connected part
-- of this calculation. 


-- Prepare some storage space!

-- All pbp's.
local pbps_unimp = {}
local pbps_isimp = {}
	
-- pbp per timeslice.
local pbps_unimp_t = {}
local pbps_isimp_t = {}
for i=1,num_src do
	pbps_unimp_t[i] = {};
	pbps_isimp_t[i] = {};
end

-- correlators
local corr_scalar = {};
	
-- Two quark objects---one to hold the source, one the prop.
-- We name these in the spirit of the dilution literature.
local eta_dil = w:quark();
eta_dil:zero();
local phi = w:quark();
phi:zero();

-- Find the total amount of dilution we're doing.
local total_dil = time_dil*color_dil*space_dil;
	
-- Time it!
local t1 = qopqdp.dtime();
	
-- Loop over all phi's.
for i=1,(num_src*total_dil) do -- right index!
	-- Prepare some memory.
	local op = 0;
	
	-- Get a dilution pattern.
	phi_dilution = get_dilution_pattern_colorin(i, num_src, time_dil, color_dil, space_dil);
	
	-- Time it.
	local t0 = qopqdp.dtime();
	
	-- Dilute and get the corresponding phi!
	eta_dil:dilute(eta[phi_dilution.stoch+1],phi_dilution.timeslice, phi_dilution.color, phi_dilution.eo)
	solve_printinfo(w, phi, eta_dil, mass, cg_prec, "all", {prec = prec, restart = cg_max});
		
	-- Now perform the right contraction!
	
	-- Now loop over all eta's. 
	-- ESW 10-23-2014 Optimizing trick
	--for j=1,(num_src*total_dil) do -- left index!
	for j=1,(num_src*color_dil*space_dil) do
	
		-- Figure out what dilution pattern this corresponds to. We're kind of hijacking this for the optimize change!
		eta_dilution = get_dilution_pattern_colorin(j, num_src, time_dil, color_dil, space_dil);
		
		-- ESW 10-23-2014: Fix up what we really mean.
		-- We looped over a smaller amount than we expected... so we re-interpret the outputs.
		-- Remember, we don't dilute in time anymore. We do the timeslice contraction to do it for us.
		eta_dilution.stoch = eta_dilution.stoch*time_dil + eta_dilution.timeslice;
		eta_dilution.timeslice = -1; -- Don't care about time dilution anymore!
		
		-- Get this diluted eta.
		-- We dilute inside now!
		-- eta_dil:dilute(eta[eta_dilution.stoch+1],eta_dilution.timeslice, eta_dilution.color, eta_dilution.eo)
		
		--if (i == j) then -- This contributes to pbp.
		if (phi_dilution.stoch == eta_dilution.stoch and phi_dilution.color == eta_dilution.color and phi_dilution.eo == eta_dilution.eo) then
			-- Dilute!
			-- Remember, the timeslice=-1 means we don't dilute there.
			eta_dil:dilute(eta[eta_dilution.stoch+1],eta_dilution.timeslice, eta_dilution.color, eta_dilution.eo)
			
			tmp_unimp = eta_dil:reDot(phi,"timeslices");
			tmp_isimp = phi:norm2("timeslices"); --"timeslices");
			pbps_unimp_t[phi_dilution.stoch+1][phi_dilution.timeslice+1] = tmp_unimp[phi_dilution.timeslice+1]/vol + (pbps_unimp_t[phi_dilution.stoch+1][phi_dilution.timeslice+1] or 0);
			pbps_unimp[phi_dilution.stoch+1] = tmp_unimp[phi_dilution.timeslice+1]/vol + (pbps_unimp[phi_dilution.stoch+1] or 0);
			-- Distribute it around for the improved.
			for t=1,#tmp_isimp do
				pbps_isimp_t[phi_dilution.stoch+1][t] = mass*tmp_isimp[t]/vol + (pbps_isimp_t[phi_dilution.stoch+1][t] or 0);
				pbps_isimp[phi_dilution.stoch+1] = mass*tmp_isimp[t]/vol + (pbps_isimp[phi_dilution.stoch+1] or 0)
			end
			--pbps_isimp_t[phi_dilution.stoch+1][phi_dilution.timeslice+1] = mass*tmp_isimp/vol + (pbps_isimp_t[phi_dilution.stoch+1][phi_dilution.timeslice+1] or 0);
			--pbps_isimp[phi_dilution.stoch+1] = mass*tmp_isimp/vol + (pbps_isimp[phi_dilution.stoch+1] or 0);
			
			-- And then I guess we don't need the others anymore? I can check this.
			
			-- So I don't think this is all relevant anymore?
			--[[for t=1,#tmp_unimp do
				-- Okay, so... I need think about this for a bit.
				pbps_unimp_t[phi_dilution.stoch+1][t] = tmp_unimp[t]/vol + (pbps_unimp_t[phi_dilution.stoch+1][t] or 0);
				pbps_unimp[phi_dilution.stoch+1] = tmp_unimp[t]/vol + (pbps_unimp[phi_dilution.stoch+1] or 0)
				pbps_isimp_t[phi_dilution.stoch+1][t] = mass*tmp_isimp[t]/vol + (pbps_isimp_t[phi_dilution.stoch+1][t] or 0);
				pbps_isimp[phi_dilution.stoch+1] = mass*tmp_isimp[t]/vol + (pbps_isimp[phi_dilution.stoch+1] or 0)
			end
			]]
			
		elseif (eta_dilution.stoch < phi_dilution.stoch) then -- Contributes to the connected part.
		
			-- Dilute!
			-- Remember, the timeslice=-1 means we don't dilute there.
			eta_dil:dilute(eta[eta_dilution.stoch+1],eta_dilution.timeslice, eta_dilution.color, eta_dilution.eo)
		
			-- Get the contractions! This is the real magic---at the expense of repeating
			-- doing dilution, we massively cut down on memory use.
			-- Get all the contractions!
			op = eta_dil:dot(phi, "timeslices");
			
			-- What's the magic?
			-- The timeslice sum is O_{(i,a,t')(j,b,t'')}(t). It's zero unless t = t'.
			-- Thus, this is O_{(i,a,t')(j,b,t'')}(t').
			-- The dagger of this is O_{(j,b,t'')(i,a,t')}(t'')*a*b*(-1)^(t'+t'')
			
			-- Thus, we have O_{(i,a,t')(j,b,t'')}(t') O_{(j,b,t'')(i,a,t')}(t'')
			
			-- Figure out the phase correction.
			phase = 1;
			if (eta_dilution.eo == 1) then -- if it's odd
				phase = -phase;
			end
			if (phi_dilution.eo == 1) then -- if it's odd
				phase = -phase;
			end
			
			-- Remember, these times are zero indexed. 
			
			for t=1,#op do
				local difference = (phi_dilution.timeslice-t+1+2*nt)%nt;
			
				-- Add in the magnitude.
				corr_scalar[difference+1] = phase*(op[t].r*op[t].r+op[t].i*op[t].i) + (corr_scalar[difference+1] or 0);
			end
		end
		
	end
	t0 = qopqdp.dtime() - t0;
	printf("Inversion %i, source %i, timeslice %i, space %i, color %i, total time: %f\n", i, phi_dilution.stoch, phi_dilution.timeslice, phi_dilution.eo, phi_dilution.color, t0);
	io.stdout:flush();
	
end

t1 = qopqdp.dtime() - t1;
printf("Total Invert+Dot Time: %f\n", t1);
io.stdout:flush();

-- All done with the calculation! Spit it out!

printf("BEGIN_PBP\n");
for i=1,#(pbps_unimp) do
   --printf("%i %.15e %.15e %.15e\n", i, pbps_other[i], pbps.unimp[i], pbps.isimp[i])
   printf("%i %.15e %.15e\n", i, pbps_unimp[i], pbps_isimp[i])
end
printf("END_PBP\n");
io.stdout:flush()

printf("BEGIN_PBPPART\n");
for i=1,#(pbps_unimp_t) do
   for t=1,#(pbps_unimp_t[i]) do
     printf("%i\t%i\t%.15e\t%.15e\n", i, t, pbps_unimp_t[i][t], pbps_isimp_t[i][t])
   end
end

printf("END_PBPPART\n");


printf("T\tScalar\n");
printf("BEGIN_SPECTRUM\n");

for i = 1,nt do
        printf("%i\t%.15e\n", i-1, corr_scalar[i]);
end

printf("END_SPECTRUM\n");

totaltime = qopqdp.dtime() - totaltime;

printf("Total time: %f seconds.\n", totaltime);


io.stdout:flush()

