-- ESW 2016-03-16
-- Example for how to create a warm field. Should be integrated into 
-- gaugeact.lua, but this suffices as a first pass.
-- Some discussion copied from below:

-- To get a "warm" configuration, we:
-- 0. Start from a unit gauge.
-- 1. Generate a random traceless antihermitian matrix, 'force'. (a.k.a. a random gauge force.)
-- 2. Multiply the unit gauge by exp(wt*force), where 'wt' is a weight. The larger
--       the weight is, the hotter the configuration is.

require 'common'
require 'smear'
require 'run'
require 'topo'

-- lattice size
local nx = nx or 8;
local nt = nt or 16;

-- Prepare lattice size.
local latsize = {nx,nx,nx,nt};
local vol = 1;
local spatvol = nx*nx*nx;
local seed = seed or os.time();
printf("latsize = ");
for k,v in ipairs(latsize) do
	vol = vol*v;
	printf(" %i",v);
end
printf("\nvolume = %i\n", vol)
printf("seed = %i\n", seed)

--set up qopqdp
L = qopqdp.lattice(latsize) --need L to make force.
qopqdp.profile(profile or 0);
qopqdp.verbosity(0);
qopqdp.seed(seed);

-- start timer
totaltime = qopqdp.dtime();

-- Print some basic information about the configuration.
function getplaq(g)
  local ps,pt = g:action{plaq=1}
  local lat = qopqdp.lattice()
  local nd,vol = #lat,1
  for i=1,nd do vol=vol*lat[i] end
  local s = 0.25*nd*(nd-1)*vol
  printf("plaq ss: %+.8e  st: %+.8e  tot: %+.8e\n", ps/s, pt/s, 0.5*(ps+pt)/s)
end

-- load the gauge field
g= qopqdp.gauge();
g:unit();


-- To get a "warm" configuration, we:
-- 0. Start from a unit gauge.
-- 1. Generate a random traceless antihermitian matrix, 'force'. (a.k.a. a random gauge force.)
-- 2. Multiply the unit gauge by exp(wt*force), where 'wt' is a weight. The larger
--       the weight is, the hotter the configuration is.

-- Create a force.
--force = L:force();
--force:randomTAH(); -- Generate a random traceless antihermitian matrix.
--printf("Norm2 of force: %15.20e\n", force:norm2())


-- Look at the unit gauge.
printf("wt %+.8e ", 0); getplaq(g);

-- Loop over various weights 'wt'.
weights = {1e-5, 1e-4, 1e-3,
           1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2,
           1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1, 1};

for i=1,#weights do
  --g:unit();                      -- Reset to the unit gauge.
  --g:update(force, {weights[i]}); -- Multiply the unit gauge by exp(wt*force).
  weakfield(g, weights[i])
  printf("wt %+.8e ", weights[i]); -- Print weight.
  getplaq(g);                      -- Print plaquettes.
end

-- As a sanity check, I made sure unitarity is preserved after updating.
-- Worked like a charm unsurprisingly, but to be thorough, the check code is here:

--[[
do
  local devavg,devmax = g:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  g:makeSU()
  devavg,devmax = g:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end
]]---

