package.path = arg[0]:gsub("[^/]*.lua","?.lua") .. ";./hmc/?.lua;" .. package.path
require 'common'
require 'gaugeact'


local nx = nx or 4
local nt = nt or 8
local xi0 = xi0 or 1
local nu  = nu or 1

_G.aniso = nu/xi0

local aniso = aniso or 1
latsize = { nx, nx, nx, nt }
prec = 1
restart = 500
resid = 1e-6
local mass = mass or 0

qopqdp.lattice(latsize)
qopqdp.profile(profile or 0)
qopqdp.verbosity(0)


function getplaq(g)
  local ps,pt = g:plaq()
  printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", ps, pt, 0.5*(ps+pt))
end

g = qopqdp.gauge()
if inlat then g:load(inlat)
else
  seed = seed or os.time()
  qopqdp.seed(seed)
  --g:unit()
  g:random()
end

w = qopqdp.wilson()
w:printcoeffs()
w:set(g, prec)

mgIsSet = true
local kappa = 0.5/(mass+1.0+3.0*aniso)
printf("mass = %f, kappa = %f\n",mass, kappa)


-------------
-- accounting
-- ----------
local tmin = 9999999
local tmax = 0
local itsmin = 99999999
local itsmax = 0
local max_scale = 1.0
local max_setup_res = 1.0
local max_npost = 10
local max_npre = 5
local max_setup_change_fac = 1.0
local max_cres = 1.0


local block = {2,2,2,2}
local lat = {}
for i=1,#latsize do lat[i] = latsize[i]/block[i] end
w:mgSet({[-1] = { nlevels=1 }})

for change_fac = 0.1, max_setup_change_fac, 0.1 do
for scale=0, max_scale, 0.1 do
for cres=0.1, max_cres, 0.1 do
for setup_res=0.1, max_setup_res, 0.1 do
for npre=0, max_npre do 
for npost=1, max_npost do 
     mgpars =  {
	        [-2] = { verbose=0 },
	        [-1] = { verbose=1, kappa=kappa, kappanv=kappa, itmax=100, ngcr=8 },
	         [0] = { lattice=lat, nvecs=24,
	                 setup_res=setup_res, setup_maxit=100, setup_change_fac=change_fac,
	                 npre=npre, npost=npost, scale=scale, cres=cres, itmax=100, ngcr=8 }
	       }

w:mgSet(mgpars)
w:mgSetup()

opts = { prec=prec, restart=restart }
src = w:quark()
dest = {}

local color = 0
local spin = 0
src:zero()
-- point({coord},color,spin,re,im)
src:point({0,0,0,0},color,spin,1,0)
printf("src norm2: %g\n", src:norm2())
src:smearGauss(g, 0.9, 10);
printf("src norm2: %g\n", src:norm2())
dest[#dest+1] = w:quark()

t0 = qopqdp.dtime()
w:solve(dest[#dest], src, mass, resid, "all", opts)
dt = qopqdp.dtime() - t0
mf = 1e-6 * w:flops() / dt
printf("its: %g  secs: %g  Mflops: %g\n", w:its(), dt, mf)
for i=1,#dest do
  printf("%i norm2: %g\n", i, dest[i]:norm2())
  dest[i]:zero()
end
printf("\n\n")

if dt < tmin then 
      	tmin = dt 
        bestt_mgpars = mgpars
end

if w:its() < itsmin then
	itsmin = w:its()
	bestits_mgpars = mgpars
end

end -- npost
end -- npre
end -- setup_res
end -- cres
end -- scale
end -- change_fac

function print_mgpars (...) 
    for i, v in pairs(...) do 
	    for key, value in pairs(v) do 
		    print(i, key, value)
            end
    end
end

printf("best time = %f secs, with mg parameters = \n", tmin)
print_mgpars(bestt_mgpars)

printf("best its = %g,       with mg parameters = \n", itsmin)
print_mgpars(bestits_mgpars)

