require 'common'
require 'gaugeact'
require 'wilsonObservables'

--trace(true)
nx = nx or 4
nt = nt or 8
latsize = { nx, nx, nx, nt }
--aniso = 1/2.38
aniso = 1
mass = mass or -0.4125
--rsmear = 0.6
--nsmear = 30
rsmear = 0
nsmear = 0
prec = 1
restart = 500
resid = 1e-12
opts = { prec=prec, restart=restart }

--qopqdp.defaultNc(4); fn = nil
Nc = Nc or qopqdp.defaultNc()
qopqdp.defaultNc(Nc)
Ns = 4

L = qopqdp.lattice(latsize)
qopqdp.profile(profile or 0)
qopqdp.verbosity(0)

seed = seed or os.time()
qopqdp.seed(seed)

--wr = qopqdp.writer("test.out", "<metadata/>")

function getplaq(g)
  local ps,pt = g:action{plaq=1}
  local lat = qopqdp.lattice()
  local nd,vol = #lat,1
  for i=1,nd do vol=vol*lat[i] end
  local s = 0.25*nd*(nd-1)*vol
  printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", ps/s, pt/s, 0.5*(ps+pt)/s)
end

g0 = qopqdp.gauge()
if fn then g0:load(fn)
else
  --g0:unit()
  g0:random()
end

getplaq(g0)
-- coulomb(j_decay, error, max iterations, overrelaxation param)
-- note that here 0->x, ..., 3->t
--g0:coulomb(3, 1e-7, 1000, 1.2)
--getplaq(g0)

-- background EM field
function setE(g, q)
  -- U_z = exp(-i 2 pi q t/(nt*nx)) U_z
  g(3):momentum({0,0,0,-q/nx},0,0,0,1)
  -- U_t = exp(i 2 pi q z delta(t,nt-1)/nx) U_t
  g(4):momentum({0,0,q,0},0,0,0,1,"timeslice"..(nt-1))
end
qu = 2
qd = -1
qn = 1
g = g0:copy()
setE(g, qn*qu)

w = qopqdp.wilson()
w:printcoeffs()
w:set(g, {aniso=aniso}, prec)

function mgSetup()
  local kappa = 0.5/(mass+1.0+3.0*aniso)
  printf("mass = %f, kappa = %f\n", mass, kappa)

  local block = {2,2,2,2}
  --local block = {3,3,3,4}
  local lat = {}
  for i=1,#latsize do lat[i] = latsize[i]/block[i] end
  w:mgSet({[-1] = { nlevels=1 }})
  w:mgSet(
    {
      [-2] = { verbose=0 },
      [-1] = { verbose=1, kappa=kappa, kappanv=kappa, itmax=100, ngcr=8 },
      [0] = { lattice=lat, nvecs=24,
	      setup_res=0.4, setup_maxit=100, setup_change_fac=0.5,
	      npre=0, npost=4, scale=1, cres=0.2, itmax=100, ngcr=8 }
    }
  )
  w:mgSetup()
end

--mgSetup()

function smear(f)
  f:smearGauss(g, 4, rsmear, nsmear)
end
--smearsrc = smear
--smeardest = smear

function spectrum()
  local pt = {0,0,0,0}
  --local pt = randomPoint(L)
  printf("src point: %i %i %i %i\n", pt[1], pt[2], pt[3], pt[4])
  dest2 = pointProp(L, w, pt, smearsrc, smeardest, mass, resid, opts)

  --wr:write(dest, "<field metadata/>")
  --wr:prop(dest, "<field metadata/>")
  --for i=1,#dest do
  --printf("%i norm2: %g\n", i, dest[i]:norm2())
  --dest[i]:zero()
  --end

  --rd,md = qopqdp.reader("test.out")
  --printf("%s\n", md)
  --md = rd:read(dest)
  --printf("%s\n", md)
  --for i=1,#dest do
  --printf("%i norm2: %g\n", i, dest[i]:norm2())
  --end
  --[[
  local d11 = dest2[1][1]
  for i=0,Nc-1 do
    for j=0,Nc-1 do
      local s = d11:point(pt,i,j)
      printf("%i %i\t%s\n", i, j, s)
    end
  end
  os.exit(1)
  ]]

  --mesons = wilsonMesons(dest)
  --mesons = wilsonMesons2(dest2)
  mesons = wilsonMesons3(dest2)

  --[[
    for g = 0,#mesons do
    myprint("meson["..tostring(g).."] = ",mesons[g],"\n")
    end
  --]]

  printf("-= mesons =-\n")
  for t=0,latsize[4]-1 do
    local t1 = (pt[4]+t)%latsize[4] + 1
    for g = 0,#mesons do
      printf("%i\t%i\t%i\t%g\n", t, g, g, mesons[g][t1])
    end
  end

  baryons = wilsonBaryons3(dest2)

  if baryons then
    printf("-= baryons =-\n")
    for t=0,latsize[4]-1 do
      local t1 = (pt[4]+t)%latsize[4] + 1
      printf("%i", t)
      for g = 0,#baryons do
	printf("\t%g\t%g", baryons[g][t1].r, baryons[g][t1].i)
      end
      printf("\n")
    end
  end
end

function scalar()
  -- random momentum source
  --local src, mom, cmom, smom = randMomSource(L)
  --printf("src momentum: %i %i %i %i : %i %i\n", mom[1], mom[2], mom[3], mom[4], cmom, smom)
  --src:point(mom, cmom, smom, 1)
  local src = L:diracFermion()
  src:randomU1()

  local tsrc = L:diracFermion()
  local dest = L:diracFermion()
  local iv3 = 1/math.sqrt(L(1)*L(2)*L(3))
  local c = {}
  for t=1,L(4) do -- source timeslice
    local sub = "timeslice"..(t-1)
    tsrc:zero()
    tsrc:set(src, sub)
    w:solve(dest, tsrc, mass, resid, "all", opts)
    c[t] = src:dot(dest, "timeslices")
    for i=1,L(4) do c[t][i] = iv3 * c[t][i] end
  end

  printf("current:\n")
  for t0=1,L(4) do -- source timeslice
    printf("%i\t%.10g\t%.10g\n", t0-1, c[t0][t0].r, c[t0][t0].i)
  end

  printf("correlator:\n")
  for t0=1,L(4) do -- source timeslice
    for tl=0,L(4)-1 do -- distance
      local t1 = (t0+tl-1)%L(4) + 1
      local conn = c[t0][t1] * c[t1][t0]
      printf("%i\t%i\t%.10g\t%.10g\n", t0-1, tl, conn.r, conn.i)
    end
  end
end

spectrum()
printf("-= scalar =-\n")
scalar()
