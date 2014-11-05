require 'common'

nx = nx or 4
nt = nt or 4
beta = beta or 5
mass = mass or 0.05
tau = tau or 1
nsteps = nsteps or 100
ntraj = ntraj or 10
first = first or -1

latsize = { nx, nx, nx, nt }
L = qcd.lattice(latsize)
setup(L)
printf("%s\n", tostring(L))
printf("volume = %i\n", L:volume())
printf("seed = %i\n", seed)
printf("beta = %g\n", beta)
printf("mass = %g\n", mass)
printf("tau = %g\n", tau)
printf("nsteps = %g\n", nsteps)
printf("ntraj = %g\n", ntraj)

basefn = string.format("l4x%i%ib%03.0fm%02.0f.", nx, nt, 100*beta, 100*mass)
last = ntraj * tau
if first >= 0 then
  infn = basefn .. tostring(first)
  U = loadGauge(infn)
  last = last + first
else
  U = {}
  for i=1,#L do
    U[i] = L:ColorMatrix(complex(1,0))
  end
end
outfn = basefn .. tostring(last)

act = {}
act.g = gaugeact(L, beta)
act.f = stagact(L, mass)
--act.f.invstats = true
--act:checkH()
--act.g.printforce = true
--act.f.printforce = true
act.f.nPhi = 1
--act.f.dynPhi = true

act.f1 = act.f
act.f.f1 = act.f1

local fields = {}
function fields.save(f)
  f.USave = f.U
  f.phi = act.f1:newphi(f.U)
  f.H = act.g:newH()
  f.Hphi = act.f1:newH()
end
function fields.action(f)
  local Sgq = act.g:gauge(f.U)
  local Sgp = act.g:normH(f.H)
  local Sfq = act.f1:action(f.U, f.phi)
  local Sfp = act.f1:normH(f.Hphi)
  local S = Sgq + Sgp + Sfq  + Sfp
  printf("Sgq: %-8.6g  Sgp: %-8.6g", Sgq, Sgp)
  if(Sfq~=0) then printf("  Sfq: %-8.6g", Sfq) end
  if(Sfp~=0) then printf("  Sfp: %-8.6g", Sfp) end
  printf("\n")
  return S
end
function fields.accept(f)
end
function fields.reject(f)
  f.U = f.USave
end
function fields.n(f)
  return 1
end
function fields.updateField(f, i, eps)
  f.U = act.g:updateU(f.U, f.H, eps)
  f.phi = act.f1:updatephi(f.phi, f.Hphi, eps)
end
function fields.updateMomentum(f, i, eps)
  f.H = act.g:updateH(f.H, f.U, eps)
  f.H = act.f1:updateH(f.H, f.U, f.phi, eps, f.Up)
  f.Hphi = act.f1:updateHphi(f.Hphi, f.U, f.phi, eps)
end

fields.U = U

eps = tau/nsteps
hmcparams = {}
hmcparams.eps = eps
hmcparams.nsteps = nsteps
--hmcparams.traceS = true
--qcd.defaults{qdpProfcontrol=1}

function measure(U,traj)
  ps,pt = plaq(U)
  printf("plaq ss: %g  st: %g  tot: %g\n", ps, pt, 0.5*(ps+pt))
  io.open("plaq.dat","a")
  io.write("%i %g %g %g\n", ps, pt, 0.5*(ps+pt))

  plp = ploop(U)
  printf("ploop: %s\n", tostring(plp))

  ccond = act.f:chiralcond(U, 4)
  printf("condensate:")
  for i=1,#ccond do
    printf(" %g", ccond[i])
  end
  printf("\n")

  io.stdout:flush()
end

totaltime = clock()

io.open("plaq.dat","w")

measure(U,0)
for traj=1,ntraj do
  t0 = clock()
  hmcstep(fields, hmcparams)
  act.g:checkSU(fields.U)
  t0 = clock() - t0

  printf("traj %i secs: %g\n", traj, t0)
  act.g:updateStats()
  printf("GF secs: %8.3f  mflops: %5.0f\n", act.g.GFtime, act.g.GFmflops)
  act.f:updateStats()
  printf("FF secs: %8.3f  mflops: %5.0f\n", act.f.FFtime, act.f.FFmflops)
  printf("LL secs: %8.3f  mflops: %5.0f\n", act.f.LLtime, act.f.LLmflops)
  printf("CG secs: %8.3f  mflops: %5.0f  iters: %i  its/step: %i\n",
	 act.f.CGtime, act.f.CGmflops, act.f.CGiters,
	 act.f.CGiters/nsteps)
  ot = t0 - act.g.GFtime - act.f.FFtime - act.f.LLtime - act.f.CGtime
  printf("?? secs: %8.3f\n", ot)
  act.g:clearStats()
  act.f:clearStats()

  measure(fields.U,traj)
end

if outfn then saveGauge(U, outfn) end

totaltime = clock() - totaltime
printf("total time: %g seconds\n", totaltime)
