require 'Lattice'
require 'GaugeAction'
require 'Integrator'

profile(0)
verbosity(0)

nx = nx or 4
nt = nt or 8
Nc = Nc or 3
beta = beta or 6
seed = seed or os.time()

L = Lattice{nx,nx,nx,nt}
L:Set{defaultNc=Nc}
L:Seed(seed)

G = L:GaugeField()
--printf("%s", G)
G:set("unit")

GA = GaugeAction{type="plaquette",beta=6,field=G}
--printf(GA._type)
--GA.printForce = true
--myprint(GA,"\n")
printf("%s", GA)

F = G:Force()

I = Integrator{type="leapfrog",action=GA,field=G,force=F,tau=1,nSteps=40}
I.tostringRecurse = true
I.tostringPrintAll = true
local function stepHook(t, fs, int)
  local a = GA:action()
  local b = 0.5*F.force:norm2()
  printf("%g: %g\t%g\t%g\t%g\n", t, fs.eps[1], a, b, a+b)
end
--I.stepHook = stepHook
print(I)

I:reset()
F.force:random()
f2 = F.force:norm2()
--f2 = 0.5*f2 - 16*L.vol
f2 = 0.5*f2
printf("f2: %g\n", f2)

ga = GA:action()
printf("ga: %g\n", ga)

I:run()

ga = GA:action()
printf("ga: %g\n", ga)


--E = Evolver()
