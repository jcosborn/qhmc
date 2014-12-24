require 'Lattice'
require 'Action'
require 'Evolver'

qopqdp.profile(0)
qopqdp.verbosity(0)
--trace(1)

nx = 4
nt = 8
--L = Lattice{nx,nx,nx,nt}
L = Lattice{nx,nx,nt}
--printf(L._type)
--myprint(L,"\n")
L:Seed(987654321)

G = L:GaugeField{nc=2}
--printf(G._type)
--myprint(G,"\n")
printf("%s", G)

G:Set("unit")
--G:Save("test.lat")

GA = Action{kind="gauge",style="plaquette",beta=6,field=G}
printf(GA._type)
--GA.printForce = true
--myprint(GA,"\n")
printf("%s", GA)

ga = GA:Action()
printf("ga: %g\n", ga)
oldga0 = GA.act0
GA.act0 = 0
ga = GA:Action()
printf("ga: %g\n", ga)
GA.act0 = oldga0

M = G:Momentum()
M.field:zero()
G.field:force(M.field, {plaq=1}, 6, 1)


MD = Evolver{kind="md",style="leapfrog",action=GA,field=G,momentum=M,
	     tau=1,nSteps=100}
MD.tostringRecurse = true
MD.tostringPrintAll = true
local function stepHook(t, fs, int)
  local a = GA:action()
  local b = 0.5*M.force:norm2()
  printf("%g: %g\t%g\t%g\t%g\n", t, fs.eps[1], a, b, a+b)
end
--I.stepHook = stepHook
print(MD)
print(M)

Am = Action{kind="momentum",momentum=M}
Am:Refresh()
am = Am:Action()
printf("am: %g\n", am)
ga = GA:Action()
printf("ga: %g\n", ga)
printf("act: %g\n", am+ga)

MD:Reset()
MD:Run()

am = Am:Action()
printf("am: %g\n", am)
ga = GA:Action()
printf("ga: %g\n", ga)
printf("act: %g\n", am+ga)


MC = Evolver{kind="mc",actions={Am,GA},fields={G},markov=MD}
MC:Reset()
MC:Run()

myprint(MC.oldActions)
myprint(MC.newActions)
