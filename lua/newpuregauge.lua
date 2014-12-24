require 'Lattice'
require 'Action'
require 'Evolver'

qopqdp.profile(0)
qopqdp.verbosity(0)

nx = nx or 4
nt = nt or 8
Nc = Nc or 3
beta = beta or 6
seed = seed or os.time()

L = Lattice{nx,nx,nx,nt}
L:Set{defaultNc=Nc}
L:Seed(seed)

G = L:GaugeField{group="SU",nc=Nc}
--printf("%s", G)
G:Set("unit")

GA = Action{kind="gauge",style="plaquette",beta=6,field=G}
--printf(GA._type)
--GA.printForce = true
--myprint(GA,"\n")
printf("%s", GA)

M = G:Momentum()

MD = Evolver{kind="md",style="leapfrog",action=GA,field=G,momentum=M,
	     tau=0.1,nSteps=20}
MD.tostringRecurse = true
MD.tostringPrintAll = true
local function stepHook(t, fs, int)
  local a = GA:Action()
  local b = 0.5*M.field:norm2()
  printf("%6g: %6g\t%8g\t%8g\t%8g\n", t, fs and fs.eps[1] or 0, a, b, a+b)
end
MD.stepHook = stepHook
print(MD)

MD:Reset()
M:Random()
f2 = M.field:norm2()
--f2 = 0.5*f2 - 16*L.vol
f2 = 0.5*f2
printf("f2: %g\n", f2)

ga = GA:Action()
printf("ga: %g\n", ga)

MD:Run()

ga = GA:Action()
printf("ga: %g\n", ga)


--E = Evolver()
