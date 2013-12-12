require 'Lattice'
require 'Action'
require 'Evolver'

L = Lattice{4,4,4,8}
L:Seed(987654321)
G = L:GaugeField{group="SU",nc=3}
G:Set("unit")
--G:Load("lattice")
GA = Action{kind="gauge",style="plaquette",beta=6,field=G}
M = G:Momentum()
MA = Action{kind="momentum",momentum=M}
--MHB = Evolver{kind="hb",action=MA,field=M}
I = Evolver{kind="md",style="leapfrog",action=GA,
	    field=G,momentum=M,tau=1,nSteps=40}
--MC = Evolver{kind="sequence",evolvers={MHB,I}}
E = Evolver{kind="mc",markov=I,actions={MA,GA},fields={G}}
--E = Evolver{kind="metropolis",markov=I,actions={MA,GA},fields={G}}
printf("action: %g\n", GA:Action())
E:Run()
printf("action: %g\n", GA:Action())
myprint(E.oldActions,"\n")
myprint(E.newActions,"\n")
myprint(E.mcRand,"\n")
