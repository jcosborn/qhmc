L = Lattice{4,4,4,8}
G = L:GaugeField{group="SU",nc=3}
G:Load("lattice")
GA = GaugeAction{kind="plaquette",beta=6,field=G}
M = G:Momentum()
I = Integrator{kind="leapfrog",action=GA,field=G,momentum=M,tau=1,nSteps=40}
E = Evolver{kind="HMC",integrator=I}
E:Run()
