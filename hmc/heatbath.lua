require 'Lattice'
require 'Action'
require 'Evolver'
require 'topo'

--profile(1)
--verbosity(1)
--trace(1)

nx = 4
nt = 8
beta = 6

--L = Lattice{nx,nx,nt}
L = Lattice{nx,nx,nx,nt}
L:Seed(987654321)

G = L:GaugeField{nc=3}
G:Set("unit")

GA = Action{kind="gauge",style="plaquette",beta=beta,field=G}

E = Evolver{kind="heatbath",actions={GA},
	    nRepetitions=10,nHeatbath=1,nOverrelax=1}

function plaq()
  local a = GA:Action()
  local b = a/beta
  local p = GA.act0 - b
  local s = 1/GA.act0
  return p*s
end

function topo()
  local tr,ti = pathDo(G.field, #L.latsize, paths0, coeffs0)
  return tr
end

printf("initial plaq: %g\n", plaq())

for i=1,10 do
  E:Run()
  printf("plaq: %g\n", plaq())
  printf("topo: %g\n", topo())
end
