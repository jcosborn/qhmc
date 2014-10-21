nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)

TESTON()

function test()
  local sa = L:subset("all")
  local ra = f1:dot(f2,sa)
  local se = L:subset("even")
  local re = f1:dot(f2,se)
  local so = L:subset("odd")
  local ro = f1:dot(f2,so)
  --print("a: "..tostring(ra).."  e: "..tostring(re).."  o: "..tostring(ro))
  printf("a: %s  e: %s  o: %s\n", tostring(ra), tostring(re), tostring(ro))
end

w = qopqdp.wilson()
f1 = w:quark()
f2 = w:quark()
f1:zero()
f2:zero()
test()

-- p coords, p color, p spin, factor, field factor
f1:momentum({0,0,0,0},0,0,1,0)
f2:momentum({0,0,0,0},0,0,1,0)
test()

--f2:momentum({1,2,3,1},1,2,1,0)
--test()

nc = f1:nc()
ns = 4
t = 2*math.pi*((6/nx)+(1/nt)+(1/nc)+(2/ns))
s = complex(math.cos(t),-math.sin(t))
f1:momentum({1,2,3,1},1,2,1,0)
f2:momentum({1,2,3,1},1,2,s,0)
test()

f1:zero()
f2:zero()
f1:momentum({1,2,3,1},1,2,1,0,"timeslice2")
f2:momentum({1,2,3,1},1,2,s,0,"timeslice2")
test()

TESTOFF()
