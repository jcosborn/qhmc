nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)

function test(x)
  local sa = L:subset("all")
  local ra = f:norm2(sa)
  local se = L:subset("even")
  local re = f:norm2(se)
  local so = L:subset("odd")
  local ro = f:norm2(so)
  --print("a: "..tostring(ra).."  e: "..tostring(re).."  o: "..tostring(ro))
  printf("a: %i  e: %i  o: %i\n", ra, re, ro)
end

f = L:real()
f:zero()
test()

f:point({0,0,0,0},1)
test()

f:point({0,0,0,1},2)
test()

print(f:point({0,0,0,1}))
