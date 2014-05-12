nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)

f = L:real()
f:zero()
f:unit()

sa = L:subset("all")
r = f:norm2(sa)
print(r)

se = L:subset("even")
r = f:norm2(se)
print(r)

so = L:subset("odd")
r = f:norm2(so)
print(r)

f:zero(se)
print(f:norm2(se))
print(f:norm2(so))
print(f:norm2(sa))
