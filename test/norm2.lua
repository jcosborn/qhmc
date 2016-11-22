-- This example is part of the regression test suite.  Functions starting
-- with TEST are part of the test framework and can be ignored.

nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)
L:seed(98765321)

TESTZEROTOL(1e-14)
TESTON()

ltypes = { "real", "complex", "colorVector", "diracFermion", "colorMatrix" }

for k,t in ipairs(ltypes) do
  f = L[t](L)
  f:random()
  printf("%g\n", f:norm2())
end

r = f:lnorm2()
n = r:sum()
printf("%g\n", n)

r:zero()
f:lnorm2(r)
n = r:sum()
printf("%g\n", n)

c = f:ldot(f)
n = c:sum()
printf("%g\t%g\n", n.r, n.i)

c:zero()
f:ldot(f, c)
n = c:sum()
printf("%g\t%g\n", n.r, n.i)

a = f:norm2("timeslices")
myprint(a,"\n")

subs = L:subset("timeslices")
a = f:norm2(subs)
myprint(a,"\n")

subs = {L:subset("even"),L:subset("odd")}
a = f:norm2(subs)
myprint(a,"\n")

TESTOFF()
