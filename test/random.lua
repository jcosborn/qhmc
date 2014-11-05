nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)

rs = L:rstate()

rs:seed(987654321)
x1 = rs:globalRand()
x2 = rs:globalRand()
x3 = rs:globalRand()

rs:seed(987654321)
r = L:real()
r:randomUniform(rs)
y1 = r:site({0,0,0,0})
y2 = rs:globalRand()
r:randomUniform(rs)
y3 = r:site({0,0,0,0})

TESTOUT("x1: %g\n", x1)
TESTOUT("x2: %g\n", x2)
TESTOUT("x3: %g\n", x3)
TESTOUT("y1: %g\n", y1)
TESTOUT("y2: %g\n", y2)
TESTOUT("y3: %g\n", y3)

L:setRstate(rs)

m = L:colorMatrix()

nn = 100
s = 1/L:volume()

for g=0,15 do
  gn = qopqdp.groupName(g)

  t = 0
  for i=1,nn do
    m:random(rs)
    m:makeGroup(g)
    t = t + m:norm2()
  end
  printf("group %4s: %10.4f\n", gn, s*t/nn)

  t = 0
  for i=1,nn do
    m:randomGroup(g)
    t = t + m:norm2()
  end
  printf("group %4s: %10.4f\n", gn, s*t/nn)
end

g = L:gauge()
g:random()
t = g(1):norm2()
printf("g: %g\n", s*t)

f = L:force()
f:random()
t = f:norm2()
printf("f: %g\n", 0.25*s*t)
