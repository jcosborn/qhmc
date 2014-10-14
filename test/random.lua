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
