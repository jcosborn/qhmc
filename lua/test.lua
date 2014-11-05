qopqdp.lattice{4,4,4,4}
qopqdp.profile(0)
qopqdp.verbosity(0)
for k,v in ipairs(qopqdp.lattice()) do print(k,v) end

qopqdp.seed(54321)

g = qopqdp.gauge()
g:unit()
g2 = g:copy()
g:set(g2)

coef = { plaq=1, rect=0, pgm=0 }
acts,actt = g:action(coef)
print(acts, actt)

f = qopqdp.force()
g:force(f, coef)
print(f:norm2())
print(f:infnorm())

g:random()

acts,actt = g:action(coef)
print(acts, actt)

g:force(f, coef)
print(f:norm2())
print(f:infnorm())

h = qopqdp.hisq()
q = h:quark()
q2 = h:quark()
q3 = h:quark()
u0 = 0.9
h:set(g, u0)

q:random()
print(q:norm2())
mass = 0.1
h:D(q2, q, mass)
print(q2:norm2())
h:D(q2, q, mass, "even", "all")
print(q2:norm2())
h:Ddag(q2, q, mass, "even", "all")
print(q2:norm2())

h:solve({q}, q2, {mass}, 1e-5, "even")
print("solve time: ", h:time(), " mflops: ", 1e-6*h:flops()/h:time())
h:D(q3, q, 0.5, "all", "even")
h:force(f, {q3}, {0.5/mass})
print("force time: ", f:time(), " mflops: ", 1e-6*f:flops()/f:time())

print(q2:Re_dot(q)*0.25/mass)
