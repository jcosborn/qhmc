nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)
L:seed(98765321)

printf("BEGIN_TEST\n")

--f = L:real()
f = L:cmatrix()
f:random()
printf("%g\n", f:norm2())

r = f:lnorm2()
n = r:sum()
printf("%g\n", n)

a = f:norm2("timeslices")
myprint(a,"\n")

subs = L:subset("timeslices")
a = f:norm2(subs)
myprint(a,"\n")

subs = {L:subset("even"),L:subset("odd")}
a = f:norm2(subs)
myprint(a,"\n")

printf("END_TEST\n")
