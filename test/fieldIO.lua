nx = 4
nt = 4
ls = { nx, nx, nx, nt }
L = qopqdp.lattice(ls)

fn = "test.lat"
writer = L:writer(fn, "test")

f = L:rstate()
--f = L:real()
--f = L:cmatrix()
if f.zero then f:zero() end

--f:write(writer, "test2", "D")
f:write(writer, "test2")

reader,md = L:reader(fn)
print(md)

md = f:read(reader)
print(md)

os.remove(fn)
