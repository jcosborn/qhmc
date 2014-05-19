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
TESTOUT("file metadata: %s\n", md)

md = f:read(reader)
TESTOUT("record metadata: %s\n", md)

os.remove(fn)
