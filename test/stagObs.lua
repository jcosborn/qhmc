--doTrace = true
seed = 987654321
nx = 4
nt = 8
gfix_prec = 1e-8
cg_prec = 1e-8

TESTON()
require 'staggeredObservablesConn'
qopqdp.verbosity(1);
require 'staggeredObservablesDisc'
TESTOFF()

TESTZEROTOL(1e-15)
TESTIGNORE("time:")
TESTIGNORE("Time:")
TESTIGNORE("secs:")
TESTPATFMT(".*","%.6g")
TESTPATFMT("^0","%.5g")
