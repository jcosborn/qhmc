--doTrace = true
seed = 3190
nx = 4
nt = 8

TESTON()
require 'staggeredObservablesConn'
require 'staggeredObservablesDisc'
TESTOFF()

TESTZEROTOL(1e-16)
TESTIGNORE("time:")
TESTIGNORE("Time:")
TESTIGNORE("secs:")
TESTPATFMT(".*","%g")
