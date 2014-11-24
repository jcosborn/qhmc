--doTrace = true
seed = 3190
nx = 4
nt = 8
beta = 6
beta_a = -0.25
nf = 8

hmcmasses = {0.1,0.2}
nfsteps = {20,20}
ngsteps = 100
lambdaG = 0.2
lambdaF = 0.2
gintalg = {type='2MNV',lambda=lambdaG}
fintalg = {type='2MNV',lambda=lambdaF}
prec = 2
faresid = 1e-8
mdresid = { 1e-6, 1e-7 }
use_prev_soln = 1
mixedRsq = 0
restart = 2000
tau = 1
warmup = 1
ntraj = 1
nlats = 2
checkReverse = true
doGfix = true
doSpectrum = true
doS4obs = true

require 'bsm'

TESTZEROTOL(1e-16)
TESTPAT("^plaq")
TESTPATFMT("^S","%.4f")
TESTPATFMT("deltaS","%.5f")
TESTPAT("^1,")
TESTPAT("MEAS")
TESTPATFMT("MEASpbp ","%.5g")
TESTPATFMT("^%d%s","%.4g")
TESTPATFMT("MEASplaq_","%.8g")
TESTPATFMT("MEASpbp_","%.4g")

TESTRANGE("^Local Pions","^MEASplaq_ss")
