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

TESTPAT("^plaq")
TESTPAT("^S")
TESTPAT("deltaS")
TESTPAT("unitarity")
TESTPAT("^1,")
TESTPAT("MEAS")

TESTRANGE("^Local Pions","^MEASplaq_ss")
