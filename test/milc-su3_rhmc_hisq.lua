local indir = arg[1]:gsub("[^/]*.lua","")
local srcfile = indir .. "milc-su3_rhmc_hisq.in"
local infile = "milc-su3_rhmc_hisq.in.tmp"
local infh = io.open(infile, "w+")
for l in io.lines(srcfile) do
  if l:match("load_rhmc_params") then
    l = "load_rhmc_params " .. indir .. "milc-su3_rhmc_hisq.rat"
  end
  infh:write(l .. "\n")
end
infh:close()

_G.infile = infile
require 'milc'

TESTZEROTOL(1e-16)
TESTPAT("^plaq")
TESTPATFMTSPACE("^S","%.4f")
TESTPATFMT("deltaS","%.5f")
TESTPATFMT("^1","%.5g")
TESTPAT("MEAS")
