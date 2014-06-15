local testfn = arg[1]
local tmpfn = testfn:gsub(".*/","") .. ".tmp"
printf("#runtest: using temp file %s\n", tmpfn)
function TESTOUT(...)
  printf("!")
  printf(...)
end
function TESTON()
  printf("!TESTON\n")
end
function TESTOFF()
  printf("!TESTOFF\n")
end
local pats = {}
function TESTPAT(s)
  pats[#pats+1] = s
end
local patrange = {}
function TESTRANGE(b,e)
  patrange[#patrange+1] = {b,e}
end

qhmc.remapout(tmpfn)

dofile(testfn)

qhmc.restoreout()

local istest = 0
for l in io.lines(tmpfn) do
  if l:sub(1,1) == '!' then
    if l == "!TESTON" then istest = istest + 1
    elseif l == "!TESTOFF" then istest = istest - 1
    else printf("%s\n", l:sub(2)) end
  else
    for i=1,#patrange do
      if l:match(patrange[i][2]) then istest = istest - 1 end
    end
    local dotest = (istest>0)
    for i=1,#pats do
      if l:match(pats[i]) then dotest = true end
    end
    if dotest then
      printf("%s\n", l)
    else
      printf("#%s\n", l)
    end
    for i=1,#patrange do
      if l:match(patrange[i][1]) then istest = istest + 1 end
    end
  end
end

os.remove(tmpfn)
