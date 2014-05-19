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

qhmc.remapout(tmpfn)

dofile(testfn)

qhmc.restoreout()

local istest = false
for l in io.lines(tmpfn) do
  if l:sub(1,1) == '!' then
    if l == "!TESTON" then istest = true
    elseif l == "!TESTOFF" then istest = false
    else printf("%s\n", l:sub(2)) end
  else
    local dotest = istest
    for i=1,#pats do
      if l:match(pats[i]) then dotest = true end
    end
    if dotest then
      printf("%s\n", l)
    else
      printf("#%s\n", l)
    end
  end
end

os.remove(tmpfn)
