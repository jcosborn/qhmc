local testfn = arg[1]
local rawfn = arg[2]
local tmpfn = rawfn
if not rawfn then
 tmpfn = testfn:gsub(".*/",""):gsub("%.lua$","") .. ".tmp"
end
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
local numpat = "[%+%-]?%d+%.?%d*[eE]?[%+%-]?%d*"
local numpatspace = "[%+%- ]?%d+%.?%d*[eE]?[%+%-]?%d*"
local zerotol = 0
local function zerotolfunc(x)
  local n = tonumber(x)
  if n and math.abs(n) < zerotol then x = " 0" end
  return x
end
function TESTZEROTOL(t)
  zerotol = t
end
local pats = {}
local gsubs = {}
function TESTPAT(s, ...)
  local n = #pats + 1
  pats[n] = s
  if select('#', ...) > 0 then
    gsubs[n] = {...}
  end
end
function TESTIGNORE(s)
  TESTPAT(s, ".*", "#%0")
end
function TESTPATFMT(s, f)
  local n = 1
  local function fmt()
    local r = f[n] or f
    if f[n+1] then n=n+1 end
    return r
  end
  TESTPAT(s, numpat,
	  function(x) return string.format(fmt(),x) end)
end
function TESTPATFMTSPACE(s, f)
  local n = 1
  local function fmt()
    local r = f[n] or f
    if f[n+1] then n=n+1 end
    return r
  end
  TESTPAT(s, numpatspace,
	  function(x) return string.format(fmt(),x) end)
end
local patrange = {}
function TESTRANGE(b,e)
  patrange[#patrange+1] = {b,e}
end

if rawfn then
  local t = ""
  for l in io.lines(testfn) do
    if l:match("^TEST[A-NP-Z]") then t = t..l end
  end
  load(t)()
else
  --os.remove(tmpfn)
  qhmc.remapout(tmpfn)
  for l in io.lines(tmpfn:gsub(".tmp$",".out")) do printf("%s",l) end
  dofile(testfn)
  qhmc.restoreout()
end

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
      if l:match(pats[i]) then
	dotest = true
	if gsubs[i] then
	  l = l:gsub(unpack(gsubs[i]))
	end
      end
    end
    if dotest then
      l = l:gsub(numpat, zerotolfunc)
      printf("%s\n", l)
    else
      printf("#%s\n", l)
    end
    for i=1,#patrange do
      if l:match(patrange[i][1]) then istest = istest + 1 end
    end
  end
end

--os.remove(tmpfn)
