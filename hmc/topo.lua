require 'Util'

paths0 = {
  { 1, 2,-1,-2, 3, 4,-3,-4}
}
coeffs0 = { {1} } -- coeff[path][rep]

local paths1 = {
  { 1, 2, 3,-2,-1, 4, 1,-4,-1,-3},
  { 1, 2, 3,-1, 4,-3, 1,-4,-1,-2}
}
--local coeffs1 = { {0.09029782,0.4884618},{-0.1786298,0.4126987} } -- SU(2)
local coeffs1 = { {0.07872507,0.3173630},{-0.1888383,0.2854577} } -- SU(3)

function permGen(a, n)
  if n == 0 then
    coroutine.yield(a)
  else
    local n1 = #a
    local n0 = 1 + n1 - n
    for i=n0,n1 do
      a[n0], a[i] = a[i], a[n0]
      permGen(a, n-1)
      a[n0], a[i] = a[i], a[n0]
    end
  end
end

function permIter(n)
  local p = {}
  for i=1,n do p[i] = i end
  return coroutine.wrap(function () permGen(p, n) end)
end

function binaryIter(n)
  local p = {}
  for i=1,n do p[i] = 0 end
  p[1] = -1
  return function()
	   for i=1,n do
	     p[i] = p[i] + 1
	     if p[i]==2 then p[i] = 0 else return p end
	   end
	   return nil
	 end
end

function load(lt, l, nd)
  for p in permIter(nd) do
    for r in binaryIter(nd) do
      local t = {}
      for i=1,#l do
	local d = l[i]
	if d>0 then
	  t[i] = p[d] * (1-2*r[d])
	else
	  t[i] = -p[-d] * (1-2*r[-d])
	end
      end
      lt[#lt+1] = t
    end
  end
end

function pathDirOrder(x) return 0.25+math.abs(x-0.25) end
function pathDirBefore(x,y) return pathDirOrder(x)<pathDirOrder(y) end
function pathBefore(x,y)
  for i=1,#x do
    if pathDirOrder(x[i])<pathDirOrder(y[i]) then return true end
    if pathDirOrder(y[i])<pathDirOrder(x[i]) then return false end
  end
  return false
end
function pathCanonicalStart(p0,i)
  local n = #p0
  local p = {}
  local s = (p0[i]<0) and -1 or 1
  for l=1,n do
    p[l] = s*p0[(i+s*(l-1)+n-1)%n + 1]
  end
  return p
end
function pathCanonical(p0)
  local n = #p0
  local p = pathCanonicalStart(p0,1)
  for i=2,n do
    local t = pathCanonicalStart(p0,i)
    if pathBefore(t,p) then p = t end
  end
  --myprint(p0,"\n")
  --myprint(p,"\n")
  return p
end
function pathAdd(t, p)
  local n = #t
  local i0 = 0
  local i1 = n+1
  while true do
    if i1==i0+1 then
      for j=n,i1,-1 do t[j+1] = t[j] end
      t[i1] = p
      break
    end
    local i = math.floor(0.5*(i0+i1+1))
    if pathBefore(p,t[i]) then i1 = i else
      if pathBefore(t[i],p) then i0 = i else
	break
      end
    end
  end
end

function pathTrace(p0)
  local p1 = {}
  for i=1,#p0 do
    local t = pathCanonical(p0[i])
    pathAdd(p1, t)
  end
  return p1
end

--[[
local paths = {}
local coeffs = {}
--load(paths, {1,2,-1,-2}, 4)
load(paths, paths0[1], 4)
myprint(paths, "\n")
paths = pathTrace(paths)
myprint(paths, "\n")
--]]

function pathDo(g, nd, paths, coeffs)
  local pr,pi = 0,0
  local np = 0
  for i=1,#paths do
    local p = {}
    load(p, paths[i], nd)
    np = np + #p
    for k=1,#p do
      -- FIXME: need other reps
      local lr,li = g:loop(p[k])
      pr = pr + coeffs[i][1]*lr
      pi = pi + coeffs[i][1]*li
    end
  end
  pr = pr/np
  pi = pi/np
  return pr,pi
end
