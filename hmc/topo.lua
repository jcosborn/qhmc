require 'Util'

local paths0 = {
  { 1, 2,-1,-2, 3, 4,-3,-4}
}
local coeffs0 = { {1} }
local paths1 = {
  { 1, 2, 3,-2,-1, 4, 1,-4,-1,-3},
  { 1, 2, 3,-1, 4,-3, 1,-4,-1,-2}
}
--local coeffs1 = { {0.09029782,0.4884618},{-0.1786298,0.4126987} } -- SU(2)
local coeffs1 = { {0.07872507,0.3173630},{-0.1888383,0.2854577} } -- SU(3)

local paths = {}
local coeffs = {}

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

load(paths, {1,2}, 2)
myprint(paths, "\n")
