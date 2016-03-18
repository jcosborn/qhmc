require 'Util'

function weakfield(g, w)
  local f = g:lattice():force()
  f:randomTAH()     -- Generate a random traceless antihermitian matrix.
  g:unit()          -- Reset to the unit gauge.
  g:update(f, w);   -- Multiply the unit gauge by exp(wt*force).
end

--[[ old stuff, not working yet
-- changed sign convention on paths
paths0 = {
  {-1,-2, 1, 2,-3,-4, 3, 4}
}
coeffs0 = { {1} } -- coeff[path][rep]

local paths1 = {
  {-1,-2,-3, 2, 1,-4,-1, 4, 1, 3},
  {-1,-2,-3, 1,-4, 3,-1, 4, 1, 2}
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
  local s = 0
  local np = 0
  for i=1,#paths do
    local pt = {}
    load(pt, paths[i], nd)
    np = np + #pt
    for k=1,#pt do
      -- FIXME: need other reps
      local l = g:loop(pt[k])
      s = s + coeffs[i][1]*l
    end
  end
  s = s/np
  return s
end
--]]

function plaq(g)
  local ss,st = g:action({plaq=1})
  local s = g:lattice():volume()*g:nc()
  return ss/s, st/s
end

function wflow(u, coeffs, eps, nsteps)
  -- u1 = exp((eps/4)f0) u0
  -- u2 = exp((eps*8/9)f1-(eps*17/36)f0) u1
  -- u3 = exp((eps*3/4)f2-(eps*8/9)f1+(eps*17/36)f0) u2
  local f = qopqdp.force()
  local ft = qopqdp.force()
  eps = u:nc() * eps -- compensate for 1/Nc convention in force
  for i=1,nsteps do
    --[[
    u:force(f, coeffs)
    u:update(f, eps)
    --]]
    --[[
    local alpha = 0.5
    u:force(ft, coeffs)
    u:update(ft, eps*alpha)
    u:force(f, coeffs)
    f:update(ft, 2*alpha-1-2*alpha*alpha)
    u:update(f, eps*0.5/alpha)
    --]]
    --f:zero()
    u:force(f, coeffs)
    --for i=1,4 do
    --printf("f[%i]: %s\n", i, (17*eps/36)^2*f(i):norm2())
    --end
    u:update(f, eps/4)
    ft:zero()
    u:force(ft, coeffs)
    f:fupdate(ft, -32/17)
    u:update(f, eps*-17/36)
    u:force(ft, coeffs)
    f:fupdate(ft, 27/17)
    u:update(f, eps*17/36)
  end
end

function plaqE(g)
  local ss,st = plaq(g)
  local nd = #(g:lattice())
  local e = 0.5*nd*(nd-1)*(2-ss-st)
  return e
end

-- get F_{mu,nu}
-- f: (out) F_{mu,nu}
-- g: (in) gauge field
-- mu: (in) mu
-- nu: (in) nu
-- order: (in) 0 (a^2) or 1 (a^4)
function fmunu(f, g, mu, nu, order)
  order = order or 0
  local ps = { { -mu, -nu,  mu,  nu },
	       {  nu, -mu, -nu,  mu },
	       {  mu,  nu, -mu, -nu },
	       { -nu,  mu,  nu, -mu } }
  --local ps = { { -mu,  -nu,  mu, nu } }
  local f0 = f:clone()
  f:zero()
  for i=1,#ps do
    f0:unit()
    f0:transport(f0, g, ps[i])
    f:combine({f,f0},{1,0.25})
    if order>0 then
      -- 1 loop: 135/90 = 1.5
      -- 2 loop: -27/180 = -3/20 = -0.15 = -0.1*1.5
      -- 3 loop: 1/90 = 1.5*1/135
      local p2 = {}; for j=1,#ps[i] do p2[#p2+1] = ps[i][j]; p2[#p2+1] = ps[i][j] end
      f0:unit()
      f0:transport(f0, g, p2)
      f:combine({f,f0},{1,-0.025})
      local p3 = {}; for j=1,#ps[i] do p3[#p3+1] = ps[i][j]; p3[#p3+1] = ps[i][j]; p3[#p3+1] = ps[i][j] end
      f0:unit()
      f0:transport(f0, g, p3)
      f:combine({f,f0},{1,1/540})
    end
  end
  --local t = f:trace():sum()/f:lattice():volume()
  --printf("trace: %s\n", t)
  --local AH = 3 -- should be globally set
  local TAH = 11 -- should be globally set
  f:makeGroup(TAH)
  if order>0 then
    f0:combine({f},{complex(0,-1.5)})
  else
    f0:combine({f},{complex(0,-1)})
  end
  f:set(f0)
end

-- comment from original C code
-- ESW addition 8-21-2014
-- Get symmE and topology at O(a^2) and optionally also O(a^4)
-- The O(a^2) is equivalent to Luscher's measurement of symmE
-- defined in his Wilson flow paper. symmQ is defined with
-- the same definition of the field-strength tensor.
-- The O(a^4) improved is based on arXiv:0203008. I use the 3-loop
-- definition, using the 1x1, 2x2, and 3x3 clover loop to build the
-- improved op. 
-- 1: gauge field
-- 2: 0 for just O(a^2) calculation, 1 for O(a^4) calculation.
-- 3: optional subset support.
function symmEQ(g, order, subsets)
  local f = {}
  local n = 0
  -- 1:1,2 2:1,3 3:2,3 4:1,4 5:2,4 6:3,4
  for nu=2,4 do
    for mu=1,nu-1 do
      n = n + 1
      f[n] = L:colorMatrix()
      fmunu(f[n], g, mu, nu, order)
    end
  end
  local function contract(r, a, f1, f2, subs)
    local z
    if subs==nil then
      z = f1:contract(f2)
    else
      z = f1:contract(f2, subs)
    end
    --printf("%i: %s\n", i, z)
    --z = f[i]:dot(f[i])
    if type(z)~="table" then z = {z} end
    for i=1,#z do
      r[i] = (r[i] or 0) + a * z[i].r
    end
  end
  -- First, get symmetric E!
  -- E = 1/4 F_{mu nu}^a F_{mu nu}^a
  local ses,set = {},{}
  local sef = 1/L:volume()
  for i=1,n do
    if i<=3 then
      contract(ses, sef, f[i], f[i], subsets)
    else
      contract(set, sef, f[i], f[i], subsets)
    end
  end
  -- Next, get symmQ!
  -- Q = 1/(32 pi^2) eps_uvrs F_uv^a F_rs^a
  local sq = {}
  local sqf = 1/(4*math.pi^2)
  contract(sq,  sqf, f[1], f[6], subsets)
  contract(sq, -sqf, f[2], f[5], subsets)
  contract(sq,  sqf, f[3], f[4], subsets)
  if #ses == 1 then ses = ses[1] end
  if #set == 1 then set = set[1] end
  if #sq == 1 then sq = sq[1] end
  return ses, set, sq
end
