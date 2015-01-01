require 'matrix'

local Ns = 4
local gammaelem = {}
-- gamma_ij -> gammaelem_i
gammaelem[0] = { {1,0},{2,0},{3,0},{4,0} }
gammaelem[1] = { {4,1},{3,1},{2,3},{1,3} }
gammaelem[2] = { {4,2},{3,0},{2,0},{1,2} }
gammaelem[4] = { {3,1},{4,3},{1,3},{2,1} }
gammaelem[8] = { {3,0},{4,0},{1,0},{2,0} }
local function gtimes(a, b)
  local c = {}
  for i=1,4 do
    c[i] = { b[a[i][1]][1], (a[i][2]+b[a[i][1]][2])%4 }
  end
  return c
end
for g=0,15 do
  if not gammaelem[g] then
    local m = gammaelem[0]
    local h,b = g,1
    repeat
      if h%2==1 then
	--printf("%i\t%i\t%i\n", g, h, b)
	m = gtimes(m, gammaelem[b])
      end
      h,b = math.floor(h/2),2*b
    until h == 0
    gammaelem[g] = m
  end
  --myprint("g[",tostring(g),"]=",gammaelem[g],"\n")
end
local gammamatrix = {}
local I = complex(0,1)
gammamatrix[0] = matrix(4,4,{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1})
gammamatrix[1] = matrix(4,4,{0,0,0,I, 0,0,I,0, 0,-I,0,0, -I,0,0,0})
gammamatrix[2] = matrix(4,4,{0,0,0,-1, 0,0,1,0, 0,1,0,0, -1,0,0,0})
gammamatrix[4] = matrix(4,4,{0,0,I,0, 0,0,0,-I, -I,0,0,0, 0,I,0,0})
gammamatrix[8] = matrix(4,4,{0,0,1,0, 0,0,0,1, 1,0,0,0, 0,1,0,0})
for g=0,15 do
  if not gammamatrix[g] then
    local m = gammamatrix[0]
    local h,b = g,1
    repeat
      if h%2==1 then
	--printf("%i\t%i\t%i\n", g, h, b)
	m = m * gammamatrix[b]
      end
      h,b = math.floor(h/2),2*b
    until h == 0
    gammamatrix[g] = m
  end
  --myprint("g[",tostring(g),"]=",gammaelem[g],"\n")
end
local function gamma(g)
  return gammamatrix[g]
end
local function gammaDagger(g)
  local m = gammamatrix[g]
  local ng = 0.5*(g%2) + 0.25*(g%4) + 0.125*(g%8) + 0.125*(g%16)
  if ng==2 or ng==3 then m = -1*m end
  return m
end
--print(gamma(15))

function pointProp(L, w, pt, smearsrc, smeardest, mass, resid, opts, oldprop)
  local src = w:quark()
  local dest = w:quark()
  local Nc = src:nc()
  local Ns = 4
  local dest2 = oldprop
  if not dest2 then
    dest2 = {}
    for spin2 = 1,Ns do
      dest2[spin2] = {}
      for spin = 1,Ns do
	local tcm = L:colorMatrix()
	dest2[spin2][spin] = tcm
	tcm:zero()
      end
    end
  end
  local cvCS = {}
  local cvSC = {}
  for spin = 1,Ns do
    cvSC[spin] = {}
    for color = 1,Nc do
      if spin==1 then cvCS[color] = {} end
      cvCS[color][spin] = L:colorVector()
      cvSC[spin][color] = cvCS[color][spin]
    end
  end
  for spin = 1,Ns do
    for spin2 = 1,Ns do
      dest2[spin2][spin]:splitColor(cvSC[spin2])
    end
    for color = 1,Nc do
      printf("source spin: %i  color: %i\n", spin, color)
      src:zero()
      -- point({coord},color,spin,value)
      src:point(pt,color-1,spin-1,1)

      if smearsrc then
	printf("src norm2: %g\n", src:norm2())
	smearsrc(src)
      end
      printf("src norm2: %g\n", src:norm2());
      io.stdout:flush();

      dest:combineSpin(cvCS[color])
      t0 = qopqdp.dtime()
      w:solve(dest, src, mass, resid, "all", opts) 
      dt = qopqdp.dtime() - t0
      mf = 1e-6 * w:flops() / dt
      printf("its: %g  secs: %g  Mflops: %g\n", w:its(), dt, mf)
      if smeardest then
	printf("dest norm2: %.16g\n", dest:norm2())
	smeardest(dest)
      end
      printf("dest norm2: %.16g\n", dest:norm2())
      dest:splitSpin(cvCS[color])
    end
    for spin2 = 1,Ns do
      dest2[spin2][spin]:combineColor(cvSC[spin2])
    end
  end
  return dest2
end

function randomPoint(L)
  local x = {}
  for i=1,#L do
    local li = L(i)
    x[i] = math.floor(li * qopqdp.random())
    assert(x[i]>=0 and x[i]<li)
  end
  return x
end

function randMomSource(L, chiral)
  local src = L:diracFermion()
  local Nc = src:nc()
  local Ns = 4
  local mom = randomPoint(L)
  local cmom = Nc * qopqdp.random()
  assert(cmom>=0 and cmom<Nc)
  local smom = Ns * qopqdp.random()
  assert(smom>=0 and smom<Ns)
  src:momentum(mom, cmom, smom)
  local ch = 0
  if chiral then
    ch = (qopqdp.random()<0.5) and 1 or -1
    local t = src:clone()
    t:gamma(15, src)
    src:combine({src,t},{0.5,0.5*ch})
  end
  return src, mom, cmom, smom, ch
end

function mydot(v, g, color, spin, qt)
  local spin2,phase = table.unpack(gammaelem[g][spin])
  local k = (color-1)*Ns + spin
  local k2 = (color-1)*Ns + spin2
  local vk = v[k]
  if g ~= 0 then
    qt:gamma(g, vk)
    vk = qt
  end
  local t
  if phase==0 or phase==2 then
    if k==k2 and g==0 then
      t = v[k2]:norm2("timeslices")
    else
      t = v[k2]:reDot(vk,"timeslices")
    end
  else
    t = v[k2]:dot(vk,"timeslices")
    for i in pairs(t) do t[i] = t[i].i end
  end
  --printf("%i\t%i\t%g\n", k2, k, t[1])
  if phase==1 or phase==2 then
    for i=1,#t do t[i] = -t[i]; end
  end
  return t
end

function wilsonMesons(dest)
  local mesons = {}
  local qt = dest[1]:clone()
  local Nc = qt:nc()
  for g = 0,#gammaelem do
    mesons[g] = {}
    for color = 1,Nc do
      for spin = 1,Ns do
	local t = mydot(dest, 15-g, color, spin, qt)
	--printf("%g\n", t[1])
	for j=1,#t do
	  mesons[g][j] = (mesons[g][j] or 0) + t[j]
	end
      end
    end
  end
  return mesons
end

local scalemt = {}
scalemt.__index = scalemt
function scale(z,x)
  return setmetatable({z=z,x},scalemt)
end
local function iszero(x)
  return (x==0) or
    (type(x)~="number" and getmetatable(x)~=scalemt and x:norm2()==0)
end
function scalemt.__tostring(x)
  return string.format("%s x %s", x.z, x[1])
end
function scalemt.__add(x,y)
  if iszero(x) then return y end
  if iszero(y) then return x end
  local r = x[1]:clone()
  r:combine({x[1],y[1]},{x.z,y.z})
  return setmetatable({z=1,r},scalemt)
end
function scalemt.__mul(x,y)
  local z,v
  if getmetatable(x)~=scalemt then
    z = x * y.z
    if iszero(z) then return 0 end
    v = y[1]
  elseif getmetatable(y)~=scalemt then
    z = x.z * y
    if iszero(z) then return 0 end
    v = x[1]
  else
    z = x.z * y.z
    if iszero(z) then return 0 end
    v = x[1]:clone()
    v:mul(x[1],y[1])
  end
  return scale(z,v)
end
function scalemt.trace(x)
  local r = x[1]:trace()
  return scale(x.z, r)
end
function scalemt.sum(x,...)
  local r = x[1]:sum(...)
  return x.z * matrix(#r,1,r)
end
function scalemt.dot(x,y,...)
  local z
  if type(x.z)=="number" then z = x.z * y.z
  else z = x.z:conj() * y.z end
  local r = x[1]:dot(y[1],...)
  --printf("dot: %s\n", z)
  --printf("dotr: %s %s\n", z, r[1])
  return z * matrix(#r,1,r)
end
function scalemt.contract(x,y,...)
  local z = x.z * y.z
  local r = x[1]:contract(y[1],...)
  --printf("dot: %s\n", z)
  --printf("dotr: %s %s\n", z, r[1])
  return z * matrix(#r,1,r)
end

local function phase2z(p)
  local ps = {complex(0,1), -1, complex(0,-1), 1}
  return ps[((p+3)%4)+1]
end

local function prop(v)
  local p = matrix(4,4)
  for i=1,Ns do
    for j=1,Ns do
      p(i,j,scale(1,v[i][j]))
    end
  end
  return p
end
function ppropClone(v)
  local p = {}
  for i=1,Ns do
    local pi = {}
    p[i] = pi
    for j=1,Ns do
      local e = v[i][j]
      pi[j] = e:clone()
    end
  end
  return p
end
function ppropSet(p,v)
  for i=1,Ns do
    for j=1,Ns do
      local ep = p[i][j]
      local ev = v[i][j]
      ep:set(ev)
    end
  end
end
function ppropMomentum(p,...)
  for i=1,Ns do
    for j=1,Ns do
      local ep = p[i][j]
      ep:momentum(...)
    end
  end
end

--local function 

function mydot2(v, g, spin, spin2)
  local jspin,phase = table.unpack(gammaelem[g][spin])
  local jspin2,phase2 = table.unpack(gammaelem[g][spin2])
  local t = v[spin2][jspin]:dot(v[jspin2][spin],"timeslices")
  --printf("%i\t%i\t%g\n", k2, k, t[1])
  local z = phase2z(phase+phase2)
  --printf("%i\t%i\t%i\t%s\t%s\n", g, spin, spin2, t[1], z)
  return t,z
end

function wilsonMesons2(dest)
  local mesons = {}
  for g = 0,#gammaelem do
    mesons[g] = {}
    for spin = 1,Ns do
      for spin2 = 1,Ns do
	local t,z = mydot2(dest, 15-g, spin, spin2)
	--printf("%g\n", t[1])
	for j=1,#t do
	  mesons[g][j] = (mesons[g][j] or 0) + z*t[j]
	end
      end
    end
    for j=1,#mesons[g] do
      mesons[g][j] = mesons[g][j].r
    end
  end
  return mesons
end

function wilsonMesons3(destu, destd)
  local mesons = {}
  local pu = prop(destu)
  local pd = pu
  if destd then pd = prop(destd) end
  for g = 0,15 do
    mesons[g] = {}
    local p2 = gammaDagger(15-g) * pu * gamma(15-g)
    local t = pd:dot(p2,"timeslices")
    for i=1,t.nr do mesons[g][i] = t(i,1).r end
  end
  return mesons
end

function contract24(p1,p2)
  local p0 = matrix(4,4)
  local r = p1(1,1)[1]:clone()
  for i1=1,4 do
    for i2=1,4 do
      local s = r:clone()
      s:zero()
      for j=1,4 do
	local e1 = p1(i1,j)
	local e2 = p2(i2,j)
	local z = e1.z * e2.z
	r:transcross({e1[1],e2[1]})
	--r:cross({e1[1],e2[1]})
	s:combine({s,r},{1,z})
      end
      p0(i1,i2, scale(1,s))
    end
  end
  return p0
end

-- add quark cross term
-- take spin trace, then add back along spin diagonal
function add_cross(p)
  local s = p(1,1)[1]:clone()
  local t,ones = {},{}
  for i=1,4 do
    t[i] = s:clone()
    local e = p(i,i)
    t[i]:combine({e[1]},{e.z})
    ones[i] = 1
  end
  s:combine(t,ones)
  for i=1,4 do
    local e = p(i,i)
    e[1]:combine({t[i],s},{1,1})
    e.z = 1
  end
end

function wilsonBaryons3(destu, destd)
  local Nc = destu[1][1]:nc()
  if Nc ~= 3 then return nil end
  local baryons = {}
  local pu = prop(destu)
  local pd = pu
  if destd then pd = prop(destd) end
  -- gamma(5) = gamma_x gamma_z = C gamma_5
  local p2 = pd * gamma(5)
  local p3 = contract24(pu, p2)
  p2 = nil
  collectgarbage()
  local p4 = p3 * gamma(5)
  add_cross(p4)
  p3 = nil
  collectgarbage()
  local p5 = p4 * pu
  p4 = nil
  collectgarbage()
  p5 = p5:mapMethod("trace")
  p5 = p5:mapMethod("sum", "timeslices")
  for g = 0,15 do
    local p6 = p5 * gamma(g)
    --local c = p6:rtrace()
    --local t = c:sum("timeslices")
    local t = p6:trace()
    baryons[g] = {}
    for i=1,t.nr do baryons[g][i] = t(i,1) end
  end
  p5 = nil
  collectgarbage()
  return baryons
end
