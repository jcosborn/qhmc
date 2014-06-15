local matrix_mt = {}
matrix_mt.__index = matrix_mt
function matrix(nr,nc,v)
  local m = {nr=nr,nc=nc,n=nc*nr,e={}}
  if v then for i=1,m.n do m.e[i] = v[i] end
  else for i=1,m.n do m.e[i] = 0 end end
  return setmetatable(m,matrix_mt)
end
function matrix_mt.set(m,i,j,v)
  m.e[(i-1)*m.nc+j] = v
end
function matrix_mt.__call(m,i,j,v)
  if v==nil then
    return m.e[(i-1)*m.nc+j]
  end
  m.e[(i-1)*m.nc+j] = v
end
function matrix_mt.__tostring(m)
  local s = ""
  s = s .. string.format("%i x %i matrix:\n", m.nr, m.nc)
  for i=1,m.nr do
    for j=1,m.nc do
      s = s .. string.format(" %i\t%i\t%s\n", i, j, tostring(m(i,j)))
    end
  end
  return s
end
function matrix_mt.__add(m1,m2)
  if getmetatable(m1)~=matrix_mt then
    local n1,n2 = m2.nr,m2.nc
    local r = matrix(n1,n2)
    for i=1,n1 do
      for j=1,n2 do
	local k = (i-1)*n2 + j
	r.e[k] = m1 + m2.e[k]
      end
    end
    return r
  end
  if getmetatable(m2)~=matrix_mt then
    local n1,n2 = m1.nr,m1.nc
    local r = matrix(n1,n2)
    for i=1,n1 do
      for j=1,n2 do
	local k = (i-1)*n2 + j
	r.e[k] = m1.e[k] + m2
      end
    end
    return r
  end
  if m1.nr ~= m2.nr or m1.nc ~= m2.nc then
    print("error incompatible dimensions")
    return nil
  end
  local r = matrix(m1.nr,m1.nc)
  for i=1,r.n do r.e[i] = m1.e[i] + m2.e[i] end
  return r
end
function matrix_mt.__sub(m1,m2)
  if m1.nr ~= m2.nr or m1.nc ~= m2.nc then
    print("error incompatible dimensions")
    return nil
  end
  local r = matrix(m1.nr,m1.nc)
  for i=1,r.n do r.e[i] = m1.e[i] - m2.e[i] end
  return r
end
function matrix_mt.__mul(m1,m2)
  if getmetatable(m1)~=matrix_mt then
    local n1,n2 = m2.nr,m2.nc
    local r = matrix(n1,n2)
    for i=1,n1 do
      for j=1,n2 do
	local k = (i-1)*n2 + j
	r.e[k] = m1 * m2.e[k]
      end
    end
    return r
  end
  if getmetatable(m2)~=matrix_mt then
    local n1,n2 = m1.nr,m1.nc
    local r = matrix(n1,n2)
    for i=1,n1 do
      for j=1,n2 do
	local k = (i-1)*n2 + j
	r.e[k] = m1.e[k] * m2
      end
    end
    return r
  end
  local n1,n2,n3 = m1.nr,m1.nc,m2.nc
  if n2 ~= m2.nr then
    print("error incompatible dimensions")
    return nil
  end
  local r = matrix(n1,n3)
  for i=1,n1 do
    for j=1,n3 do
      local t = m1.e[(i-1)*n2+1] * m2.e[j]
      for k=2,n2 do
	t = t + m1.e[(i-1)*n2+k] * m2.e[(k-1)*n3+j]
      end
      r.e[(i-1)*n3+j] = t
    end
  end
  return r
end
function matrix_mt.trace(m)
  if m.nr ~= m.nc then
    print("error incompatible dimensions")
    return nil
  end
  local r = 0
  for i=1,m.nr do
    r = r + m(i,i)
  end
  return r
end
function matrix_mt.rtrace(m)
  if m.nr ~= m.nc then
    print("error incompatible dimensions")
    return nil
  end
  local r = 0
  for i=1,m.nr do
    r = r + m(i,i):trace()
  end
  return r
end
function matrix_mt.transpose(m)
  local r = matrix(m.nc,m.nr)
  for i=1,m.nr do
    for j=1,m.nc do
      r(j,i, m(i,j))
    end
  end
  return r
end
function matrix_mt.dot(m1,m2,...)
  if m1.nr ~= m2.nr or m1.nc ~= m2.nc then
    print("error incompatible dimensions")
    return nil
  end
  local s = 0
  for i=1,m1.nr do
    for j=1,m1.nc do
      s = s + m1(i,j):dot(m2(i,j),...)
    end
  end
  return s
end
function matrix_mt.clone(m)
  local r = matrix(m.nr,m.nc)
  for i=1,m.n do r.e[i] = m.e[i] end
  return r
end
function mysqrt(x)
  if type(x)=="number" then return math.sqrt(x) end
  return x:sqrt()
end
function mynorm2(x)
  if type(x)=="number" then return x*x end
  return x:norm2()
end
local function hh_project(m, row, col, b)
  local cn = 0
  for i=row,m.nr do cn = cn + m(i,col)*m(i,col) end
  local alpha = mysqrt(cn)
  if mynorm2(m(row,col))>=0 then alpha = -alpha end
  local ri = 1/mysqrt(cn-m(row,col)*alpha)
  local v = {}
  v[0] = (m(row,col)-alpha)*ri
  for i=row+1,m.nr do v[i-row] = m(i,col)*ri end
  local function hhproj(a)
    for j=1,a.nc do
      local d = 0
      for i=row,a.nr do
	d = d + v[i-row]*a(i,j)
      end
      for i=row,m.nr do
	a(i,j, a(i,j)-d*v[i-row])
      end
    end
  end
  hhproj(m)
  if b then for k,v in pairs(b) do hhproj(v) end end
end
local function backsub(m, b)
  local x = matrix(m.nc,b.nc)
  for k=1,b.nc do
    for i=m.nc,1,-1 do
      local s = b(i,k)
      for j=i+1,m.nc do
	s = s - m(i,j)*x(j,k)
      end
      x(i,k, s/m(i,i))
    end
  end
  return x
end
function matrix_mt.solve(m, b)
  local a = m:clone()
  local c,bIsMatrix
  if getmetatable(b)==matrix_mt then
    c = b:clone()
    bIsMatrix = true
  else
    c = matrix(m.nr,1)
    for i=1,m.nr do c.e[i] = b[i] end
    bIsMatrix = false
  end
  local c0 = c:clone()
  for i=1,m.nc do
    hh_project(a,i,i,{c})
  end
  --print(m)
  local x = backsub(a,c)
  local r = c0 - m*x
  --print(c0)
  --print(m*x)
  --print(r)
  --print(m:transpose()*r)
  if not bIsMatrix then
    local r = {}
    for i=1,x.n do r[i] = x.e[i] end
    return r
  end
  return x
end
