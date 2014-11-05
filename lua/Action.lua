require 'Util'

local actionmt = {}
actionmt._type = "Action"
actionmt.__index = actionmt
local gaugecoeffs={}

-- kind =
---- "gauge"
---- "momentum"
function Action(opts)
  local self = {}
  tableCopyTo(self, opts)
  if not self.name then
    self.name = "A" .. self.kind
  end
  if self.kind == "gauge" then
    self.style = self.style or "custom"
    local gcfunc = gaugecoeffs[self.style]
    if not gcfunc then
      abort("unknown gauge action style %s\n", self.style)
    end
    self.coeffs = gcfunc(self)
    local plaq = self.coeffs.plaq or 0
    local rect = self.coeffs.rect or 0
    local pgm = self.coeffs.pgm or 0
    local adjplaq = self.coeffs.adjplaq or 0
    local vol = self.field.lattice.volume
    local xi0 = self.xi0 or 1
    local nd = #self.field.lattice.latsize
    local nd1 = nd-1
    self.xi0 = xi0
    self.act0 = vol*nd*(0.5*nd1*plaq + nd1*rect + 4*pgm + 0.5*nd1*adjplaq)
    self.act0 = self.act0*(1+xi0*xi0)/(2*xi0)
  elseif self.kind == "momentum" then
    local lat = self.momentum.lattice
    local nd = #lat.latsize
    local nc = self.momentum.nc
    local ngen = nc*nc-1 -- FIXME: assume SU for now
    self.act0 = -0.5*nd*lat.volume*ngen
  else
    error(string.format("unknown Action kind: %s\n", self.kind))
  end
  self.nupdate = 0
  clearStats(self, "getForce")
  return setmetatable(self,actionmt)
end

function actionmt.__tostring(self)
  local s = {}
  if self.tostringRecurse then
    local sf = self.field
    local sfr = sf.tostringRecurse
    sf.tostringRecurse = true
    s[#s+1] = tostring(sf)
    sf.tostringRecurse = sfr
  end
  s[#s+1] = self.name
  s[#s+1] = " = Action {\n"
  s[#s+1] = string.format("  kind = \"%s\",\n",self.kind)
  s[#s+1] = string.format("  style = \"%s\",\n",self.style)
  s[#s+1] = string.format("  field = %s,\n",self.field.name)
  s[#s+1] = "}\n"
  return table.concat(s)
end

function actionmt.Clone(self)
  local ga = tableCopy(self)
  return setmetatable(ga,actionmt)
end

function actionmt.Set(self, opts)
  if type(opts)=="table" then
    local t = opts._type
    if t=="Action" then
      tableCopyTo(self, opts)
    elseif t==nil then
      tableCopyTo(self, opts)
    else
      abort("unknown type for opts: %s\n", t)
    end
  else
    abort("unknown type for opts: %s\n", type(opts))
  end
end

function actionmt.Refresh(self)
  if self.kind == "momentum" then
    self.momentum:Random(self.variance)
  end
end

function actionmt.Action(self)
  local f = self.field
  local a = 0
  if self.kind == "gauge" then
    local ss,st = f.field:action(self.coeffs, self.xi0)
    a = self.beta*(self.act0-ss-st)
  end
  if self.kind == "momentum" then
    if self.variance then
      local _,fn = self.momentum:Norm2()
      for i=1,#fn do
	local si = 1/self.variance[i]
	a = a + 0.5*si*si*fn[i]
      end
    else
      a = 0.5*self.momentum:Norm2()
    end
    a = a + self.act0
  end
  return a
end

function actionmt.Heatbath(self, opts)
  if self.kind == "gauge" then
    local f = self.field
    local nrep = opts.nRepetitions or 1
    local nhb = opts.nHeatbath or 1
    local nor = opts.nOverrelax or 1
    local beta = self.beta
    local coeffs = self.coeffs
    f.field:heatbath(nrep, nhb, nor, beta, coeffs)
  end
end

function actionmt.GetForce(self, momentum, field, eps)
  local t0 = clock()
  assert(field==self.field)
  if not self.momentum then
    self.momentum = self.field:Momentum()
  end
  local m = momentum:GetField()
  local g = field:GetField()
  local f = self.momentum:GetField()
  f:zero()
  g:force(f, self.coeffs, self.beta, self.xi0)
  --momentum.field:update(f, eps)
  self.momentum:IncrementUpdateNumber()
  --m:combine({m,f}, {1,eps})
  for i=1,#m do 
    local mi = m(i)
    mi:combine({mi,f(i)},{1,eps})
  end
  momentum:IncrementUpdateNumber()
  local gf2 = eps*eps*f:norm2()
  local gfi = eps*f:infnorm()
  updateStats(self, "getForce",
	      {seconds=(clock()-t0),nUpdate=1,rms=gf2,max=gfi})
  if self.printForce then
    local ff2 = m:norm2()
    printf("%s:GetForce: norm %g\tinf %g\tff2: %g\n", self.name, gf2, gfi, ff2)
  end
end

function gaugecoeffs.custom(p)
  return { plaq=(p.plaq or 0), rect=(p.rect or 0),
	   pgm=(p.pgm or 0), adjplaq=(p.adjPlaq or 0) }
end
function gaugecoeffs.plaquette(p)
  return { plaq=1, rect=0, pgm=0, adjplaq=0 }
end
function gaugecoeffs.plaquetteAdj(p)
  return { plaq=1, rect=0, pgm=0, adjplaq=p.adjPlaq }
end
function gaugecoeffs.symanzikTree(p)
  local u0 = p.u0 or 1
  local u2 = u0*u0
  local c = { plaq=1, pgm=0, adjplaq=0 }
  c.rect = -1/(20*u2)
  return c
end
function gaugecoeffs.iwasaki(p)
  local u0 = p.u0 or 1
  local u2 = u0*u0
  local c1 = -0.331
  local c = { pgm=0, adjplaq=0 }
  c.plaq = 1 - 8*c1
  c.rect = c1/u2
  return c
end
function gaugecoeffs.dbw2(p)
  local u0 = p.u0 or 1
  local u2 = u0*u0
  local c1 = -1.4067
  local c = { pgm=0, adjplaq=0 }
  c.plaq = 1 - 8*c1
  c.rect = c1/u2
  return c
end
function gaugecoeffs.symanzik1Loop(p)
  local u0 = p.u0 or 1
  local nf = p.nf or 0
  local u2 = u0*u0
  local lu0 = math.log(u0)
  local c = { plaq=1, adjplaq=0 }
  c.rect = -(1 - (0.6264-0.4742*nf)*lu0 ) / (20*u2)
  c.pgm = (0.0433-0.012*nf)*lu0 / u2
  return c
end
function gaugecoeffs.symanzik1LoopHisq(p)
  local u0 = p.u0 or 1
  local nf = p.nf or 0
  local u2 = u0*u0
  local lu0 = math.log(u0)
  local c = { plaq=1, adjplaq=0 }
  c.rect = -(1 - (0.6264-1.1746*nf)*lu0 ) / (20*u2)
  c.pgm = (0.0433-0.0156*nf)*lu0 / u2
  return c
end
