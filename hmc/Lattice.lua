require 'Util'
require 'Field'

local latticemt = {}
latticemt._type = "Lattice"

function Lattice(latsize)
  local self = {}
  self.latsize = tableCopy(latsize)
  self.name = "L"
  self.defaultPrecision = 'D'
  self.defaultGroup = "SU"
  self.defaultNc = 3
  self.volume = 1
  for i=1,#latsize do self.volume=self.volume*latsize[i] end
  self.qdplat = qopqdp.lattice(latsize)
  return setmetatable(self,latticemt)
end

function latticemt.__tostring(self)
  local s = {}
  s[#s+1] = self.name
  s[#s+1] = " = Lattice{ "
  for i=1,#self.latsize do
    if i>1 then s[#s+1] = ", " end
    s[#s+1] = self.latsize[i]
  end
  s[#s+1] = " }\n"
  return table.concat(s)
end

function latticemt.__index(self, k)
  if type(k)=="number" then
    return self.latsize[k]
  end
  return latticemt[k]
end

function latticemt.__len(self)
  return #self.latsize
end

function latticemt.Set(self, opts)
  tableCopyTo(self, opts)
end

function latticemt.Seed(self, seed)
  qopqdp.seed(seed)
end

function latticemt.GaugeField(self, opts)
  local opts = opts or {}
  local opt = {}
  opt.precision = opts.precision or self.defaultPrecision
  opt.group = opts.group or self.defaultGroup
  opt.nc = opts.nc or self.defaultNc
  opt.kind = "gauge"
  return Field(self, opt)
end
