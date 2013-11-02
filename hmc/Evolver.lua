require 'Util'

local evolvermt = {}
evolvermt._type = "Evolver"
evolvermt.__index = evolvermt
local intpat = {}
local getIntPat

-- kind = md
-- style =
--- recursive
--- simultaneous
--- leapfrog, OMF, custom
function Evolver(opts)
  local self = {}
  tableCopyTo(self, opts)
  if not self.name then
    self.name = "Evol" .. self.kind
  end
  if self.kind == "md" then
    self.fields = {self.field}
    self.momenta = {self.momentum}
    getIntPat(self)
  elseif self.kind == "mc" then
    --self.fields = {self.field}
    --self.momenta = {self.momentum}
    self.globalRand = qopqdp.random
  end
  clearStats(self, "run")
  return setmetatable(self,evolvermt)
end

function evolvermt.__tostring(self)
  local s = {}
  tostringRecurse(s, self, {self.field,self.action})
  if s.done then return "" end
  s[#s+1] = self.name
  s[#s+1] = " = Evolver {\n"
  s[#s+1] = string.format("  kind = \"%s\",\n",self.kind)
  s[#s+1] = string.format("  style = \"%s\",\n",self.style)
  s[#s+1] = string.format("  tau = %g,\n",self.tau)
  s[#s+1] = string.format("  nSteps = %g,\n",self.nSteps)
  s[#s+1] = string.format("  field = %s,\n",self.field.name)
  s[#s+1] = string.format("  action = %s\n",self.action.name)
  s[#s+1] = "}\n"
  return table.concat(s)
end

function evolvermt.Clone(self)
  local int = tableCopy(self)
  return setmetatable(int,evolvermt)
end

function evolvermt.Set(self, opts)
  if type(opts)=="table" then
    local t = opts._type
    if t=="Evolver" then
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

function evolvermt.Reset(self, opts)
  if self.kind == "md" then
    self.step = 0
    self.patStep = 0
  end
end

-- fs.t - time for next update
-- fs.fieldIndex - fields that need to be updated to fs.t
-- fs.actions - actions whose force needs calculating
-- fs.actionFieldIndex - which field of the action force to calculate
-- fs.eps - step for this action/force
function evolvermt.GetNextForceStruct(self)
  -- step is how many steps done so far
  local step = self.step
  local patStep = self.patStep
  local pat = self.pattern
  local eps = 0
  repeat
    patStep = patStep + 1
    if patStep > #pat.momentumStep then
      if step >= self.nSteps-1 then patStep=patStep-1; break end
      patStep = 1
      step = step + 1
    end
    eps = eps + pat.momentumStep[patStep]
  until eps ~= 0 and pat.fieldStep[patStep] ~= 0
  if eps==0 then return nil end
  --printf("step: %i\tpatStep: %i\n", step, patStep)
  local fs = {}
  local dt = self.tau/self.nSteps
  fs.t = (step+pat.fieldTime[patStep])*dt
  fs.fieldIndex = {1}
  fs.actions = {self.action}
  fs.actionFieldIndex = {1}
  fs.eps = {eps*dt}
  self.step = step
  self.patStep = patStep
  return fs
end

local function runmd(self, opts)
  local t = 0
  local tau = self.tau
  local fields = self.fields
  local momenta = self.momenta
  local fieldTimes = {}
  for i=1,#fields do fieldTimes[i] = 0 end
  --self:reset()
  while(t<tau) do
    local fs = self:GetNextForceStruct()
    if self.stepHook then self.stepHook(t, fs, self) end
    if not fs then break end
    --myprint(fs,"\n")
    -- update fields to time fs.t
    for i=1,#fs.fieldIndex do
      local fi = fs.fieldIndex[i]
      local eps = fs.t - fieldTimes[fi]
      fields[fi]:Update(momenta[fi], eps)
      fieldTimes[fi] = fs.t
    end
    for i=1,#fs.actions do
      local fi = fs.actionFieldIndex[i]
      fs.actions[i]:GetForce(momenta[fi], fields[fi], fs.eps[i])
    end
    t = fs.t
  end
  for i=1,#fields do
    local eps = tau - fieldTimes[i]
    fields[i]:Update(momenta[i], eps)
  end
  if self.stepHook then self.stepHook(tau, nil, self) end
end

function evolvermt.Run(self, opts)
  local t0 = clock()
  if self.kind == "md" then
    runmd(self, opts)
  elseif self.kind == "mc" then
    for i=1,#self.actions do
      self.actions[i]:Refresh()
    end
    local act0 = 0
    for i=1,#self.actions do
      act0 = act0 + self.actions[i]:Action()
    end
    local fSave = {}
    for i=1,#self.fields do
      fSave[i] = self.fields[i]:CopyField()
    end
    self.markov:Reset()
    self.markov:Run()
    local act1 = 0
    for i=1,#self.actions do
      act1 = act1 + self.actions[i]:Action()
    end
    local ds = act1 - act0
    local r = self.globalRand()
    local p = math.exp(-ds)
    if r > p  then -- reject
      for i=1,#self.fields do
	self.fields[i]:SetField(fSave[i])
      end
    end
    pushArray(self, "oldActions", act0)
    pushArray(self, "newActions", act1)
    pushArray(self, "mcRand", r)
  end
  updateStats(self, "run", {seconds=(clock()-t0)})
end

function getIntPat(self)
  local ipfunc = intpat[self.style]
  if not ipfunc then
    abort("unknown integrator type %s\n", self.style)
  end
  local pat = ipfunc(self)
  self.pattern = pat
  pat.fieldTime = {}
  local t = 0
  for i=1,#pat.fieldStep do
    t = t + pat.fieldStep[i]
    pat.fieldTime[i] = t
  end
end

function intpat.leapfrog(opts)
  local ip = {}
  ip.fieldStep    = { 0.5, 0.5 }
  ip.momentumStep = {   1,   0 }
  return ip
end
function intpat.omelyan(opts)
  local lambda = opts.lambda or 0.1932
  local s0 = lambda
  local s1 = (1-2*s0)
  local ip = {}
  ip.fieldStep    = {  s0,  s1, s0 }
  ip.momentumStep = { 0.5, 0.5,  0 }
  return ip
end
intpat["2MNV"] =
  function(opts)
    local lambda = opts.lambda or 0.1932
    local s0 = lambda
    local s1 = (1-2*s0)
    local ip = {}
    ip.fieldstep    = {  0, 0.5, 0.5 }
    ip.momentumstep = { s0,  s1,  s0 }
    return ip
  end
