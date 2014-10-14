function clock()
  --return os.clock()
  --return os.time()
  return qopqdp.dtime()
end

function profcontrol(...)
  return qopqdp.profile(...)
end

function verbcontrol(...)
  return qopqdp.verbosity(...)
end

function clearStats(object, name)
  if not object.stats then object.stats={} end
  if not object.stats[name] then object.stats[name]={} end
  local p = object.stats[name]
  for k,v in pairs(p) do p[k] = 0 end
end

function updateStats(object, name, stats)
  local p = object.stats[name]
  for k,v in pairs(stats) do
    if k:match("max") then
      p[k] = math.max((p[k] or 0), v)
    else
      p[k] = (p[k] or 0) + v
    end
  end
end

function pushArray(base, name, value)
  if base==nil then base = _G end
  local t = base[name]
  if t==nil then t = {}; base[name] = t end
  t[#t+1] = value
end

function tostringRecurse(s, self, others)
  if self.tostringRecurse then
    self.tostringRecurseCounter = (self.tostringRecurseCounter or 0) + 1
    if self.tostringRecurseCounter>1 then
      s.done = true
      return
    end
    local r = {}
    for i=1,#others do
      r[i] = others[i].tostringRecurse
      others[i].tostringRecurse = true
      others[i].tostringRecurseCounter = 0
    end
    for i=1,#others do
      s[#s+1] = tostring(others[i])
    end
    for i=1,#others do
      others[i].tostringRecurse = r[i]
    end
  end
end
