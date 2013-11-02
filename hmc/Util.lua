function printf(...)
  io.write(string.format(...))
end

function clock()
  --return os.clock()
  --return os.time()
  return qopqdp.dtime()
end

function profile(...)
  return qopqdp.profile(...)
end

function verbosity(...)
  return qopqdp.verbosity(...)
end

do
  local spaces = ""
  local function printline(event, ln)
    if(event == "call") then 
      spaces = spaces .. "  "
      elseif(event == "return") then
      spaces = string.sub(spaces, 1, -3)
    else
      if ln==nil then ln = "" end
      print(spaces .. (debug.getinfo(2,"n").name or "main") ..  " : " .. ln)
    end
  end
  function trace(on)
    if(on) then
      debug.sethook(printline, "crl")
    else
      debug.sethook(nil, "crl")
    end
  end
end

function tableCopy(t)
  local r = {}
  for k,v in pairs(t) do r[k]=v end
  return r
end

function tableCopyTo(r, t)
  for k,v in pairs(t) do r[k]=v end
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

function abort(...)
  printf(...)
  printf(debug.traceback())
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

do
  local function simple(t)
    for k,v in pairs(t) do
      if type(v) == "table" then return false end
    end
    return true
  end
  local function numeric(t)
    for k,v in pairs(t) do
      if type(k) ~= "number" then return false end
      if k ~= math.floor(k) then return false end
      if k<1 then return false end
    end
    return true
  end
  local function sortkeys(a, b)
    local ta,tb = type(a),type(b)
    if ta==tb then return a<b end
    return ta>tb
  end
  local function myprintr(x, indent, seen)
    if type(x) == "table" then
      if seen[x] then printf("<loop>"); return end
      seen[x] = true
      printf("{")
      if simple(x) and numeric(x) then
        local n = table.maxn(x)
        printf(" ")
        for i=1,n do
          myprintr(x[i], indent, tableCopy(seen))
          if i<n then printf(", ")
          else printf(" ") end
        end
      else
        printf("\n")
        local keys = {}
        for k in pairs(x) do keys[#keys+1] = k end
        table.sort(keys, sortkeys)
        for _,k in ipairs(keys) do
          local v = x[k]
          printf("  %s%s = ", indent, tostring(k))
          myprintr(v, indent.."  ", tableCopy(seen))
          printf("\n")
        end
        if(indent~="") then printf("%s", indent) end
      end
      printf("}")
    else
      local s = tostring(x)
      if x==nil then s = "nil" end
      printf("%s", s)
    end
  end
  function myprint(...)
    local t = {...}
    local n = table.maxn(t)
    for i=1,n do
      myprintr(t[i], "", {})
    end
  end
end
