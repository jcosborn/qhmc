-- this file is loaded at startup by default

function printfAll(...)
  io.write(string.format(...))
end

function eprintfAll(...)
  io.stderr:write(string.format(...))
end

if qopqdp and qopqdp.master then
  function printf(...)
    if(qopqdp.master()) then
      printfAll(...)
    end
  end
  function eprintf(...)
    if(qopqdp.master()) then
      eprintfAll(...)
    end
  end
else
  printf = printfAll
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

function abort(...)
  printf(...)
  printf(debug.traceback())
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

Class = require 'pl.class'
