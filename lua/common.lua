function printf(...)
  io.write(string.format(...))
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
	  myprintr(x[i], indent, copy(seen))
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
	  myprintr(v, indent.."  ", copy(seen))
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

function readfile(fn, levelkeys)
  local p = {}
  local pl = {}
  pl[0] = p
  local c = p
  local cw
  local first
  local function insert(word)
    if first then
      first = false
      local lk = levelkeys[word]
      if lk then -- change level
	if lk.l~=0 then
	  local x = pl[lk.l-1]
	  local xt = x[lk.t]
	  if not xt then
	    xt = {}
	    x[lk.t] = xt
	  end
	  pl[lk.l] = {}
	  xt[#xt+1] = pl[lk.l]
	end
	c = pl[lk.l]
      end
      cw = c[word]
      if not cw then
	cw = {}
	c[word] = cw
      end
    else
      cw[#cw+1] = word
    end
  end
  for line in io.lines(fn) do
    if line:sub(1,1)~='#' then
      first = true
      line:gsub("%S+", insert)
    end
  end
  return p
end

function rep(x, n)
  local t = {}
  for i=1,n do t[i] = x end
  return t
end

function repelem(t, n)
  local r = {}
  for i=1,#t do
    for j=1,n do
      r[#r+1] = t[i]
    end
  end
  return r
end

function paste(t1,t2)
  local t = {}
  for i=1,#t1 do
    t[i] = { t1[i], t2[i] }
  end
  return t
end

function sum(t)
  local s = 0
  for k,v in pairs(t) do s = s + v end
  return s
end

local function copyr(t, seen)
  if type(t) ~= "table" then return t end
  if seen[t] then return seen[t] end
  local r = {}
  seen[t] = r
  for k,v in pairs(t) do
    r[k] = seen[v] or copy(v)
  end
  return r
end
function copy(t)
  return copyr(t, {})
end

local function copytor(d, s, seen)
  if seen[s] then return end
  seen[s] = true
  for k,v in pairs(s) do
    if type(v)=="table" and type(d[k])=="table" then
      copytor(d[k], v, seen)
    else
      d[k] = copyr(v, seen)
    end
  end
end
function copyto(d, s)
  copytor(d, s, {})
end

local function wrap(r)
  return { __index =
	   function(t,k)
	     local v = rawget(t,k)
	     if not v then
	       v = r[k]
	       if type(v) == "table" then
		 v = setmetatable({}, wrap(v))
		 rawset(t,k,v)
	       end
	     end
	     return v
	   end,
	 __newindex = 
	   function(t,k,v)
	     if r[k] then r[k]=v
	     else rawset(t,k,v) end
	   end }
end


if stdoutfile then
  stdoutfile = stdoutfile:gsub("%%job",jobnum)
  qopqdp.remapout(stdoutfile)
end

-- translate args
if arg then
  for i=1,#arg do
    arg[i] = arg[i]:gsub("%%job",jobnum)
  end
end
