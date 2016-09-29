require 'smear'
require 'cgms2'

local vecmt = {}
local function vector(v)
  return setmetatable(v,vecmt)
end
function vecmt.__index(v, key)
  local result = vecmt[key]
  if result then return result end
  local n = #v
  return function(...)
	   local na = select('#',...)
	   local args = { ... }
	   local r = {}
	   for iv=1,n do
	     local a = {}
	     for ia=1,na do
	       if getmetatable(args[ia])==vecmt then a[ia] = args[ia][iv]
	       else a[ia] = args[ia] end
	     end
	     r[iv] = v[iv][key](unpack(a))
	   end
	   return vector(r)
	 end
end
function vecmt.sum(v)
  local r = v[1]
  for i=2,#v do r = r + v[i] end
  return r
end
function vecmt.norm2(v)
  local r = v[1]:norm2()
  for i=2,#v do r = r + v[i]:norm2() end
  return r
end
function vecmt.reDot(v, a)
  local r = v[1]:reDot(a[1])
  for i=2,#v do r = r + v[i]:reDot(a[i]) end
  return r
end
function vecmt.dot(v, a)
  local r = v[1]:dot(a[1])
  for i=2,#v do r = r + v[i]:dot(a[i]) end
  return r
end
function vecmt.combine(v, a, c)
  for i=1,#v do
    local ai = {}
    for j=1,#a do
      ai[j] = a[j][i]
    end
    v[i]:combine(ai, c)
  end
end
function vecmt:nc()
  return self[1]:nc()
end

local actmt = {}
actmt.__index = actmt

local function getqt(a, i)
  if not a.qt[i] then
    a.qt[i] = a:quark()
  end
  return a.qt[i]
end

local setGR, setFA, setMD -- definitions moved down near functions using them

function wilson2fact(ga, params)
  local rhmc = params.rhmc
  local smear = params.smear
  local a = setmetatable({}, actmt)
  a.ga = ga
  a.coeffs = params.coeffs or {}
  a.coeffs.clov_s = a.coeffs.clov_s or 0
  a.coeffs.clov_t = a.coeffs.clov_t or 0
  a.coeffs.aniso = a.coeffs.aniso or 1
  a.w = qopqdp.wilson()
  a.qt = {}
  a.ncg = 0
  a.smear = smear
  if rhmc then
    a.npseudo = #rhmc 
    a.pseudo = {}
    for i=1,a.npseudo do a.pseudo[i] = a:quark() end
    a.ff = ga:forceNew()
    a.nff = a.npseudo
    a.rhmc = rhmc
    setGR(a)
    setFA(a)
    setMD(a)
  end
  actmt.clearStats(a)
  a.gnupdate = -1
  return a
end

function actmt.quark(a)
  --return vector{ a.w:quark() }
  return vector{ a.w:quark(), a.w:quark() }
end

function actmt.clearStats(a)
  a.LLtime = 0
  a.LLflops = 0
  a.LLn = 0
  if not a.FFtime then a.FFtime = {} end
  if not a.FFflops then a.FFflops = {} end
  --if not a.FFn then a.FFn = {} end
  --if not a.FFnorm2 then a.FFnorm2 = {} end
  --if not a.FFmax then a.FFmax = {} end
  a.FFn = {}
  a.FFnorm2 = {}
  a.FFmax = {}
  for i=1,a.nff do
    a.FFtime[i] = 0
    a.FFflops[i] = 0
    a.FFn[i] = 0
    a.FFnorm2[i] = 0
    a.FFmax[i] = 0
  end
  if not a.CGtime then a.CGtime = {} end
  if not a.CGflops then a.CGflops = {} end
  if not a.CGits then a.CGits = {} end
  if not a.CGmaxits then a.CGmaxits = {} end
  if not a.CGresid then a.CGresid = {} end
  if not a.CGmaxresid then a.CGmaxresid = {} end
  if not a.CGn then a.CGn = {} end
  for i=1,a.ncg do
    a.CGtime[i] = 0
    a.CGflops[i] = 0
    a.CGits[i] = 0
    a.CGmaxits[i] = 0
    a.CGresid[i] = 0
    a.CGmaxresid[i] = 0
    a.CGn[i] = 0
  end
end

function actmt.updateStats(a)
  a.LLmflops = 1e-6 * a.LLflops / a.LLtime
  if not a.FFmflops then a.FFmflops = {} end
  --if not a.FFrms then a.FFrms = {} end
  a.FFrms = {}
  for i=1,a.nff do
    if a.FFn[i]>0 then
      a.FFmflops[i] = 1e-6 * a.FFflops[i] / a.FFtime[i]
      --a.FFrms[i] = math.sqrt(a.FFnorm2[i]/(4*a.FFn[i]*a.ga.vol))
    else
      a.FFmflops[i] = 0
      --a.FFrms[i] = 0
    end
  end
  for k in pairs(a.FFn) do
    if a.FFn[k]>0 then
      a.FFrms[k] = math.sqrt(a.FFnorm2[k]/(4*a.FFn[k]*a.ga.vol))
    else
      a.FFrms[k] = 0
    end
  end
  if not a.CGmflops then a.CGmflops = {} end
  for i=1,a.ncg do
    if a.CGn[i]>0 then
      a.CGmflops[i] = 1e-6 * a.CGflops[i] / a.CGtime[i]
    else
      a.CGmflops[i] = 0
    end
  end
end

function actmt.set(a, g, prec)
  local nup = g:nupdates()
  if nup~=a.gnupdate then
    local t0 = clock()
    local sg = smearGauge(g, a.smear)
    a.w:set(sg, a.coeffs, prec)
    if mgSetup then mgSetup(a) end
    a.LLtime = a.LLtime + clock() - t0
    a.LLflops = a.LLflops + a.w:flops()
    a.LLn = a.LLn + 1
    a.gnupdate = nup
  end
end

function actmt.D(a, dest, src, mass)
  for i=1,#dest do
    a.w:D(dest[i], src[i], mass)
  end
end

function actmt.Ddag(a, dest, src, mass)
  for i=1,#dest do
    a.w:Ddag(dest[i], src[i], mass)
  end
end

function actmt.precD(a, dest, src, mass)
  for i=1,#dest do
    a.w:precD(dest[i], src[i], mass)
  end
end

function actmt.precDdag(a, dest, src, mass)
  for i=1,#dest do
    a.w:precDdag(dest[i], src[i], mass)
  end
end


-- multishift solver
-- dest[i] = [D D^+ + shifts[i]]^-1 src
local function solveNEshift(a, dest, src, mass, res, shifts, opts, n)
  local t = src[i]:clone()
  a.w:solve(t, src[i], mass, res, "all", opts)
  dest[i]:gamma(15, t)
  a.w:solve(t, dest[i], mass, res, "all", opts)
  dest[i]:gamma(15, t)
end

function actmt.solve(a, dest, src, mass, res, sub, opts, n, xtra)
  if sub=="NEshift" then
    local t0 = clock()
    local shifts = opts
    opts = n
    n = xtra
    local t = src:clone()
    local function op(d, s)
      a:Ddag(t, s, mass)
      a:D(d, t, mass)
    end
    local info = cgms(dest, src, op, shifts, res, opts)
    local resid = math.sqrt(info.rsq)
    if resid > res then
      printf("warning resid: %g > %g  ratio: %g  its: %i\n",
  	     resid, res, resid/res, a.w:its())
    end
    if n>0 then
      a.CGtime[n] = a.CGtime[n] + clock() - t0
      a.CGflops[n] = a.CGflops[n] + info.flops
      a.CGits[n] = a.CGits[n] + info.its
      a.CGmaxits[n] = math.max(a.CGmaxits[n], info.its)
      a.CGresid[n] = a.CGresid[n] + resid
      a.CGmaxresid[n] = math.max(a.CGmaxresid[n], resid)
      a.CGn[n] = a.CGn[n] + 1
    end
    return
  end
  for i=1,#dest do
    local t0 = clock()
    if sub=="NE" or sub=="NEshift" then
      local t = src[i]:clone()
      a.w:solve(t, src[i], mass, res, "all", opts)
      dest[i]:gamma(15, t)
      a.w:solve(t, dest[i], mass, res, "all", opts)
      dest[i]:gamma(15, t)
    else
      a.w:solve(dest[i], src[i], mass, res, sub, opts)
    end
    local resid = math.sqrt(a.w:rsq())
    if resid > res then
      printf("warning resid: %g > %g  ratio: %g  its: %i\n",
  	     resid, res, resid/res, a.w:its())
    end
    if n>0 then
      a.CGtime[n] = a.CGtime[n] + clock() - t0
      a.CGflops[n] = a.CGflops[n] + a.w:flops()
      a.CGits[n] = a.CGits[n] + a.w:its()
      a.CGmaxits[n] = math.max(a.CGmaxits[n], a.w:its())
      a.CGresid[n] = a.CGresid[n] + resid
      a.CGmaxresid[n] = math.max(a.CGmaxresid[n], resid)
    end
  end
  if n>0 then a.CGn[n] = a.CGn[n] + 1 end
end

-- Gaussian random vectors to set up the pseudofermion field
-- needs: GR{mass, shifts{}, pt{}, pt2{}, qt, coeffs{}}
function setGR(a)
  for i,r in ipairs(a.rhmc) do -- loop over pseudofermions
    local gr = r.GR
    local dq2 = gr.diquark2 or 0 -- diquark mass squared
    -- expect gr.mass, gr[1..nterms]{c,sigma} already set
    local qti = 1
    gr.qt = getqt(a,qti); qti = qti + 1
    gr.pt = {}
    gr.pt2 = {}
    gr.shifts = {}
    gr.coeffs = {}
    local nshifts = 0
    for j,t in ipairs(gr) do -- loop over rational function terms
      gr.coeffs[j] = t[1]
      if t[2] then
	nshifts = nshifts + 1
	gr.shifts[nshifts] = t[2] + dq2
	gr.pt[nshifts] = getqt(a,qti); qti = qti + 1
	gr.pt2[j] = gr.pt[nshifts]
      else
	gr.pt2[j] = gr.qt
      end
    end
    gr.resid = 1e-5 or gr.resid
    a.ncg = a.ncg + 1; gr.cgnum = a.ncg
  end
end

-- phi = [ c_0 + c_i (D D^+ + sigma_i)^-1 ] eta
-- needs: GR{mass, shifts{}, pt{}, pt2{}, qt, coeffs{}}
function actmt.refresh(a, g)
  a:set(g, 2)
  for i,r in ipairs(a.rhmc) do
    local gr = r.GR
    --gr.qt:random("all")
    gr.qt:random(math.sqrt(0.5), "all")
    a:solve(gr.pt, gr.qt, gr.mass, gr.resid, "NEshift", gr.shifts, gr.solveopts, gr.cgnum)
    a.pseudo[i]:combine(gr.pt2, gr.coeffs)
  end
end

-- Set up fermion action
-- needs: FA{mass, shifts{}, pt{}, pt2{}, qt, coeffs{}}
function setFA(a)
  for i,r in ipairs(a.rhmc) do -- loop over pseudofermions
    local fa = r.FA
    local dq2 = fa.diquark2 or 0 -- diquark mass squared
    -- expect fa.mass, fa[1..nterms]{c,sigma} already set
    local qti = 1
    fa.pt = {}
    fa.pt2 = {}
    fa.qt = getqt(a,qti); qti = qti + 1
    fa.shifts = {}
    fa.coeffs = {}
    local nshifts = 0
    for j,t in ipairs(fa) do -- loop over rational function terms
      fa.coeffs[j] = t[1]
      if t[2] then
	nshifts = nshifts + 1
	fa.shifts[nshifts] = t[2] + dq2
	fa.pt[nshifts] = getqt(a,qti); qti = qti + 1
	fa.pt2[j] = fa.pt[nshifts]
      else
	fa.pt2[j] = a.pseudo[i]
      end
    end
    fa.resid = 1e-5 or fa.resid
    a.ncg = a.ncg + 1; fa.cgnum = a.ncg
  end
end

-- q^+ q
-- q = [ c_0 + c_i (D D^+ + sigma_i)^-1 ] phi
-- needs: FA{mass, shifts{}, pt{}, pt2{}, qt, coeffs{}}
function actmt.action(a, g)
  a:set(g, 2)
  local act = 0
  for i,r in ipairs(a.rhmc) do
    local fa = r.FA
    a:solve(fa.pt, a.pseudo[i], fa.mass, fa.resid, "NEshift", fa.shifts, fa.solveopts, fa.cgnum)
    fa.qt:combine(fa.pt2, fa.coeffs)
    act = act + fa.qt:norm2()
  end
  return act
end

function actmt.force(a, f, pt, qt, ms, c, opts)
  for i=1,#pt[1] do
    local pti,qti = {},{}
    for j=1,#pt do
      pti[j] = pt[j][i]
      qti[j] = qt[j][i]
    end
    a.w:force(f, pti, qti, ms, c, opts)
  end
end

-- set up the force calculation
-- needs: MD{mass, shifts{}, pt{}, qt{}, coeffs{}}
function setMD(a)
  local qti = 1
  local nshifts = 0
  for i,r in ipairs(a.rhmc) do
    local md = r.MD
    local dq2 = md.diquark2 or 0 -- diquark mass squared
    md.ff = a.ff
    md.shifts = {}
    md.coeffs = {}
    md.pt = {}
    md.qt = {}
    for j,t in ipairs(md) do -- loop over rational function terms
      if t[2] then
	nshifts = nshifts + 1
	md.coeffs[nshifts] = t[1]
	md.shifts[nshifts] = t[2] + dq2
	md.pt[nshifts] = getqt(a,qti); qti = qti + 1
	md.qt[nshifts] = getqt(a,qti); qti = qti + 1
      end
    end
    md.resid = 1e-5 or md.resid
    a.ncg = a.ncg + 1; md.cgnum = a.ncg
  end
end

-- (d/dU) c_i phi^+ [D D^+ + sigma_i]^-1 phi
-- c_i p^+ dD q
-- q = [D D^+ + sigma_i]^-1 phi
-- p = D^+ q
-- needs: MD{mass, shifts{}, pt{}, qt{}, coeffs{}}
function actmt.updateMomentum(a, f, g, teps, ti)
  if type(teps) ~= "table" then teps = rep(teps, a.npseudo) end
  if not ti then
    ti = {}
    for k=1,a.npseudo do ti[#ti+1] = k end
  end
  a:set(g, 2)

  local pt,qt,ms,c = {},{},{},{}
  local imin = a.npseudo
  for k,i in ipairs(ti) do
    local t = a.rhmc[i].MD
    -- t.pt = (M M^+ + sigma)^{-1} phi 
    a:solve(t.pt, a.pseudo[i], t.mass, t.resid, "NEshift", t.shifts, t.solveopts, t.cgnum)
    for j=1,#t.pt do
      a:Ddag(t.qt[j], t.pt[j], t.mass) -- t.qt = M^+ Y
      pt[#pt+1] = t.pt[j]
      qt[#qt+1] = t.qt[j]
      ms[#ms+1] = t.mass
      c[#c+1] = teps[k]*t.coeffs[j]
    end
    if i<imin then imin = i end
  end
  local tmin = a.rhmc[imin].MD

  local t0 = clock()
  tmin.ff.f:zero()
  if a.smear and #a.smear>0 then
    --a:force(tmin.ff.f, pt, qt, ms, c, {prec=tmin.ffprec,deriv=1})
    a:force(tmin.ff.f, pt, qt, ms, c, {prec=tmin.ffprec,deriv=1,all=1})
    smearForce(tmin.ff.f, g, a.smear)
  else
    --a:force(tmin.ff.f, pt, qt, ms, c, {prec=tmin.ffprec})
    a:force(tmin.ff.f, pt, qt, ms, c, {prec=tmin.ffprec,all=1})
    --a.w:force(tmin.ff.f, pt, qt, ms, c, {prec=tmin.ffprec,deriv=1})
    --tmin.ff.f:derivForce(g.g)
  end
  a.FFtime[imin] = a.FFtime[imin] + clock() - t0
  --a.FFtime[imin] = a.FFtime[imin] + a.ff.f:time()
  a.FFflops[imin] = a.FFflops[imin] + a.ff.f:flops()
  a.FFn[imin] = a.FFn[imin] + 1
  local ff2 = tmin.ff.f:norm2()
  local ffi = tmin.ff.f:infnorm()
  if a.printforce then
    printf("fermion force: norm2 %g  inf %g\n", ff2, ffi)
  end
  a.FFnorm2[imin] = a.FFnorm2[imin] + ff2
  if a.FFmax[imin] < ffi then a.FFmax[imin] = ffi end

  table.sort(ti)
  local key = table.concat(ti, ",")
  if not a.FFn[key] then a.FFn[key] = 0 end
  a.FFn[key] = a.FFn[key] + 1
  if not a.FFnorm2[key] then a.FFnorm2[key] = 0 end
  a.FFnorm2[key] = a.FFnorm2[key] + ff2
  if not a.FFmax[key] then a.FFmax[key] = 0 end
  if a.FFmax[key] < ffi then a.FFmax[key] = ffi end

  f.f:fupdate(tmin.ff.f, 1)
end

function actmt.pbp(a, g, mass, resid, opts)
  a:set(g, 2)
  local x = getqt(a, 1)
  local y = getqt(a, 2)
  x:randomU1()
  a:solve(y, x, mass, resid, "all", opts, 0)
  return x:reDot(y)/(4*x:nc()*a.ga.vol)
end

function actmt.makeprop(a, g, mass, resid, opts)
  local nqt = 0
  local src = getqt(a, nqt); nqt=nqt+1
  local Nc = src:nc()
  local Ns = 4
  local dest = {}
  for color = 0,Nc-1 do
    for spin = 0,Ns-1 do
      src:zero()
      -- point({coord},color,spin,re,im)
      src:point({0,0,0,0},color,spin,1,0)
      printf("src norm2: %g\n", src:norm2())
      if(opts.ngauss) then src:smearGauss(g, opts.gaussalpha, opts.ngauss) end
      printf("src norm2: %g\n", src:norm2())
      dest[#dest+1] = a.w:quark()
      t0 = qopqdp.dtime()
      w:solve(dest[#dest], src, mass, resid, "all", opts)
      dt = qopqdp.dtime() - t0
      mf = 1e-6 * w:flops() / dt
      printf("its: %g  secs: %g  Mflops: %g\n", a.w:its(), dt, mf)
    end
  end
  return dest
end

function actmt.pions(prop)
  local Nc = prop[1]:nc()
  local Ns = 4
  local pions = {}
  local t = {}
  for g = 0,0 do
    pions[g] = {}
    for color = 0,Nc-1 do
      for spin = 0,Ns-1 do
	local k = color*Ns + spin
	t[k] = prop[k]:norm2("timeslices")
	for j=1,#t[k] do
	  pions[g][j] = t[k][j] + (pions[g][j] or 0)
	end
      end
    end
  end
  return pions
end
