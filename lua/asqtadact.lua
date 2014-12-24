require 'smear'
require 'staggeredObservables'

local actmt = {}
actmt.__index = actmt

local function getqt(a, i)
  if not a.qt[i] then
    a.qt[i] = a.h:quark()
    a.qt[i]:zero()
  end
  return a.qt[i]
end

local function setGR(a)
  local qti = 1
  for i,r in ipairs(a.rhmc) do -- loop over pseudofermions
    local gr = r.GR
    gr.pt0 = {}
    gr.coeffs = {}
    for j,t in ipairs(r.GR) do -- loop over random sources
      local fac = 1
      if t.allfacodd and t.allfacodd ~= 0 then
	t.qt2 = getqt(a,qti); qti = qti + 1
	fac = 2*t.allfacodd
	t.allmass = t.allfaceven / fac
      else
	fac = t.allfaceven or fac
      end
      t.qt = getqt(a,qti); qti = qti + 1
      t.pt = {}
      t.masses = {}
      for k=1,#t do
	if(t[k][2]) then
	  t.pt[#t.pt+1] = getqt(a,qti); qti = qti + 1
	  gr.pt0[#gr.pt0+1] = t.pt[#t.pt]
	  t.masses[#t.masses+1] = math.sqrt(0.25*t[k][2])
	  gr.coeffs[#gr.coeffs+1] = fac*t[k][1]/(4*t.masses[#t.masses])
	else
	  gr.pt0[#gr.pt0+1] = t.qt
	  gr.coeffs[#gr.coeffs+1] = fac*t[k][1]
	end
      end
      if #t.pt == 0 then
	t.pt = nil
	t.masses = nil
      else
	t.resid = t.resid or 1e-5
	a.ncg = a.ncg + 1; t.cgnum = a.ncg
      end
    end
  end
end

local function setFA(a)
  local qti = 1
  for i,r in ipairs(a.rhmc) do
    local fa = r.FA
    fa.pt = {}
    fa.masses = {}
    for j,t in ipairs(fa) do
      local fac = 1
      if t.allfacodd and t.allfacodd ~= 0 then
	t.qt2 = getqt(a,qti); qti = qti + 1
	fac = 2*t.allfacodd
	t.allmass = t.allfaceven / fac
      else
	fac = t.allfaceven or fac
      end
      t.qt = getqt(a,qti); qti = qti + 1
      t.pt0 = {}
      t.coeffs = {}
      for k=1,#t do
	if(t[k][2]) then
	   fa.pt[#fa.pt+1] = getqt(a,qti); qti = qti + 1
	   t.pt0[#t.pt0+1] = fa.pt[#fa.pt]
	   fa.masses[#fa.masses+1] = math.sqrt(0.25*t[k][2])
	   t.coeffs[#t.coeffs+1] = fac*t[k][1]/(4*fa.masses[#fa.masses])
	else
	   t.pt0[#t.pt0+1] = a.pseudo[i]
	   t.coeffs[#t.coeffs+1] = fac*t[k][1]
	end
      end
      printf("%i %i : %s\n", i, j, tostring(t.resid))
      fa.resid = fa.resid or 1e-5
    end
    if #fa.pt == 0 then
      fa.pt = nil
      fa.masses = nil
    else
      a.ncg = a.ncg + 1; fa.cgnum = a.ncg
    end
  end
end

local function setMD(a)
  local qti = 1
  for i=1,a.npseudo do
    local t = a.rhmc[i].MD
    t.resid = t.resid or 1e-5
    a.ncg = a.ncg + 1; t.cgnum = a.ncg
    t.ff = a.ff
    t.pt = {}
    t.masses = {}
    t.coeffs = {}
    for j=1,#t do
      if(t[j][2]) then
	t.pt[#t.pt+1] = getqt(a,qti); qti = qti + 1
	local m = math.sqrt(0.25*t[j][2])
	t.masses[#t.masses+1] = m
	t.coeffs[#t.coeffs+1] = t[j][1]/(8*m*m)
      end
    end
  end
end

function asqtad_coeffs(u0)
  local u2 = 1/(u0*u0);
  local u4 = u2*u2;
  return {
    one_link = 5.0/8.0,
    three_staple = -u2/16.0,
    five_staple = u4/64.0,
    seven_staple = -(u4*u2)/384.0,
    lepage = -u4/16.0,
    naik = -u2/24.0
  }
end

function asqtadact(ga, params)
  local a = {}
  a.ga = ga
  a.h = qopqdp.asqtad()
  a.smear = params.smear
  a.u0 = params.u0
  if params.coeffs then
    a.h:coeffs(params.coeffs)
  else
    a.h:coeffs(asqtad_coeffs(params.u0))
  end
  local cfs = a.h:coeffs()
  myprint("asqtad_coeffs = ", cfs, "\n")
  a.npseudo = #params.rhmc
  a.pseudo = {}
  for i=1,a.npseudo do a.pseudo[i] = a.h:quark() end
  a.ff = ga:forceNew()
  a.qt = {}
  a.ncg = 0
  a.nff = a.npseudo
  a.rhmc = params.rhmc
  setGR(a)
  setFA(a)
  setMD(a)
  actmt.clearStats(a)
  a.gnupdate = -1
  return setmetatable(a, actmt)
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
  if not a.CGn then a.CGn = {} end
  for i=1,a.ncg do
    a.CGtime[i] = 0
    a.CGflops[i] = 0
    a.CGits[i] = 0
    a.CGmaxits[i] = 0
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
    a.h:set(sg, prec)
    a.LLtime = a.LLtime + clock() - t0
    a.LLflops = a.LLflops + a.h:flops()
    a.LLn = a.LLn + 1
    a.gnupdate = nup
  end
end

function actmt.solve(a, dest, src, m, res, sub, opts, n)
  local t0 = clock()
  a.h:solve(dest, src, m, res, sub, opts)
  if a.h:rsq() > res*res then
    printf("warning rsq: %g > %g  ratio: %g  its: %i\n",
	   a.h:rsq(), res*res, a.h:rsq()/(res*res), a.h:its())
  end
  if n>0 then
    a.CGtime[n] = a.CGtime[n] + clock() - t0
    --a.CGtime[n] = a.CGtime[n] + a.h:time()
    a.CGflops[n] = a.CGflops[n] + a.h:flops()
    a.CGits[n] = a.CGits[n] + a.h:its()
    a.CGmaxits[n] = math.max(a.CGmaxits[n], a.h:its())
    a.CGn[n] = a.CGn[n] + 1
  end
end

function actmt.refresh(a, g)
  a:set(g, 2)
  for i,r in ipairs(a.rhmc) do
    for j,t in ipairs(r.GR) do
      if t.allmass then
	if t.allmass2 then
	  t.qt:random(math.sqrt(0.5), "all")
	  a.h:Ddag(t.qt2, t.qt, t.allmass2, "all", "all")
        else
	  t.qt2:random(math.sqrt(0.5), "all")
        end
	--printf("t.qt2 norm2: %g\n", t.qt2:norm2("all"))
        a.h:D(t.qt, t.qt2, t.allmass, "even", "all")
	--printf("t.qt norm2: %g\n", t.qt:norm2("even"))
      else
	t.qt:random(math.sqrt(0.5), "even")
      end
      if t.pt then
	a:solve(t.pt, t.qt, t.masses, t.resid, "even", t.solveopts, t.cgnum)
      end
    end
    a.pseudo[i]:combine(r.GR.pt0, r.GR.coeffs, "even")
  end
end

function actmt.action(a, g)
  a:set(g, 2)
  local act = 0
  for i,r in ipairs(a.rhmc) do
    --printf("a.pseudo[%i] norm2: %g\n", i, a.pseudo[i]:norm2("even"))
    local fa = r.FA
    if fa.pt then
      a:solve(fa.pt, a.pseudo[i], fa.masses, fa.resid, "even", fa.solveopts, fa.cgnum)
      --for ii=1,#fa.pt do
      --printf("fa.pt[%i] norm2: %g\n", ii, fa.pt[ii]:norm2("even"))
      --end
    end
    for j,t in ipairs(fa) do
      t.qt:combine(t.pt0, t.coeffs, "even")
      if t.allmass then
	a.h:D(t.qt2, t.qt, t.allmass, "all", "even")
	act = act + t.qt2:norm2("all")
      else
	act = act + t.qt:norm2("even")
      end
    end
  end
  --act = a.pseudo[1]:Re_dot(a.pseudo[2]);
  return act
end

function actmt.updateMomentum(a, f, g, teps, ti)
  if type(teps) ~= "table" then teps = rep(teps, a.npseudo) end
  if not ti then
    ti = {}
    for k=1,a.npseudo do ti[#ti+1] = k end
  end
  a:set(g, 2)

  local pt,c = {},{}
  local imin = a.npseudo
  for k,i in ipairs(ti) do
    local t = a.rhmc[i].MD
    a:solve(t.pt, a.pseudo[i], t.masses, t.resid, "even", t.solveopts, t.cgnum)
    for j=1,#t.pt do
      a.h:D(t.pt[j], t.pt[j], 0.5, "all", "even")
      pt[#pt+1] = t.pt[j]
      c[#c+1] = teps[k]*t.coeffs[j]
    end
    if i<imin then imin = i end
  end
  local tmin = a.rhmc[imin].MD

  local t0 = clock()
  if a.smear then
    a.h:force(tmin.ff.f, pt, c, {prec=tmin.ffprec,deriv=1})
    smearForce(tmin.ff.f, g, a.smear)
  else
    a.h:force(tmin.ff.f, pt, c, {prec=tmin.ffprec})
    --a.h:force(tmin.ff.f, pt, c, {prec=tmin.ffprec,deriv=1})
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

------------ start observables ------------------

function actmt.pbp(a, g, mass, resid, opts)
  a:set(g, 2)
  local x = getqt(a, 1)
  local y = getqt(a, 2)
  --x:random()
  x:randomU1()
  a:solve({y}, x, {mass}, resid, "all", opts, 0)
  return mass*y:norm2()/a.ga.vol, x:norm2()/a.ga.vol
end

function actmt.pions(a, g, wallslice, mass, resid, opts)
  a:set(g, 2)
  return staggeredPions(a.h, g, wallslice, mass, resid, opts)
end

-- Fully reproduces MILC's nl_spectrum function.
function actmt.pions_wall(a, g, wallslice, mass, resid, opts)
  a:set(g, 2) -- performs the nHYP smearing.
  return staggeredPionsWall(a.h, g, wallslice, mass, resid, opts)
end

function actmt.s4_broken_observe(a, g, mass, resid, opts, npbp)
  a:set(g, 2) -- performs the nHYP smearing.
  -- Get smeared gauge now.
  local sg = smearGauge(g, a.smear)
  return staggeredS4Broken(a.h, sg, mass, resid, opts, npbp)
end
