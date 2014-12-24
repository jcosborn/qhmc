require 'smear'

local actmt = {}
actmt.__index = actmt

local function getqt(a, i)
  if not a.qt[i] then
    a.qt[i] = a.w:quark()
  end
  return a.qt[i]
end

-- Gaussian random vectors to set up the pseudofermion field
local function setGR(a)
  local qti = 1
  for i,r in ipairs(a.rhmc) do -- loop over pseudofermions
    local gr = r.GR
    gr.qt2 = getqt(a,qti); qti = qti + 1
    gr.qt = a.pseudo[i]
    gr.pt = a.pseudo[i]
    gr.mass = gr[1]
    gr.mass2 = nil
    if gr[2] then
      gr.mass2 = gr[2]
      gr.qt = getqt(a,qti); qti = qti + 1
      a.ncg = a.ncg + 1; gr.cgnum = a.ncg
      gr.resid = gr.resid or 1e-5
    end
  end
end

-- Set up fermion action
local function setFA(a)
  local qti = 1
  for i,r in ipairs(a.rhmc) do
    local fa = r.FA
    fa.qt = a.pseudo[i]
    fa.pt = getqt(a,qti); qti = qti + 1
    fa.mass = fa[1]
    fa.mass2 = nil
    if fa[2] then
      fa.mass2 = fa[2]
      fa.qt = getqt(a,qti); qti = qti + 1
    end
    fa.resid = fa.resid or 1e-5
    a.ncg = a.ncg + 1; fa.cgnum = a.ncg
  end
end

-- set up the force calculation
local function setMD(a)
  local qti = 1
  for i,r in ipairs(a.rhmc) do
    local md = r.MD
    md.ff = a.ff
    md.pt = getqt(a,qti); qti = qti + 1
    md.qt = getqt(a,qti); qti = qti + 1
    md.mass = md[1]
    md.mass2 = nil
    if md[2] then
      md.mass2 = md[2]
      md.qt = getqt(a,qti); qti = qti + 1
    end
    md.resid = md.resid or 1e-5
    a.ncg = a.ncg + 1; md.cgnum = a.ncg
  end
end

function wilsonact(ga, params)
  local rhmc = params.rhmc
  local smear = params.smear
  local a = {}
  a.ga = ga --gauge action? 
  a.coeffs = params.coeffs or {}
  a.coeffs.clov_s = a.coeffs.clov_s or 0
  a.coeffs.clov_t = a.coeffs.clov_t or 0
  a.coeffs.aniso = a.coeffs.aniso or 1
  a.w = qopqdp.wilson() -- fermion action?
  a.npseudo = #rhmc 
  a.pseudo = {}
  for i=1,a.npseudo do a.pseudo[i] = a.w:quark() end
  a.ff = ga:forceNew()
  a.qt = {}
  a.ncg = 0
  a.nff = a.npseudo
  a.rhmc = rhmc
  a.smear = smear
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
    a.w:set(sg, a.coeffs, prec)
    if mgSetup then mgSetup(a) end
    a.LLtime = a.LLtime + clock() - t0
    a.LLflops = a.LLflops + a.w:flops()
    a.LLn = a.LLn + 1
    a.gnupdate = nup
  end
end

function actmt.solve(a, dest, src, mass, res, sub, opts, n)
  local t0 = clock()
  dest:zero()
  a.w:solve(dest, src, mass, res, sub, opts)
  if n>0 then
    a.CGtime[n] = a.CGtime[n] + clock() - t0
    a.CGflops[n] = a.CGflops[n] + a.w:flops()
    a.CGits[n] = a.CGits[n] + a.w:its()
    a.CGmaxits[n] = math.max(a.CGmaxits[n], a.w:its())
    a.CGn[n] = a.CGn[n] + 1
  end
end

function actmt.refresh(a, g)
  a:set(g, 2)
  for i,r in ipairs(a.rhmc) do
    local t = r.GR
    --t.qt2:random("even") -- t.qt2 holds the random Gaussian vector
    t.qt2:random(math.sqrt(0.5), "even") -- t.qt2 holds random Gaussian vector
    --printf("t.qt2 norm2: %g\n", t.qt2:norm2("even"))
    a.w:precD(t.qt, t.qt2, t.mass) -- t.qt holds the pseudofermion vector
    --printf("t.qt norm2: %g\n", t.qt:norm2("even"))
    if t.mass2 then
      a:solve(t.pt, t.qt, t.mass2, t.resid, "prec", t.solveopts, t.cgnum)
    --[[ not sure if this is needed
    else 
      t.pt = t.qt
    -- not sure ended ]]
    end
  end
end

function actmt.action(a, g)
  a:set(g, 2)
  local act = 0
  for i,r in ipairs(a.rhmc) do
    local fa = r.FA
    if fa.mass2 then
      a.w:precD(fa.qt, a.pseudo[i], fa.mass2) --fa.qt: fermion action quark temporary
    end
    -- fa.pt holds the \chi vector 
    a:solve(fa.pt, fa.qt, fa.mass, fa.resid, "prec", fa.solveopts, fa.cgnum)
    act = act + fa.pt:norm2("even")
  end
  -- print("act = ", act)
  return act
end

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
    if t.mass2 then
       a.w:precD(t.qt, a.pseudo[i], t.mass2) -- tempphi is now W \phi
       -- (M M^+)^{-1} W \phi = Y
       a:solve(t.pt, t.qt, t.mass, t.resid, "precNE", t.solveopts, t.cgnum)	
    else
       -- t.pt = (M M^+)^{-1} \phi 
       a:solve(t.pt, a.pseudo[i], t.mass, t.resid, "precNE", t.solveopts, t.cgnum)
    end
    a.w:precDdag(t.qt, t.pt, t.mass) -- t.qt = M^+ Y

    if t.mass2 then
      local kappa = 0.5/(t.mass+1.0+3.0*a.coeffs.aniso)
      local kappa2 = 0.5/(t.mass2+1.0+3.0*a.coeffs.aniso)
      local coeffs = {1, -kappa2*kappa2/(kappa*kappa)}
      t.qt:combine({t.qt, a.pseudo[i]}, coeffs, "even") -- M^+Y - \tilde{k}^2/k^2 \phi
    end

    pt[#pt+1] = t.pt
    qt[#qt+1] = t.qt
    ms[#ms+1] = t.mass
    c[#c+1] = teps[k] -- integration step size
    if i<imin then imin = i end
  end
  local tmin = a.rhmc[imin].MD

  local t0 = clock()
  tmin.ff.f:zero()
  if a.smear then
    a.w:force(tmin.ff.f, pt, qt, ms, c, {prec=tmin.ffprec,deriv=1})
    smearForce(tmin.ff.f, g, a.smear)
  else
    a.w:force(tmin.ff.f, pt, qt, ms, c, {prec=tmin.ffprec})
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
