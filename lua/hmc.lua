local setupint, integrate

function hmcstep(fields, params)
  fields:save()
  local Sold = fields:action()

  if type(params.forceparams)=="function" then
    params.forceparams(fields, params)
  else
    local intparams = setupint(fields, params)
    integrate(fields, intparams)
  end

  local Snew = fields:action()

  local ds = Snew - Sold
  printf("Sold: %20.6f  Snew: %20.6f\n", Sold, Snew)

  if params.checkReverse then
    fields:reverse()
    if type(params.forceparams)=="function" then
      params.forceparams(fields, params)
    else
      local intparams = setupint(fields, params)
      integrate(fields, intparams)
    end
    local Srev = fields:action()
    printf("Sold: %20.6f Srev: %20.6f dS: %20.6f\n", Sold, Srev, Srev-Sold)
    fields:endReverse()
  end

  local r = 0
  if not params.md then
    r = globalRand()
  end
  local p = math.exp(-ds)
  if( r > p ) then -- reject
    printf("REJECT deltaS = %g  p = %g  r = %g\n", ds, p, r)
    fields:reject()
  else
    printf("ACCEPT deltaS = %g  p = %g  r = %g\n", ds, p, r)
    fields:accept()
  end
end

local function sym(a, n)
  local t,s = {},0
  local h = #a
  for i=1,h do
    s = s + a[i]
    t[i] = a[i]
  end
  local r = n - 2*h
  s = (1-2*s)/r
  for i=1,r do
    t[h+i] = s
  end
  for i=1,h do
    t[n-h+i] = a[h+1-i]
  end
  return t
end

local function fn(eps, nforce, as, bs, gs)
  local nsteps = nforce + 1
  local e = nforce*eps
  local fld = sym(as, nsteps)
  local frc = sym(bs, nforce)
  local ip = {}
  ip.nsteps = nsteps
  ip.nforcesteps = nforce
  ip.sstep = {}
  ip.fieldstep = {}
  ip.forcestep = {}
  ip.fgstep = {}
  for i=1,nsteps do
    local a = e*fld[i]
    local b = e*(frc[i] or 0)
    local g = gs[i] or gs[nsteps-i] or 0
    ip.sstep[i] = a
    ip.fieldstep[i] = a
    ip.forcestep[i] = b
    ip.fgstep[i] = (b==0 and 0) or e*e*e*g/(0.5*b)
  end
  return ip
end

local function fnv(eps, nforce, as, bs, gs)
  local nsteps = nforce + 1
  local e = nforce*eps
  local frc = sym(as, nsteps)
  local fld = sym(bs, nforce)
  local ip = {}
  ip.nsteps = nsteps
  ip.nforcesteps = nforce
  ip.sstep = {}
  ip.fieldstep = {}
  ip.forcestep = {}
  ip.fgstep = {}
  for i=1,nsteps do
    local a = e*frc[i]
    local b = e*(fld[i-1] or 0)
    local g = gs[i] or gs[nsteps+1-i] or 0
    ip.sstep[i] = b
    ip.fieldstep[i] = b
    ip.forcestep[i] = a
    ip.fgstep[i] = (a==0 and 0) or e*e*e*g/(0.5*a)
  end
  return ip
end

local intpat = {}
function intpat.f1(eps, p) -- TVT
  local g0 = p.g0 or 0
  return fn(eps, 1, {}, {}, {g0})
end
function intpat.f1v(eps, p) -- VTV
  local g0 = p.g0 or 0
  return fnv(eps, 1, {}, {}, {g0})
end
function intpat.f2(eps, p) -- TVTVT
  local a0 = p.a0 or 0.2
  local g0 = p.g0 or 0
  return fn(eps, 2, {a0}, {}, {g0})
end
function intpat.f2v(eps, p) -- VTVTV
  local a0 = p.a0 or 0.2
  local g0 = p.g0 or 0
  local g1 = p.g1 or 0
  return fnv(eps, 2, {a0}, {}, {g0,g1})
end
function intpat.f3(eps, p) -- TVTVTVT
  local a0 = p.a0 or 1/6
  local b0 = p.b0 or 3/8
  local g0 = p.g0 or 0
  local g1 = p.g1 or 0
  return fn(eps, 3, {a0}, {b0}, {g0,g1})
end
function intpat.f3v(eps, p) -- VTVTVTV
  local a0 = p.a0 or 1/6
  local b0 = p.b0 or 3/8
  local g0 = p.g0 or 0
  local g1 = p.g1 or 0
  return fnv(eps, 3, {a0}, {b0}, {g0,g1})
end
function intpat.f4(eps, p) -- TVTVTVTVT
  local a0 = p.a0 or 0.1
  local a1 = p.a1 or (0.5-2*a0)
  local b0 = p.b0 or 0.25
  local g0 = p.g0 or 0
  local g1 = p.g1 or 0
  return fn(eps, 4, {a0,a1}, {b0}, {g0,g1})
end
function intpat.f4v(eps, p) -- VTVTVTVTV
  local a0 = p.a0 or 0.1
  local a1 = p.a1 or (0.5-2*a0)
  local b0 = p.b0 or 0.25
  local g0 = p.g0 or 0
  local g1 = p.g1 or 0
  local g2 = p.g2 or 0
  return fnv(eps, 4, {a0,a1}, {b0}, {g0,g1,g2})
end
function intpat.f2g(eps, p) -- VTGTV
  local a0 = p.a0 or 1/6
  local g0 = p.g0 or 0
  local g1 = p.g1 or 1/72
  return fnv(eps, 2, {a0}, {}, {g0,g1})
end
function intpat.f2h(eps, p) -- GTVTG
  local a0 = p.a0 or 1/6
  local g0 = p.g0 or 1/144
  local g1 = p.g1 or 0
  return fnv(eps, 2, {a0}, {}, {g0,g1})
end
function intpat.f3g(eps, p) -- TVTGTVT
  local a0 = p.a0 or 0.1
  local b0 = p.b0 or 1/(6*(2*a0-1)^2)
  local g0 = p.g0 or 0
  local g1 = p.g1 or (48*a0^3-48*a0^2+12*a0-1)/(72*(2*a0-1)^3)
  return fn(eps, 3, {a0}, {b0}, {g0,g1})
end
function intpat.f3h(eps, p) -- GTVTVTG
  local b0 = p.b0 or 0.18
  local a0 = p.a0 or (6*b0^2-6*b0+1)/(12*(b0-1)*b0)
  local g0 = p.g0 or -(6*b0^3-12*b0^2+6*b0-1)/(288*(b0-1)^2*b0)
  local g1 = p.g1 or 0
  return fnv(eps, 3, {a0}, {b0}, {g0,g1})
end
function intpat.f4g(eps, p) -- VTVTGTVTV
  local a0 = p.a0 or 0.07
  local b0 = p.b0 or 0.2
  local a1 = p.a1 or (1-6*a0)/(6*(2*b0-1)^2)
  local g0 = p.g0 or 0
  local g1 = p.g1 or 0
  local g2 = p.g2 or (576*a1^2*b0^4-576*a1^2*b0^3+96*a1*b0^2+72*a1^2*b0-48*a1*b0+1)/72
  return fnv(eps, 4, {a0,a1}, {b0}, {g0,g1,g2})
end
function intpat.f4h(eps, p) -- GTVTVTVTG
  local a0 = p.a0 or 0.05
  local b0 = p.b0 or 0.2
  local a1 = p.a1 or (1-6*a0)/(6*(2*b0-1)^2)
  local g0 = p.g0 or (576*a1^2*b0^4-576*a1^2*b0^3+96*a1*b0^2+72*a1^2*b0-48*a1*b0+1)/144
  local g1 = p.g1 or 0
  local g2 = p.g2 or 0
  return fnv(eps, 4, {a0,a1}, {b0}, {g0,g1,g2})
end

function intpat.leapfrog(eps, p)
  return fn(eps, 1, {}, {}, {})
end
function intpat.omelyan(eps, p)
  local lambda = p.lambda or 0.1932
  return fn(eps, 2, {lambda}, {}, {})
end
intpat["2MNV"] = function(eps, p)
  local lambda = p.lambda or 0.1932
  return fnv(eps, 2, {lambda}, {}, {})
end

function intpat.fgv(eps, p)
  local sigma = p.sigma or 0
  local s0 = 2*eps/6
  local s1 = 2*(eps-s0)
  local e3 = 8*eps*eps*eps
  local f0 = e3*sigma
  local f1 = e3/72 - 2*f0
  f0 = f0/(0.5*s0)
  f1 = f1/(0.5*s1)
  local ip = {}
  ip.nsteps = 3
  ip.nforcesteps = 2
  ip.sstep =     {  0, eps, eps }
  ip.fieldstep = {  0, eps, eps }
  ip.forcestep = { s0,  s1,  s0 }
  ip.fgstep =    { f0,  f1,  f0 }
  return ip
end
function intpat.fgp(eps, p)
  local s0 = 2*eps*(3-math.sqrt(3))/6
  local s1 = 2*(eps-s0)
  local e3 = 8*eps*eps*eps
  local f0 = e3*(2-math.sqrt(3))/48
  f0 = f0/(0.5*eps)
  local ip = {}
  ip.nsteps = 3
  ip.nforcesteps = 2
  ip.sstep =     {  s0,  s1, s0 }
  ip.fieldstep = {  s0,  s1, s0 }
  ip.forcestep = { eps, eps,  0 }
  ip.fgstep =    {  f0,  f0,  0 }
  return ip
end
function intpat.mn4f3(eps, p)
  local e = 3*eps
  local t = math.pow(4,1/3)
  local a0 = e*(4+t*(1+t))/12
  local a1 = 0.5*e - a0
  local b0 = 2*a0
  local b1 = e - 2*b0
  local s1 = 1e-6*e
  local s0 = 0.5*e - s1
  local ip = {}
  ip.nsteps = 4
  ip.nforcesteps = 3
  ip.sstep =     { s0, s1, s1, s0 }
  ip.fieldstep = { a0, a1, a1, a0 }
  ip.forcestep = { b0, b1, b0,  0 }
  return ip
end
function intpat.mn4f3v(eps, p)
  local e = 3*eps
  local t = math.pow(4,1/3)
  local a0 = e*(4+t*(1+t))/12
  local a1 = 0.5*e - a0
  local b0 = 2*a0
  local b1 = e - 2*b0
  local s1 = 1e-6*e
  local s0 = 0.5*(e - s1)
  local ip = {}
  ip.nsteps = 4
  ip.nforcesteps = 3
  ip.sstep =     {  0, s0, s1, s0 }
  ip.fieldstep = {  0, b0, b1, b0 }
  ip.forcestep = { a0, a1, a1, a0 }
  return ip
end

function getintpat(tau, nsteps, p)
  local eps = tau/nsteps
  local ip = intpat[p.type](eps, p)
  local nfs = ip.nforcesteps
  local nreps = nsteps/nfs
  if nreps ~= math.floor(nreps) then
    printf("getintpat(%s): nsteps (%i) not multiple of %i\n", p.type, nsteps, nfs)
    exit(1)
  end
  --myprint("ip: ", ip, "\n")
  return ip, nreps, eps
end

local function mergeQ(q, qy)
  local qx = {}
  for i=1,#q do qx[i] = q[i] end
  local i,ix,iy = 1,1,1
  local x = qx[ix]
  local y = qy[iy]
  while true do
    if y==nil then
      if x==nil then break end
      q[i] = x
      ix = ix + 1
      x = qx[ix]
    else
      if x==nil or y.s<x.s then
        q[i] = y
        iy = iy + 1
        y = qy[iy]
      else
        q[i] = x
        ix = ix + 1
        x = qx[ix]
      end
    end
    i = i + 1
  end
end

local function addQ(q, ifield, iforce, pat, nreps, tau)
  local q2 = {}
  for i=0,nreps-1 do
    local s = (tau*i)/nreps
    local t = s
    for j=1,pat.nsteps do
      s = s + pat.sstep[j]
      t = t + pat.fieldstep[j]
      if pat.forcestep[j] ~= 0 then
        local z = {}
        z.s = s
        z.t = t
        z.ifield = ifield
        z.iforce = iforce
        z.eps = pat.forcestep[j]
        z.fgeps = (pat.fgstep and pat.fgstep[j]) or 0
        q2[#q2+1] = z
      end
    end
  end
  mergeQ(q, q2)
end

function setupint(f, p)
  local ip = {}
  local nf = f:nfields()
  local nforces = f:nforces()
  local tau = p.tau

  ip.fields = {}
  ip.tau = tau
  ip.q = {}
  for i=1,nf do
    ip.fields[i] = {}
    local ipf = {}
    ip.fields[i].forces = ipf
    for j=1,nforces[i] do
      local fpij = p.forceparams[i][j]
      local nsteps = fpij.nsteps
      local intalg = fpij.intalg or {type="omelyan"}
      local pat,nreps,eps = getintpat(tau, nsteps, intalg)
      ipf[j] = {pat=pat,nreps=nreps,eps=eps}
      addQ(ip.q, i, j, pat, nreps, tau)
    end
  end
  --myprint("q: ", ip.q, "\n")
  return ip
end

function integrate(f, p)
  local seps = 1e-12
  local teps = 1e-12
  local pf = p.fields
  local nf = #pf
  local fieldtime = {}
  local q = p.q
  for i=1,nf do
    fieldtime[i] = 0
  end

  local iq = 1
  while iq<=#q do
    local z = q[iq]
    local s = z.s
    local t = z.t
    local ifield = z.ifield
    local ff,fe,fg = {z.iforce},{z.eps},{z.fgeps}
    local iq1 = iq + 1
    while q[iq1] do
      q1 = q[iq1]
      if separateForces then break end
      if math.abs(s-q1.s)>seps then break end
      if q1.ifield ~= ifield then break end
      if math.abs(t-q1.t)>teps then
        printf("inteagrate: t difference too large %g %g\n", t, q1.t)
        exit(1)
      end
      local k = 1
      while ff[k] and ff[k] ~= q1.iforce do k=k+1 end
      if ff[k] then
        local fek = fe[k]
        fe[k] = fe[k] + q1.eps
        fg[k] = (fek*fg[k] + q1.eps*q1.fgeps)/fe[k]
      else
        ff[k] = q1.iforce
        fe[k] = q1.eps
        fg[k] = q1.fgeps
      end
      iq = iq + 1
      iq1 = iq + 1
    end
    iq = iq + 1

    local d = t - fieldtime[ifield]
    if math.abs(d)>teps then
      --printf("updateField[%i]: %g\t%g\n", ifield, t, d)
      f:updateField(ifield, d)
      fieldtime[ifield] = t
    end
    f:updateMomentum(ifield, ff, fe, fg)
  end
  for i=1,nf do
    local d = p.tau - fieldtime[i]
    if math.abs(d)>teps then
      --printf("updateField[%i]: %g\t%g\n", ifield, t, d)
      f:updateField(i, d)
    end
  end
end

function integrateOld(f, p)
  local pf = p.fields
  local nf = #pf
  local fieldtime = {}
  local q = {}
  local teps = math.huge
  for i=1,nf do
    fieldtime[i] = 0
    local forces = pf[i].forces
    for j=1,#forces do
      local pff = forces[j]
      q[#q+1] = {t=0, i=i, j=j, eps=0, rep=1, n=0}
      if pff.eps < teps then teps = pff.eps end
    end
  end
  teps = 1e-6*teps

  local fi, fj, fe = 0
  local tfinal = 0
  while true do
    local x = q[1]
    if not x then break end
    local i = x.i
    local j = x.j
    local pff = pf[i].forces[j]
    local pat = pff.pat
    --printf("q %i %i t  %-8g  eps  %-8g\n", i, j, x.t, x.eps)
    if x.eps ~= 0 then
      local s = x.t - fieldtime[i]
      if fi==i and s<teps then
	local k = 1
	while(true) do
	  if fj[k]==j then
            local fek = fe[k]
	    fe[k] = fe[k] + x.eps
	    fg[k] = (fek*fg[k] + x.eps*x.fge)/fe[k]
	    break
	  end
	  k = k + 1
	  if k>#fj then
	    fj[k] = j
	    fe[k] = x.eps
	    fg[k] = x.fge
	    break
	  end
	end
      else
	if fi>0 then
	  --myprint("fi: ",fi,"  fj: ",fj,"  fe: ",fe,"\n")
	  f:updateMomentum(fi, fj, fe, fg)
	end
	if s > teps then
	  f:updateField(i, s)
	  fieldtime[i] = x.t
	end
	fi = i
	fj = { j }
	fe = { x.eps }
	fg = { x.fge }
      end
    end
    local n = x.n + 1
    if n > pat.nsteps then
      n = 1
      x.rep = x.rep + 1
    end
    if x.rep <= pff.nreps then
      x.t = x.t + pat.fieldstep[n]
      x.eps = pat.forcestep[n]
      x.fge = (pat.fgstep and pat.fgstep[n]) or 0
      x.n = n
      sort(q, 1)
    else
      if tfinal>0 and math.abs(tfinal-x.t)>teps then
	printf("warning: different tfinals %g %g\n", tfinal, x.t)
      end
      tfinal = x.t
      remove(q, 1)
    end
  end
  if fi>0 then
    f:updateMomentum(fi, fj, fe, fg)
  end
  for i=1,nf do
    local s = tfinal - fieldtime[i]
    if s > teps then
      f:updateField(i, s)
    end
  end
end

--[[
function integrate.4fxx(f, p)
  local u0 = 4*eps*0.1786
  local u1 = 4*eps*-0.0663
  --local h0 = 4*eps*0.7123
  local h0 = 4*eps*0.68
  for ups = 1,nsteps/4 do
    U = act.g:updateU(U, H, u0)
    H = act.g:updateH(H, U, h0)
    H = act.f:updateH(H, U, phi, h0)
    U = act.g:updateU(U, H, u1)
    H = act.g:updateH(H, U, 2*eps-h0)
    H = act.f:updateH(H, U, phi, 2*eps-h0)
    U = act.g:updateU(U, H, 4*eps-2*u0-2*u1)
    H = act.f:updateH(H, U, phi, 2*eps-h0)
    H = act.g:updateH(H, U, 2*eps-h0)
    U = act.g:updateU(U, H, u1)
    H = act.f:updateH(H, U, phi, h0)
    H = act.g:updateH(H, U, h0)
    U = act.g:updateU(U, H, u0)
  end
end
--]]
