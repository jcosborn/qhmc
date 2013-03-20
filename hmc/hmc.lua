local setupint, integrate

function hmcstep(fields, params)
  fields:save()
  local Sold = fields:action()

  local intparams = setupint(fields, params)
  integrate(fields, intparams)

  local Snew = fields:action()

  local ds = Snew - Sold
  printf("Sold: %-12.10g  Snew: %-12.10g\n", Sold, Snew)

  if params.checkReverse then
    fields:reverse()
    integrate(fields, intparams)
    local Srev = fields:action()
    printf("Sold: %-12.10g  Srev: %-12.10g  dS: %-12.10g\n", Sold, Srev, Srev-Sold)
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

local intpat = {}
function intpat.leapfrog(eps, p)
  local s0 = 0.5*eps
  local ip = {}
  ip.nsteps = 2
  ip.nforcesteps = 1
  ip.fieldstep = {  s0, s0 }
  ip.forcestep = { eps,  0 }
  return ip
end
function intpat.omelyan(eps, p)
  local lambda = p.lambda or 0.1932
  local s0 = 2*eps*lambda
  local s1 = 2*(eps-s0)
  local ip = {}
  ip.nsteps = 3
  ip.nforcesteps = 2
  ip.fieldstep = {  s0,  s1, s0 }
  ip.forcestep = { eps, eps,  0 }
  return ip
end
intpat["2MNV"] = function(eps, p)
  local lambda = p.lambda or 0.1932
  local s0 = 2*eps*lambda
  local s1 = 2*(eps-s0)
  local ip = {}
  ip.nsteps = 3
  ip.nforcesteps = 2
  ip.fieldstep = {  0, eps, eps }
  ip.forcestep = { s0,  s1,  s0 }
  return ip
end

function getintpat(tau, nsteps, p)
  local eps = tau/nsteps
  local ip = intpat[p.type](eps, p)
  local nfs = ip.nforcesteps
  local nreps = nsteps/nfs
  if nreps ~= math.floor(nreps) then
    printf("intpat.omelyan: nsteps (%i) not multiple of %i\n", nsteps, nfs)
    exit(1)
  end
  return ip, nreps, eps
end

function setupint(f, p)
  local ip = {}
  local nf = f:nfields()
  local nforces = f:nforces()
  local tau = p.tau

  ip.fields = {}
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
    end
  end
  return ip
end

local function remove(q, k)
  for i=k,#q do
    q[i] = q[i+1]
  end
end

local function sort(q, k)
  while q[k] and q[k+1] and q[k].t > q[k+1].t do
    q[k],q[k+1] = q[k+1],q[k]
    k = k + 1
  end
  while q[k] and q[k-1] and q[k].t < q[k-1].t do
    q[k],q[k-1] = q[k-1],q[k]
    k = k - 1
  end
end

function integrate(f, p)
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
	    fe[k] = fe[k] + x.eps
	    break
	  end
	  k = k + 1
	  if k>#fj then
	    fj[k] = j
	    fe[k] = x.eps
	    break
	  end
	end
      else
	if fi>0 then
	  --myprint("fi: ",fi,"  fj: ",fj,"  fe: ",fe,"\n")
	  f:updateMomentum(fi, fj, fe)
	end
	if s > teps then
	  f:updateField(i, s)
	  fieldtime[i] = x.t
	end
	fi = i
	fj = { j }
	fe = { x.eps }
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
    f:updateMomentum(fi, fj, fe)
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
