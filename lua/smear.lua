local setupSmear = {}
--local setupForce = {}
local dosmear = {}
local dochain = {}

function smearGauge(g, params)
  local sg = g.g
  if params then
    for i=1,#params do
      if not params[i].sg then
        params[i].sg = qopqdp.gauge()
      end
      dosmear[params[i].type](params[i].sg, sg, params[i])
      sg = params[i].sg
    end
  end
  return sg
end

function smearForce(f, g, params)
  if params then
    for i=#params,1,-1 do
      local g0
      if i>1 then g0 = params[i-1].sg
      else g0 = g.g end
      dochain[params[i].type](f, params[i].sg, g0, params[i])
    end
  end
  --f:leftMultAdj(g.g)
  --f:makeAntiHerm()
  f:derivForce(g.g)
end

local function execSmear(sg, g, t)
  if not t.g then
    t.g = {sg,g}
    for i=3,t.ng do t.g[i] = qopqdp.gauge() end
  else
    t.g[1] = sg
    t.g[2] = g
  end
  --myprint("t = ",t,"\n")
  for k,v in ipairs(t.steps) do
    local inn,out = {},{}
    for i=1,#v.inn do inn[i] = t.g[v.inn[i]] end
    for i=1,#v.out do out[i] = t.g[v.out[i]] end
    --myprint("steps[",k,"] = ",v,"\n")
    qopqdp.smear(out, inn, v)
  end
end

local function execForce(f, sg, g, t)
  t.f[1] = f
  t.g[1] = sg
  t.g[2] = g
  --myprint("t = ",t,"\n")
  for k=#t.steps,1,-1 do
    v = t.steps[k]
    local inn,out,fin,fout = {},{},{},{}
    for i=1,#v.inn do inn[i] = t.g[v.inn[i]] end
    for i=1,#v.out do out[i] = t.g[v.out[i]] end
    for i=1,#v.fin do fin[i] = t.f[v.fin[i]] end
    for i=1,#v.fout do
      fout[i] = t.f[v.fout[i]]
      if v.zeroFout[i] then fout[i]:zero() end
    end
    --myprint("steps[",k,"] = ",v,"\n")
    qopqdp.smearChain(fout, fin, out, inn, v)
  end
  local lastf = t.steps[1].fout[1]
  if lastf ~= 1 then
    f:set(t.f[lastf])
  end
  --os.exit()
end

local function setupForce(f, sg, g, t)
  local g2f,fUsed,fmax = {},{},0
  g2f[1] = 1
  fUsed[1] = true
  for k=#t.steps,1,-1 do
    v = t.steps[k]
    v.fin,v.fout,v.zeroFout = {},{},{}
    for i,n in ipairs(v.inn) do
      if g2f[n]==nil then
	local j = 1
	if n~=2 or g2f[2]~=nil or fUsed[j] then
	  repeat j=j+1 until not fUsed[j]
	end
	if j>fmax then fmax = j end
	g2f[n] = j
	fUsed[j] = true
	v.zeroFout[i] = true
      end
      v.fout[i] = g2f[n]
    end
    for i,n in ipairs(v.out) do
      if g2f[n]~=nil then
	-- assume no longer used
	fUsed[g2f[n]] = false
      else
	printf("error: unused force output %i %i %i\n", k, i, n)
      end
      v.fin[i] = g2f[n]
    end
  end
  t.f = { f }
  for i=1,fmax do t.f[i+1] = qopqdp.force() end
end


function setupSmear.staples(t, isg, ig, p)
  local s = {}
  t.steps[#t.steps+1] = s
  s.inn = {ig}
  s.out = {isg}
  s.type = "staples"
  s.topdir = {}
  s.sidedir = {}
  s.toplinknum = {}
  s.sidelinknum = {}
  s.coef = {}
  for mu=1,4 do
    local ns = 0
    s.topdir[mu] = {}
    s.sidedir[mu] = {}
    s.toplinknum[mu] = {}
    s.sidelinknum[mu] = {}
    s.coef[mu] = {}
    if p.coeffs[mu][0] and p.coeffs[mu][0] ~= 0 then
      ns = ns + 1
      s.topdir[mu][ns] = mu
      s.sidedir[mu][ns] = 0
      s.toplinknum[mu][ns] = mu-1
      s.sidelinknum[mu][ns] = -1
      s.coef[mu][ns] = p.coeffs[mu][0]
    end
    for nu=1,4 do
      if p.coeffs[mu][nu] and p.coeffs[mu][nu] ~= 0 then
	ns = ns + 1
	s.topdir[mu][ns] = mu
	s.sidedir[mu][ns] = nu
	s.toplinknum[mu][ns] = mu-1
	s.sidelinknum[mu][ns] = nu-1
	s.coef[mu][ns] = p.coeffs[mu][nu]
      end
      if p.coeffs[mu][-nu] and p.coeffs[mu][-nu] ~= 0 then
	ns = ns + 1
	s.topdir[mu][ns] = mu
	s.sidedir[mu][ns] = -nu
	s.toplinknum[mu][ns] = mu-1
	s.sidelinknum[mu][ns] = nu-1
	s.coef[mu][ns] = p.coeffs[mu][-nu]
      end
    end
  end
end







function dosmear.staples(sg, g, p)
  if not p.steps then
    p.steps = {}
    p.ng = 2
    setupSmear.staples(p, 1, 2, p)
  end
  execSmear(sg, g, p)
end

function dochain.staples(f, sg, g, p)
  if not p.f then 
    setupForce(f, sg, g, p)
  end
  execForce(f, sg, g, p)
end

function dosmear.fat7(sg, g, p)
  --sg:set(g)
  qopqdp.smear({sg}, {g}, p)
end

function dochain.fat7(f, sg, g, p)
  if not p.f then p.f = qopqdp.force() end
  --p.f:set(f)
  p.f:zero()
  qopqdp.smearChain({p.f}, {f}, {sg}, {g}, p)
  f:set(p.f)
end

function dosmear.projectU(sg, g, p)
  qopqdp.smear({sg}, {g}, p)
end

function dochain.projectU(f, sg, g, p)
  if not p.f then p.f = qopqdp.force() end
  p.f:zero()
  qopqdp.smearChain({p.f}, {f}, {sg}, {g}, p)
  f:set(p.f)
end

function dosmear.projectRat(sg, g, p)
  if not p.smearg then
    p.smearg = qopqdp.gauge()
  end
  if not p.plaq then
    p.plaq = {type="product",adj={false,true}}
    p.plaqg = qopqdp.gauge()
  end
  if not p.ah then
    p.ah = {type="antiherm"}
    p.ahg = qopqdp.gauge()
  end
  if not p.rat then
    p.rat = {type="mobius",coeffs={1,0.5*p.rho,1,-0.5*p.rho}}
    p.ratg = qopqdp.gauge()
  end
  if not p.proj then
    p.proj = {type="product",adj={false,false}}
  end
  qopqdp.smear({p.smearg}, {g}, p.smear)
  qopqdp.smear({p.plaqg}, {p.smearg, g}, p.plaq)
  qopqdp.smear({p.ahg}, {p.plaqg}, p.ah)
  qopqdp.smear({p.ratg}, {p.ahg}, p.rat)
  qopqdp.smear({sg}, {p.ratg, g}, p.proj)
end

function dochain.projectRat(f, sg, g, p)
  if not p.f then
    p.f = qopqdp.force()
    p.fc = qopqdp.force()
  end
  p.f:zero()
  p.fc:zero()
  qopqdp.smearChain({p.fc, p.f}, {f}, {sg}, {p.ratg, g}, p.proj)
  f:set(p.f)
  p.f:zero()
  qopqdp.smearChain({p.f}, {p.fc}, {p.ratg}, {p.ahg}, p.rat)
  p.fc:zero()
  qopqdp.smearChain({p.fc}, {p.f}, {p.ahg}, {p.plaqg}, p.ah)
  p.f:zero()
  qopqdp.smearChain({p.f, f}, {p.fc}, {p.plaqg}, {p.smearg, g}, p.plaq)
  qopqdp.smearChain({f}, {p.f}, {p.smearg}, {g}, p.smear)
end


-- Cayley-Hamilton projection

local function set_projectCH(t, isg, ig, ig0, p)
  local n = #t.steps + 1
  local ng = t.ng + 1
  t.steps[n]   = {type="product",adj={false,true},inn={ig,ig0},out={ng}}
  t.steps[n+1] = {type="mobius",coeffs={1,0,1,1},inn={ng},out={ng+1}}
  t.steps[n+2] = {type="antiherm",inn={ng+1},out={ng+2}}
  t.steps[n+3] = {type="mobius",coeffs={1,-2,1,2},inn={ng+2},out={ng+3}}
  t.steps[n+4] = {type="product",adj={false,false},inn={ng+3,ig0},out={isg}}
  t.ng = ng + 3
end

function dosmear.projectCH(sg, g, p)
  if not p.steps then
    p.steps = {}
    p.ng = 2
    p.g = { sg, g }
    set_projectCH(p, 1, 2)
  end
end

--sg: smeared gauge field
--g : original gauge field
--p : smearing parameters
function dosmear.stout(sg, g, p)
--if fat7 is not set
  if not p.fat7 then
--just include the staples in the fat7 construction
    p.fat7 = {type="fat7",coeffs={three_staple=1}}
    p.fat7g = qopqdp.gauge()
  end
--if plaq is not set
  if not p.plaq then
-- U U^+
    p.plaq = {type="product",adj={false,true}}
    p.plaqg = qopqdp.gauge()
  end
-- if ah is not set
  if not p.ah then
--P_{taH}
    p.ah = {type="tracelessAntiherm"}
    --p.ah = {type="mobius",coeffs={0,1,1,0}}
    p.ahg = qopqdp.gauge()
  end
-- if exp is not set
  if not p.exp then
-- exponentiate
    p.exp = {type="exp",rho={p.rho,p.rho,p.rho,p.rho}}
    --p.exp = {type="mobius",coeffs={1,0.5*p.rho,1,-0.5*p.rho}}
    --p.exp = {type="mobius",coeffs={0.1,p.rho,1,0}}
    p.expg = qopqdp.gauge()
  end
-- if stout is not set
  if not p.stout then
-- U U
    p.stout = {type="product",adj={false,false}}
  end
-- do the stout smearing. follows the convention in Morningstar&Peardon 2004.
-- C_\mu(x)
  qopqdp.smear({p.fat7g}, {g}, p.fat7) 
-- \Omega = C_\mu(x) U_\mu^+(x) 
  qopqdp.smear({p.plaqg}, {p.fat7g, g}, p.plaq) 
-- i Q_\mu(x) 
  qopqdp.smear({p.ahg}, {p.plaqg}, p.ah) 
-- exp(i\rho_\mu Q_\mu)
  qopqdp.smear({p.expg}, {p.ahg}, p.exp) 
-- exp(i\rho Q) U
  qopqdp.smear({sg}, {p.expg, g}, p.stout) 
end

function dochain.stout(f, sg, g, p)
  if not p.f then
    p.f = qopqdp.force()
    p.fc = qopqdp.force()
  end
  p.f:zero()
  p.fc:zero() 
  -- same as "product": exp(iQ) U
  qopqdp.smearChain({p.fc, p.f}, {f}, {sg}, {p.expg, g}, p.stout)
  -- now p.fc = f sg, p.f = expg f
  f:set(p.f) -- f = p.f = expg f
  p.f:zero()
  qopqdp.smearChain({p.f}, {p.fc}, {p.expg}, {p.ahg}, p.exp)
  -- now p.f = the exponential derivative 
  p.fc:zero()
  qopqdp.smearChain({p.fc}, {p.f}, {p.ahg}, {p.plaqg}, p.ah)
  p.f:zero()
  qopqdp.smearChain({p.f, f}, {p.fc}, {p.plaqg}, {p.fat7g, g}, p.plaq)
  qopqdp.smearChain({f}, {p.f}, {p.fat7g}, {g}, p.fat7)
end

local function appendUniq(t, a)
  for i=1,#t do
    local same = true
    for j=1,#a do
      if t[i][j] ~= a[j] then same=false; break end
    end
    if same then return i end
  end
  t[#t+1] = a
  return #t
end

function hypStaples(ad, alpha)
  local nd = #ad
  local s = { type="staples" }
  s.topdir = {}
  s.sidedir = {}
  s.toplinknum = {}
  s.sidelinknum = {}
  s.coef = {}
  local ad2 = {}
  for mu=1,nd do ad2[mu] = {} end
  for mu=1,nd do
    for i=1,#ad[mu] do
      local io = (i-1)*nd + mu
      local a = ad[mu][i]
      local na,nl,n = 0,0,1
      for nu=1,nd do nl=nl+a[nu]; if nu~=mu and a[nu]>0 then na=na+1 end end
      local ap = alpha/(2*na)
      s.topdir[io] = { [n] = mu }
      s.sidedir[io] = { [n] = 0 }
      s.toplinknum[io] = { [n] = mu-1 }
      s.sidelinknum[io] = { [n] = -1 }
      s.coef[io] = { [n] = 1-alpha }
      for nu=1,nd do
	if nu~=mu and a[nu]>0 then
	  a[nu] = a[nu] - 1
	  local tn = mu - 1
	  local sn = nu - 1
	  if nl>1 then
	    local ti = appendUniq(ad2[mu], {table.unpack(a)})
	    local si = appendUniq(ad2[nu], {table.unpack(a)})
	    tn = tn + ti*nd
	    sn = sn + si*nd
	  end
	  n = n + 1
	  s.topdir[io][n] = mu
	  s.sidedir[io][n] = nu
	  s.toplinknum[io][n] = tn
	  s.sidelinknum[io][n] = sn
	  s.coef[io][n] = ap
	  n = n + 1
	  s.topdir[io][n] = mu
	  s.sidedir[io][n] = -nu
	  s.toplinknum[io][n] = tn
	  s.sidelinknum[io][n] = sn
	  s.coef[io][n] = ap
	  a[nu] = a[nu] + 1
	end
      end
    end
  end
  return s, ad2
end

function setupSmear.hyp(t, isg, ig, p)
  local allowedDirs = {
    [1]={{0,1,1,1}}, [2]={{1,0,1,1}}, [3]={{1,1,0,1}}, [4]={{1,1,1,0}}
  }
  --[[ local allowedDirs = {
  [1]={{0,1,0,0},{0,0,1,0},{0,0,0,1}},
  [2]={{1,0,0,0},{0,0,1,0},{0,0,0,1}},
  [3]={{1,0,0,0},{0,1,0,0},{0,0,0,1}},
  [4]={{1,0,0,0},{0,1,0,0},{0,0,1,0}},
  --}--]]
  --[[local allowedDirs = {
  [1]={{0,1,1,0}},
  [2]={{1,0,1,0}},
  [3]={{1,1,0,0}},
  [4]={{1,1,0,0}},
  --}--]]
  local nlevels = 3
  local rsteps = {}
  local ng = t.ng
  local pout = { isg }
  for i=nlevels,1,-1 do
    local sout = {}
    for j=1,#allowedDirs[1] do
      ng = ng + 1
      sout[j] = ng
      if not p.projectCH then
	rsteps[#rsteps+1] = {type="projectU",inn={ng},out={pout[j]}}
      else
	local tt = { ng=ng, steps={} }
	set_projectCH(tt, pout[j], ng, ig, nil)
	ng = tt.ng
	for k=#tt.steps,1,-1 do rsteps[#rsteps+1]=tt.steps[k] end
      end
      --sout[j] = pout[j]
    end
    local st
    st,allowedDirs = hypStaples(allowedDirs, p.alpha[i])
    st.out = sout
    st.inn = { ig }
    pout = {}
    for j=1,#allowedDirs[1] do
      ng = ng + 1
      st.inn[j+1] = ng
      pout[j] = ng
    end
    rsteps[#rsteps+1] = st
  end
  t.ng = ng
  for i=#rsteps,1,-1 do
    t.steps[#t.steps+1] = rsteps[i]
  end
end

function dosmear.hyp(sg, g, p)
  if not p.steps then
    p.steps = {}
    p.ng = 2
    setupSmear.hyp(p, 1, 2, p)
  end
  execSmear(sg, g, p)
end

function dochain.hyp(f, sg, g, p)
  if not p.f then 
    setupForce(f, sg, g, p)
  end
  execForce(f, sg, g, p)
end


function tableSetIJ(t,k,i,j,v)
  if t[k] == nil then t[k] = {} end
  if t[k][i] == nil then t[k][i] = {} end
  t[k][i][j] = v
end
function addStaple(s,i,j,td,sd,tn,sn,c)
  tableSetIJ(s,"topdir",i,j,td)
  tableSetIJ(s,"sidedir",i,j,sd)
  tableSetIJ(s,"toplinknum",i,j,tn)
  tableSetIJ(s,"sidelinknum",i,j,sn)
  tableSetIJ(s,"coef",i,j,c)
end

function dosmear.hyp2(sg, g, p)
  if not p.staples then
    p.staples = {}
    p.inn = {}
    p.out = {}

    p.inn[1] = { g }
    --p.out[1] = {}
    --p.out[1][1] = qopqdp.gauge()
    --p.out[1][2] = qopqdp.gauge()
    local s = { type="staples" }
    addStaple(s,1,1, 1,0, 0,-1, 1)
    --addStaple(s,5,1, 1,2, 0, 1, 1)
    addStaple(s,2,1, 2,0, 1,-1, 1)
    --addStaple(s,6,1, 2,1, 1, 0, 1)
    addStaple(s,3,1, 3,0, 2,-1, 1)
    --addStaple(s,7,1, 3,1, 2, 0, 1)
    addStaple(s,4,1, 4,0, 3,-1, 1)
    --addStaple(s,8,1, 4,1, 3, 0, 1)
    p.staples[1] = s

    --p.inn[2] = { p.out[1][1] }
    --p.inn[2] = { p.out[1][1], p.out[1][2] }
    p.out[1] = { sg }
    s = { type="staples" }
    local f = 0.99
    addStaple(s,1,1, 1,2, 0, 1, f)
    --addStaple(s,1,2, 1,2, 0, 1, f)
    addStaple(s,2,1, 2,1, 1, 0, f)
    --addStaple(s,2,2, 2,1, 1, 0, f)
    addStaple(s,3,1, 3,0, 2,-1, f)
    --addStaple(s,3,2, 3,1, 2, 0, f)
    addStaple(s,4,1, 4,0, 3,-1, f)
    --addStaple(s,4,2, 4,1, 3, 0, f)
    p.staples[1] = s
  end

  p.inn[1][1] = g  -- might have changed since last call
  p.out[1][1] = sg -- might have changed since last call
  qopqdp.smear(p.out[1], p.inn[1], p.staples[1])
  --qopqdp.smear(p.out[2], p.inn[2], p.staples[2])
end

function dochain.hyp2(f, sg, g, p)
  if not p.f then
    p.f = {}
    p.fc = {}

    p.f[1] = { f }
    p.fc[1] = {}
    p.fc[1][1] = qopqdp.force()
    --p.fc[1][2] = qopqdp.force()

    --p.f[2] = { p.fc[1][1] }
    --p.f[2] = { p.fc[1][1], p.fc[1][2] }
    --p.fc[2] = { f }
  end

  p.inn[1][1] = g  -- might have changed since last call
  p.out[1][1] = sg -- might have changed since last call
  p.f[1][1] = f  -- might have changed since last call
  --p.fc[2][1] = f -- might have changed since last call
  --p.f[2][1]:zero()
  --p.f[2][2]:zero()
  --qopqdp.smearChain(p.f[2], p.fc[2], p.out[2], p.inn[2], p.staples[2])
  p.fc[1][1]:set(f)
  f:zero()
  qopqdp.smearChain(p.f[1], p.fc[1], p.out[1], p.inn[1], p.staples[1])
end
