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

function dosmear.stout(sg, g, p)
  if not p.fat7 then
    p.fat7 = {type="fat7",coeffs={three_staple=1}}
    p.fat7g = qopqdp.gauge()
  end
  if not p.plaq then
    p.plaq = {type="product",adj={false,true}}
    p.plaqg = qopqdp.gauge()
  end
  if not p.ah then
    p.ah = {type="antiherm"}
    --p.ah = {type="mobius",coeffs={0,1,1,0}}
    p.ahg = qopqdp.gauge()
  end
  if not p.exp then
    --p.exp = {type="exp",rho={p.rho,p.rho,p.rho,p.rho}}
    p.exp = {type="mobius",coeffs={1,0.5*p.rho,1,-0.5*p.rho}}
    --p.exp = {type="mobius",coeffs={0.1,p.rho,1,0}}
    p.expg = qopqdp.gauge()
  end
  if not p.stout then
    p.stout = {type="product",adj={false,false}}
  end
  qopqdp.smear({p.fat7g}, {g}, p.fat7)
  qopqdp.smear({p.plaqg}, {p.fat7g, g}, p.plaq)
  qopqdp.smear({p.ahg}, {p.plaqg}, p.ah)
  qopqdp.smear({p.expg}, {p.ahg}, p.exp)
  qopqdp.smear({sg}, {p.expg, g}, p.stout)
end

function dochain.stout(f, sg, g, p)
  if not p.f then
    p.f = qopqdp.force()
    p.fc = qopqdp.force()
  end
  p.f:zero()
  p.fc:zero()
  qopqdp.smearChain({p.fc, p.f}, {f}, {sg}, {p.expg, g}, p.stout)
  f:set(p.f)
  p.f:zero()
  qopqdp.smearChain({p.f}, {p.fc}, {p.expg}, {p.ahg}, p.exp)
  p.fc:zero()
  qopqdp.smearChain({p.fc}, {p.f}, {p.ahg}, {p.plaqg}, p.ah)
  p.f:zero()
  qopqdp.smearChain({p.f, f}, {p.fc}, {p.plaqg}, {p.fat7g, g}, p.plaq)
  qopqdp.smearChain({f}, {p.f}, {p.fat7g}, {g}, p.fat7)
end
