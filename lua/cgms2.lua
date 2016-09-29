function cgms2(x, src, op, shifts, resid, opts, subset)
  local nx = #x
  local imin = 1
  local maxits = opts.restart or 500
  local a,b,bo,p,z,zn,zo,pAp = {},{},{},{},{},{},{}
  for i=1,nx do
    x[i]:zero(subset)
    p[i] = src:clone(subset)
    a[i] = 0
    bo[i] = -1
    z[i] = 1
    zo[i] = 1
  end
  local r,Ap = src:clone(subset),src:clone(subset)
  local rsq = r:norm2(subset)
  local rsqstop = resid*resid*rsq
  local its = 0
  while true do
    oldrsq = rsq
    op(Ap, p[imin])
    Ap:combine({Ap,p[imin]}, {1,shifts[imin]}, subset)
    its = its + 1

    pAp = p[imin]:reDot(Ap, subset)
    if pAp<=0 then break end  -- loss of precision in calculating pAp
    b[imin] = rsq / pAp
    zn[imin] = 1
    for i=1,nx do
      if i~=imin then
        zn[i] = z[i]*zo[i]*bo[imin]
        local c1 = b[imin]*a[imin]*(zo[i]-z[i])
        c1 = c1 + zo[i]*bo[imin]*(1+shifts[i]*b[imin])
        if c1~=0 then zn[i] = zn[i]/c1
	else zn[i] = 0 end
        if z[i]~=0 then b[i] = b[imin]*zn[i]/z[i]
	else zn[i],b[i] = 0,0 end
      end
    end

    for i=1,nx do
      x[i]:combine({x[i],p[i]},{1,b[i]}, subset)
    end

    r:combine({r,Ap},{1,-b[imin]}, subset)
    rsq = r:norm2(subset)
    --printf("%i\t%g\n", its, rsq)
    if its >= maxits or rsq < rsqstop then break end

    a[imin] = rsq / oldrsq
    for i=1,nx do
      if i~=imin then
        local c2 = z[i]*b[imin]
        if c2~=0 then a[i] = a[imin]*zn[i]*b[i]/c2
	else a[i] = 0 end
      end
    end

    p[imin]:combine({p[imin],r},{a[imin],1},subset)
    for i=1,nx do
      if i~=imin then
	p[i]:combine({p[i],r},{a[i],zn[i]},subset)
      end
      bo[i] = b[i]
      zo[i] = z[i]
      z[i] = zn[i]
    end
  end

  return {its=its, flops=0}
end

function cgms(x, src, op, shifts, resid, opts)
  --local t0 = qopqdp.dtime()
  local info = cgms2(x, src, op, shifts, resid, opts)
  local ssq = src:norm2()
  local rsqstop = resid*resid*ssq
  local r = src:clone()
  local xx = x[1]:clone()
  local its = 0
  for i=1,#x do
    local it = 0
    while(true) do
      op(r, x[i])
      r:combine({r,x[i],src},{-1,-shifts[i],1})
      local rsq = r:norm2()
      --printf("%i\t%i\t%g\n", i, it, rsq)
      --break
      if rsq<=rsqstop then break end
      local res = math.sqrt(rsqstop/rsq)
      local inf = cgms2({xx}, r, op, {shifts[i]}, res, opts)
      x[i]:combine({x[i],xx},{1,1})
      it = it + inf.its
      info.flops = info.flops + inf.flops
    end
    if it>its then its = it end
  end
  --collectgarbage()
  info.its = info.its + its
  info.rsq = 0
  return info
end
