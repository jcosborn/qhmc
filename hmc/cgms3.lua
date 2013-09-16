function cgms2(x, src, op, shifts, resid, opts, subset)
  local nx = #x
  local imin = 1
  local maxits = opts.restart or 500
  local a,b,p,pAp = {},{},{},{}
  local r,Ap = {},{}
  local Ar,Apo = src:clone(subset),src:clone(subset)
  local rsq = src:norm2(subset)
  local rs,rso = {},{}
  for i=1,nx do
    x[i]:zero(subset)
    p[i] = src:clone(subset)
    r[i] = src:clone(subset)
    Ap[i] = x[i]:clone(subset)
    a[i] = 0
    rs[i] = rsq
  end
  local rsqstop = resid*resid*rsq
  local its = 0
  while true do
    oldrsq = rsq
    Apo,Ap[imin] = Ap[imin],Apo
    op(Ap[imin], p[imin])
    Ap[imin]:combine({Ap[imin],p[imin]}, {1,shifts[imin]}, subset)
    its = its + 1

    Ar:combine({Ap[imin],Apo},{1,-a[imin]},subset)
    for i=1,nx do
      if i~=imin then
	local cr = shifts[i] - shifts[imin]
	Ap[i]:combine({Ap[i],Ar,r[imin]},{a[i],1,cr},subset)
	--op(Ap[i], p[i])
	--Ap[i]:combine({Ap[i],p[i]}, {1,shifts[i]}, subset)
      end
    end

    for i=1,nx do
      if i~=imin then
	pAp[i] = p[i]:Re_dot(Ap[i], subset)
	if pAp[i]<=0 then b[i] = 0
	else b[i] = p[i]:dot(r[i],subset)/pAp[i] end
	--else b[i] = p[i]:Re_dot(r[i],subset)/pAp[i] end
	--else b[i] = r[imin]:Re_dot(r[i],subset)/pAp[i] end
	--printf("%i\t%s\t%g\n", i, tostring(b[i]), p[i]:Re_dot(r[i],subset)/pAp[i])
      end
    end

    pAp[imin] = p[imin]:Re_dot(Ap[imin], subset)
    if pAp[imin]<=0 then break end  -- loss of precision in calculating pAp
    b[imin] = rsq / pAp[imin]

    for i=1,nx do
      x[i]:combine({x[i],p[i]},{1,b[i]}, subset)
      r[i]:combine({r[i],Ap[i]},{1,-b[i]}, subset)
    end
    --r[imin]:combine({r[imin],Ap[imin]},{1,-b[imin]}, subset)
    rsq = r[imin]:norm2(subset)
    --printf("%i\t%g\n", its, rsq)
    if its >= maxits or rsq < rsqstop then break end

    a[imin] = rsq / oldrsq
    for i=1,nx do
      if i~=imin then
	rs[i] = r[i]:norm2(subset)
	--printf("%i\t%g\n", i, rs[i])
	--local t = math.sqrt(rs[i]*rsq)/rso[i]
	a[i] = -Ap[i]:Re_dot(r[imin])/pAp[i]
	--a[i] = -Ap[i]:dot(r[imin])/pAp[i]
      end
    end

    for i=1,nx do
      p[i]:combine({p[i],r[imin]},{a[i],1},subset)
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
      printf("%i\t%i\t%g\n", i, it, rsq)
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
  return info
end
