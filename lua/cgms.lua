function cg(x, src, op, resid, opts)
  --local t0 = qopqdp.dtime()
  x:zero()
  local r,p,Ap = src:clone(),x:clone(),x:clone()
  local rsq,rsqold = r:norm2(),0
  local rsqstop = resid*resid*rsq
  local its = 0
  while true do
    p:combine({p,r},{1,1/rsq})
    op(Ap, p)
    its = its + 1
    local pAp = p:Re_dot(Ap)
    local alpha = 1/pAp
    x:combine({x,p},{1,alpha})
    r:combine({r,Ap},{1,-1*alpha})
    rsqold = rsq
    rsq = r:norm2()
    if rsq<rsqstop then break end
  end
  --local t1 = qopqdp.dtime()
  --printf("its: %i\tsecs: %g\trsq: %g\n", its, t1-t0, rsq)
  return {its=its, flops=0, rsq=0}
end

function cgwrap(x, src, op, resid, opts)
  local r = src:clone()
  local xx = x:clone()
  local ssq = src:norm2()
  local rsqstop = resid*resid*ssq
  local info = {its=0,flops=0}
  while(true) do
    op(r, x)
    r:combine({r,src},{-1,1})
    local rsq = r:norm2()
    if rsq<=rsqstop then break end
    local res = math.sqrt(rsqstop/rsq)
    local inf = cg(xx, r, op, res, opts)
    x:combine({x,xx},{1,1})
    info.its = info.its + inf.its
    info.flops = info.flops + inf.flops
  end
  info.rsq = 0
  return info
end

function cgms(x, src, op, shifts, resid, opts)
  --local t0 = qopqdp.dtime()
  local its,flops = 0,0
  local nx = #x
  for i=nx,1,-1 do
    local function op2(d,s)
      op(d, s)
      d:combine({d,s},{1,shifts[i]})
    end
    if i==nx then x[i]:zero()
    else x[i]:set(x[i+1]) end
    local info = cgwrap(x[i], src, op2, resid, opts)
    its = math.max(its,info.its)
    flops = flops + info.flops
  end
  return {its=its, flops=flops, rsq=0}
end
