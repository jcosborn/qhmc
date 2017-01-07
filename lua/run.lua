require 'common'
require 'gaugeact'
require 'asqtadact'
require 'hisqact'
require 'wilsonact'
require 'wilson2fact'
require 'fields'
require 'hmc'

local actsmt = {}
actsmt.__index = actsmt
function setupacts(p)
  local acts = {}
  p.vol = 1
  printf("latsize =")
  for k,v in ipairs(p.latsize) do p.vol=p.vol*v; printf(" %i",v); end
  printf("\n")
  printf("seed = %i\n", p.seed)
  printf("beta = %g\n", p.beta)
  printf("u0 = %g\n", p.u0)
  printf("nf = %g\n", p.nf)
  acts.g = gaugeact(p)
  if p.fermact.type == "wilson" then
    acts.f = wilsonact(acts.g, p.fermact)
  elseif p.fermact.type == "wilson2f" then
    acts.f = wilson2fact(acts.g, p.fermact)
  elseif p.fermact.type == "hisq" then
    acts.f = hisqact(acts.g, 1, p.fermact.rhmc)
  else
    acts.f = asqtadact(acts.g, p.fermact)
  end
  acts.fields = setupfields(acts, p)
  return setmetatable(acts, actsmt)
end

function actsmt.unit(a)
  a.fields.G:unit()
end

function actsmt.warm(a, w)
  a.fields.G:warm(w)
end

function actsmt.load(a, fn)
  printf("loading lattice: %s\n", fn)
  local t0 = clock()
  a.fields.G:load(fn)
  t0 = clock() - t0
  printf("load time: %g seconds\n", t0)
  local devavg,devmax = a.fields.G:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  a.fields.G:makeSU()
  devavg,devmax = a.fields.G:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end

function actsmt.save(a, fn, ...)
  printf("saving lattice: %s\n", fn)
  local t0 = clock()
  a.fields.G:save(fn, ...)
  t0 = clock() - t0
  printf("save time: %g seconds\n", t0)
end

local function measure(a, r)
  local t0 = clock()

  local ps,pt = a.fields.G:plaq()
  printf("MEASplaq ss: %-8g  st: %-8g  tot: %-8g\n", ps, pt, 0.5*(ps+pt))

  local pl = a.fields.G:ploop()
  local pls = 0
  local plt = pl[#pl]
  for i=1,#pl-1 do
    pls = pls + pl[i]
  end
  pls = pls/(#pl-1)
  --printf("ploop: %s\n", tostring(plp))
  printf("MEASploop: spatial: %g %g temporal: %g %g\n",pls.r,pls.i,plt.r,plt.i)

  for i,v in ipairs(r.pbp) do
    for j=1,v.reps do
      local cc,nrm2 = a.f:pbp(a.fields.G, v.mass, v.resid, v.opts)
      --printf("pbp nrm2 %g : %g\n", v.mass, nrm2)
      printf("MEASpbp mass %g : %g\n", v.mass, cc)
    end
  end

  --[[
  local v = r.pbp[1]
  local pions = a.f:pions(a.fields.G, v.mass, v.resid, v.opts)
  for k,v in pairs(pions) do
    for i=1,#v do
      printf("%s %i\t%g\n", k, i-1, v[i])
    end
  end
  --]]

  -- user defined measurements
  if r.meas then r.meas(a, r) end

  t0 = clock() - t0
  printf("measurement time: %g seconds\n", t0)
end

function actsmt.run(a, r)
  --[[
  printf("ntraj = %g\n", r.ntraj)
  printf("tau = %g\n", r.tau)
  for i,v in ipairs(r.forceparams[1]) do
    printf("nsteps[%i] = %g\n", i, v.nsteps)
  end
  --]]

  local totaltime = clock()
  --measure(a, r)
  for traj=1,r.ntraj do
    local trajtime = clock()
    a.g:clearStats()
    a.f:clearStats()

    local t0 = clock()
    hmcstep(a.fields, r)
    t0 = clock() - t0
    printf("traj %i secs: %g\n", traj, t0)

    a.g:updateStats()
    a.f:updateStats()
    printf("GF     secs: %8.3f %3.0f%% calls: %4.0f mflops: %5.0f rms: %6.4f max: %6.4f\n", a.g.GFtime, 100*a.g.GFtime/t0, a.g.GFn, a.g.GFmflops, a.g.GFrms, a.g.GFmax)
    printf("GU     secs: %8.3f %3.0f%% calls: %4.0f\n", a.g.GUtime, 100*a.g.GUtime/t0, a.g.GUn)
    printf("LL     secs: %8.3f %3.0f%% calls: %4.0f mflops: %5.0f\n", a.f.LLtime, 100*a.f.LLtime/t0, a.f.LLn, a.f.LLmflops)
    local fft,ffn,fff = 0,0,0
    for i=1,a.f.nff do
      printf("FF[%02i] secs: %8.3f %3.0f%% calls: %4.0f mflops: %5.0f rms: %6.4f max: %6.4f\n", i, a.f.FFtime[i], 100*a.f.FFtime[i]/t0, a.f.FFn[i], a.f.FFmflops[i], a.f.FFrms[i], a.f.FFmax[i])
      fft = fft + a.f.FFtime[i]
      ffn = ffn + a.f.FFn[i]
      fff = fff + a.f.FFflops[i]
    end
    printf("FFtot  secs: %8.3f %3.0f%% calls: %4.0f mflops: %5.0f\n", fft, 100*fft/t0, ffn, 1e-6*fff/fft)
    local cgt,cgn,cgf,cgi,cgm = 0,0,0,0,0
    for i=1,a.f.ncg do
      local ai = 0
      if a.f.CGn[i] > 0 then ai = a.f.CGits[i]/a.f.CGn[i] end
      --printf("CG[%02i] secs: %8.3f %3.0f%% calls: %4.0f mflops: %5.0f avgits: %5.0f max: %5.0f\n", i, a.f.CGtime[i], 100*a.f.CGtime[i]/t0, a.f.CGn[i], a.f.CGmflops[i], ai, a.f.CGmaxits[i])
      printf("CG[%02i] secs: %8.3f %3.0f%% calls: %4.0f mflops: %5.0f avgits: %5.0f max: %5.0f maxres: %g\n", i, a.f.CGtime[i], 100*a.f.CGtime[i]/t0, a.f.CGn[i], a.f.CGmflops[i], ai, a.f.CGmaxits[i], a.f.CGmaxresid[i])
      cgt = cgt + a.f.CGtime[i]
      cgn = cgn + a.f.CGn[i]
      cgf = cgf + a.f.CGflops[i]
      cgi = cgi + a.f.CGits[i]
      cgm = math.max(cgm, a.f.CGmaxits[i])
    end
    printf("CGtot  secs: %8.3f %3.0f%% calls: %4.0f mflops: %5.0f avgits: %5.0f max: %5.0f\n", cgt, 100*cgt/t0, cgn, 1e-6*cgf/cgt, cgi/cgn, cgm)
    local ot = t0 - a.g.GFtime - a.g.GUtime - a.f.LLtime - fft - cgt
    printf("other  secs: %8.3f %3.0f%%\n", ot, 100*ot/t0)

    local keys = {}
    local keylen = 0
    for k in pairs(a.f.FFn) do if type(k) == "string" then keys[#keys+1] = k; if #k>keylen then keylen = #k end end end
    table.sort(keys)
    -- width has to be 2 digits at most (lstrlib.c), so we give up for long keylen.
    local fs = keylen>98 and "%s" or "%-"..tostring(keylen+1).."s"
    for i,k in ipairs(keys) do
      printf(fs.." %-12g %-12g\n", k, a.f.FFrms[k], a.f.FFmax[k])
    end
    io.stdout:flush()

    measure(a, r)
    printf("trajectory time: %g seconds\n", clock()-trajtime)
    io.stdout:flush()
  end
  totaltime = clock() - totaltime
  printf("total time: %g seconds\n", totaltime)
end
