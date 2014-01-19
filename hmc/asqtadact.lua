require 'smear'

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
	t.resid = 1e-5 or t.resid
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
      fa.resid = 1e-5 or t.resid
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
    t.resid = 1e-5 or t.resid
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
	  t.qt:random("all")
	  a.h:Ddag(t.qt2, t.qt, t.allmass2, "all", "all")
        else
	  t.qt2:random("all")
        end
        a.h:D(t.qt, t.qt2, t.allmass, "even", "all")
      else
	t.qt:random("even")
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
    local fa = r.FA
    if fa.pt then
      a:solve(fa.pt, a.pseudo[i], fa.masses, fa.resid, "even", fa.solveopts, fa.cgnum)
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

  f.f:update(tmin.ff.f, 1)
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

-- Use 0 for no phase, 1 for phase.
function make_phase_term(xsign,ysign,zsign,tsign)
  return xsign + 2*ysign + 4*zsign + 8*tsign
end

function actmt.pions(a, g, wallslice, mass, resid, opts)
  a:set(g, 2)
  local x,y,z,t,p5,ri,ps = {},{},{},{},{},{},{}
  local x2,y2 = {},{}
  local x3,y3 = {},{}
  local ps = {}
  for srcnum=1,#wallslice do

    printf("Start point source %i at t=%i.\n", srcnum, wallslice[srcnum]);
    for i=1,3 do
      printf("Start Color %i.\n", i)
      io.stdout:flush()
      x[i] = getqt(a, 3*i)
      y[i] = getqt(a, 3*i+1)
      z[i] = getqt(a, 3*i+2) -- For the rho sign change.
      x[i]:zero()
      x[i]:point({0,0,0,wallslice[srcnum]},i,1,0)
      a:solve({y[i]}, x[i], {mass}, resid, "all", opts, 0)
      do -- Pion 0_A^(-+)
	t[i] = y[i]:norm2("timeslices")
	for j=1,#t[i] do
	  j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 --Compensate for shifted wall source
	  p5[j] = t[i][j_real]/(#wallslice) + (p5[j] or 0)
	end
      end
      do -- Rho 1_A^(--)
	z[i]:set(y[i]) -- Make a copy.
	-- First do the x. Our source was at 0, so the relative term is 0.
	z[i]:rephase(make_phase_term(1,0,0,0), {0,0,0,wallslice[srcnum]})
	t[i] = y[i]:Re_dot(z[i], "timeslices");
	for j = 1,#t[i] do
	  j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 --Compensate for shifted wall source
	  ri[j] = t[i][j_real]/(#wallslice) + (ri[j] or 0)
	end
	-- Next do the y. The phase undoes the x phase and applies the y phase.
	z[i]:rephase(make_phase_term(1,1,0,0), {0,0,0,wallslice[srcnum]})
	t[i] = y[i]:Re_dot(z[i], "timeslices");
	for j = 1,#t[i] do
	  j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 --Compensate for shifted wall source
	  ri[j] = t[i][j_real]/(#wallslice) + (ri[j] or 0)
	end
	-- Next do the z.
	z[i]:rephase(make_phase_term(0,1,1,0), {0,0,0,wallslice[srcnum]})
	t[i] = y[i]:Re_dot(z[i], "timeslices");
	for j = 1,#t[i] do
	  j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 --Compensate for shifted wall source
	  ri[j] = t[i][j_real]/(#wallslice) + (ri[j] or 0)
	end
      end
      printf("End Color.\n");
      io.stdout:flush();
    end

    -- The pion scalar state. (0_S^(-+))
    for i=1,3 do
      printf("Start color %i.\n", i);
      x2[i] = getqt(a, 2*i+9)
      y2[i] = getqt(a, 2*i+10)
      x3[i] = getqt(a, 2*i+15)
      y3[i] = getqt(a, 2*i+16)
      x2[i]:symshift(x[i], g.g, 1) -- D_1
      x2[i]:rephase(make_phase_term(0,0,0,0),{0,0,0,wallslice[srcnum]}) -- eta_1. Tech. useless.
      x3[i]:symshift(x2[i], g.g, 2) -- D_1
      x3[i]:rephase(make_phase_term(1,0,0,0),{0,0,0,wallslice[srcnum]}) -- eta_2.
      x2[i]:symshift(x3[i], g.g, 3)
      x2[i]:rephase(make_phase_term(1,1,0,0),{0,0,0,wallslice[srcnum]}) -- eta_3.
      a:solve({y2[i]}, x2[i], {mass}, resid, "all", opts, 0)
      y3[i]:symshift(y2[i], g.g, 1)
      y3[i]:rephase(make_phase_term(0,0,0,0),{0,0,0,wallslice[srcnum]}) -- eta_1. Tech. useless.
      y2[i]:symshift(y3[i], g.g, 2)
      y2[i]:rephase(make_phase_term(1,0,0,0),{0,0,0,wallslice[srcnum]}) -- eta_2.
      y3[i]:symshift(y2[i], g.g, 3)
      y3[i]:rephase(make_phase_term(1,1,0,0),{0,0,0,wallslice[srcnum]}) -- eta_3.
      y3[i]:rephase(make_phase_term(1,1,1,1), {0,0,0,wallslice[srcnum]}) -- for anti-staggered-fermion
      t[i] = y[i]:Re_dot(y3[i], "timeslices")
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 --Compensate for shifted wall source
	ps[j] = t[i][j_real]/(#wallslice) + (ps[j] or 0)
      end
      printf("End color.\n");
      io.stdout:flush();
    end 
    printf("End point source.\n");
  end
  return {pion5=p5, rhoi=ri, pionS=ps}
end

-- Fully reproduces MILC's nl_spectrum function.
function actmt.pions_wall(a, g, wallslice, mass, resid, opts)
  a:set(g, 2) -- performs the nHYP smearing.
  local nslice = #wallslice;
  local t = {}

  -- These variables hold the correlators and are returned at the end of the function.
  p5,p5_g4,pion_ps_ck, pion_4_ck, pion_i5,pion_ij = {},{},{},{},{},{}
  nucleon, nucleon_ck, delta = {}, {}, {} -- Yeah, yeah, I know it's not a pion.

  -- Since we have gauge fixed, we don't need to use the gauge field in the parallel transporter.
  local sg = qopqdp.gauge()
  sg:unit();

  -- Iterate over all wanted sources. 
  for srcnum=1,nslice do

    printf("Start wall source %i at t=%i.\n", srcnum, wallslice[srcnum]);
    -- Prepare some space. We look forward and know we need some propagators
    -- for all three colors for some nucleon measurements. In other cases, 
    -- we don't need to save them for nucleons, and we reuse memory in that case.

    -- These are only needed once per color.
    local even_src, odd_src, odd_soln, o_gupta, Do_gupta, Dq_gupta;

    -- These are needed for the nucleons. 
    local even_soln, q_gupta = {}, {}
    local Dq_gupta_all = {{},{},{}};

    -- These fill some temporary space.
    local temp1, temp2;

    -- Allocate everything in advance. This avoids over-riding issues.
    even_src = getqt(a, 1);
    odd_src = getqt(a, 2);
    odd_soln = getqt(a, 3);
    o_gupta = getqt(a, 4);
    Do_gupta = getqt(a, 5);
    Dq_gupta = getqt(a, 6);
    for i=1,3 do
      even_soln[i] = getqt(a, 7+2*(i-1));
      q_gupta[i] = getqt(a, 7+2*(i-1)+1);
    end
    temp1 = getqt(a, 13);
    temp2 = getqt(a, 14);
    for i=1,3 do
      for j=1,3 do
	Dq_gupta_all[i][j] = getqt(a, 15+3*(i-1)+(j-1));
      end
    end

    -- Loop over all colors. This is for meson measurements, and preparing
    -- for nucleon measurements.
    for i=1,3 do
      printf("Start Color %i.\n", i)
      io.stdout:flush()

      -- Obtain all of the solutions we need.
      -- Prepare the even source.
      even_src:zero()
      even_src:wall(wallslice[srcnum], 0, i, -0.125) -- Even wall. Norm matches MILC.

      -- Prepare the odd source.
      odd_src:zero()
      odd_src:wall(wallslice[srcnum], 1, i, -0.125) -- Odd wall.

      -- Invert on the even source. Just specifying the "even" inverter doesn't work
      -- when it comes to reproducing MILC. We note we need each individual even solution
      -- for the nucleon.
      printf("even_src norm2: %g\n", even_src:norm2())
      a:solve({even_soln[i]}, even_src, {mass}, resid, "all", opts, 0)

      -- Invert on the odd source.
      a:solve({odd_soln}, odd_src, {mass}, resid, "all", opts, 0)
      --printf("Solved odd wall.\n");

      -- First off, let's follow MILC in computing the local pions. 
      -- These are the Goldstone pion, and the extra gamma_4 pion. (gamma_4 gamma_5).

      do -- Pion 0_A^(-+)
	t[i] = even_soln[i]:norm2("timeslices")
	for j=1,#t[i] do
	  j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 --Compensate for shifted wall source
	  p5[j] = t[i][j_real]/nslice + (p5[j] or 0)
	end
      end

      do -- Pion 0_A^(-+) w/ extra gamma_4
	temp1:set(even_soln[i]); -- Make a copy. This is for rephasing.
	temp1:rephase(make_phase_term(1,1,1,0), {0,0,0,wallslice}) -- (-1)^(x+y+z)
	t[i] = even_soln[i]:Re_dot(temp1, "timeslices");
	for j=1,#t[i] do
	  j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 --Compensate for shifted wall source
	  p5_g4[j] = t[i][j_real]/nslice + (p5_g4[j] or 0)
	end
      end

      -- In the language of Gupta et al, the even_soln is "o+q".
      -- The odd_soln is "q-o". We can thus reconstruct Gupta's o, q.
      q_gupta[i]:zero()
      o_gupta:zero()
      q_gupta[i]:combine({even_soln[i], odd_soln},{1.0,1.0})
      o_gupta:combine({even_soln[i], odd_soln},{1.0,-1.0})

      -- We next prepare symmetric shifted values in the z direction.
      -- MILC additionally zeroes out all values that satisfy z%2==1. We do that
      -- in a slick way using rephasing and adding.

      -- First, symmetric shift the odd solution in the z direction.
      Do_gupta:symshift(o_gupta, sg, 3)

      -- Next, kill off odd z values.
      temp1:set(Do_gupta);
      temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,wallslice})
      temp2:set(Do_gupta);
      Do_gupta:combine({temp1, temp2}, {0.5,0.5})

      -- Do the same thing for the even solution.
      Dq_gupta:symshift(q_gupta[i], sg, 3)
      temp1:set(Dq_gupta);
      temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,wallslice})
      temp2:set(Dq_gupta);
      Dq_gupta:combine({temp1, temp2}, {0.5,0.5})

      -- Using these states, we can compute the one-link separated states
      -- denoted pion_i5 and pion_ij.

      -- First, pion_i5. This requires no phasing, just q Dq - o Do. 
      t[i] = q_gupta[i]:Re_dot(Dq_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	pion_i5[j] = t[i][j_real]/nslice + (pion_i5[j] or 0)
      end

      t[i] = o_gupta:Re_dot(Do_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	pion_i5[j] = -t[i][j_real]/nslice + (pion_i5[j] or 0)
      end

      -- Next, pion_ij. This requires rephasing by a factor of (-1)^(x+y+z),
      -- then we take o Dq - q Do.
      Dq_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, wallslice});
      Do_gupta:rephase(make_phase_term(1,1,1,0), {0,0,0,wallslice});

      t[i] = o_gupta:Re_dot(Dq_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	pion_ij[j] = t[i][j_real]/nslice + (pion_ij[j] or 0)
      end

      t[i] = q_gupta[i]:Re_dot(Do_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	pion_ij[j] = -t[i][j_real]/nslice + (pion_ij[j] or 0)
      end

      -- We can also do checks on the goldstone boson pion and the
      -- local gamma4 pion. To be more specific, these are the states that
      -- Goltermann defines! The measurements before are slick shortcuts/consistency
      -- checks MILC does.

      -- First, the Goldstone boson pion. This is q q + o o.
      t[i] = q_gupta[i]:norm2("timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	pion_ps_ck[j] = t[i][j_real]/nslice + (pion_ps_ck[j] or 0)
      end

      t[i] = o_gupta:norm2("timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	pion_ps_ck[j] = t[i][j_real]/nslice + (pion_ps_ck[j] or 0)
      end

      -- Next, the gamma4 pion. This is q o with a (-1)^(x+y+z) phase factor.
      temp1:set(o_gupta);
      temp1:rephase(make_phase_term(1,1,1,0), {0,0,0,wallslice})
      t[i] = q_gupta[i]:Re_dot(temp1, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	pion_4_ck[j] = t[i][j_real]/nslice + (pion_4_ck[j] or 0)
      end

      printf("End Color.\n");
      io.stdout:flush();
    end

    -- Now that we're done with all the mesons, we can use the even and q
    -- solutions we saved to build Nucleons. For nucleons, one just sums over
    -- the primary site in the hypercube. We're going to reuse the trick we
    -- used to zero out sites with z%2==1 to zero out all sites with x%2 or
    -- y%2 == 1 as well!

    do -- Nucleon!
      -- Perform the zeroing.
      for i=1,3 do
	temp1:set(even_soln[i]);
	temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,wallslice})
	temp2:set(even_soln[i]); -- reuse variable
	even_soln[i]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd x.

	temp1:set(even_soln[i]);
	temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,wallslice})
	temp2:set(even_soln[i]); -- reuse variable
	even_soln[i]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd y.

	temp1:set(even_soln[i]);
	temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,wallslice})
	temp2:set(even_soln[i]); -- reuse variable
	even_soln[i]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd z.
      end

      -- Now that we've zeroed out non-primary sites, we perform the epsilon
      -- contraction. This takes the determinant of a matrix where each column
      -- corresponds to one of the three color vectors. 
      t[1] = even_soln[1]:epsContract({even_soln[2],even_soln[3]}, "timeslices");
      for j=1,#t[1] do
	j_real = (j+wallslice[srcnum]-1)%(#t[1])+1 -- Compensate for shifted wall source
	-- Compensate for a baryon wrap-around effect. 
	if (((math.floor((j+wallslice[srcnum]-1)/(#t[1]))-math.floor((wallslice[srcnum]-1)/(#t[1])))%2) == 0) then
	  nucleon[j] = -t[1][j_real]/nslice + (nucleon[j] or 0)
	else
	  nucleon[j] = t[1][j_real]/nslice + (nucleon[j] or 0)
	end
      end  
    end

    do -- Delta!
      -- An operator which measures the Delta baryon can be found in the end of
      -- Golterman's Lattice baryon paper. The method here is the operator defined
      -- by equation 6.3. 

      -- We need symmetric shifts in all directions for all colors.
      -- We also need to zero out all non-primary lattice sites. 
      for i=1,3 do -- For all colors
	for j=1,3 do -- For all directions

	  -- Perform the symmetric shift.
	  Dq_gupta_all[i][j]:symshift(q_gupta[i], sg, j);

	  -- Zero out non-primary sites.
	  temp1:set(Dq_gupta_all[i][j]);
	  temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,wallslice})
	  temp2:set(Dq_gupta_all[i][j]); 
	  Dq_gupta_all[i][j]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd x.

	  temp1:set(Dq_gupta_all[i][j]);
	  temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,wallslice})
	  temp2:set(Dq_gupta_all[i][j]);
	  Dq_gupta_all[i][j]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd y.

	  temp1:set(Dq_gupta_all[i][j]);
	  temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,wallslice})
	  temp2:set(Dq_gupta_all[i][j]); 
	  Dq_gupta_all[i][j]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd z.
	end
      end

      -- While there's slicker ways to do this, I explicitly construct
      -- an antisymmetric tensor. 
      eps_symbol = {};
      for i=1,3 do
	eps_symbol[i] = {};
	for j=1,3 do
	  eps_symbol[i][j] = {};
	  for k=1,3 do
	    eps_symbol[i][j][k] = 0;
	  end
	end
      end
      eps_symbol[1][2][3] = 1;
      eps_symbol[2][3][1] = 1;
      eps_symbol[3][1][2] = 1;
      eps_symbol[2][1][3] = -1;
      eps_symbol[3][2][1] = -1;
      eps_symbol[1][3][2] = -1;
      t[1] = {}; 
      for k=1,3 do -- All elements
	for l = 1,3 do -- of the
	  for m = 1,3 do -- epsilon symbol (color)
	    -- Don't waste time with zero elements!
	    if not (eps_symbol[k][l][m] == 0) then 
	      -- Construct the proper epsilon contraction.
	      t[1] = Dq_gupta_all[k][1]:epsContract({Dq_gupta_all[l][2],Dq_gupta_all[m][3]}, "timeslices");
	      for j=1,#t[1] do
		j_real = (j+wallslice[srcnum]-1)%(#t[1])+1 -- Compensate for shifted wall source
		if (((math.floor((j+wallslice[srcnum]-1)/(#t[1]))-math.floor((wallslice[srcnum]-1)/(#t[1])))%2) == 0) then
		  delta[j] = -eps_symbol[k][l][m]*t[1][j_real]/nslice + (delta[j] or 0)
		else
		  delta[j] = eps_symbol[k][l][m]*t[1][j_real]/nslice + (delta[j] or 0)
		end
	      end  
	    end
	  end
	end
      end

    end
    do -- Nucleon check! This is similar to how local meson states have a check.
      -- Zeroing.
      for i=1,3 do
	temp1:set(q_gupta[i]);
	temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,wallslice})
	temp2:set(q_gupta[i]); 
	q_gupta[i]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd x.

	temp1:set(q_gupta[i]);
	temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,wallslice})
	temp2:set(q_gupta[i]); 
	q_gupta[i]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd y.

	temp1:set(q_gupta[i]);
	temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,wallslice})
	temp2:set(q_gupta[i]); 
	q_gupta[i]:combine({temp2, temp1}, {0.5,0.5}) -- Kills off odd z.
      end

      -- Contract.
      t[1] = q_gupta[1]:epsContract({q_gupta[2], q_gupta[3]}, "timeslices");
      for j=1,#t[1] do
	j_real = (j+wallslice[srcnum]-1)%(#t[1])+1 -- Compensate for shifted wall source
	if (((math.floor((j+wallslice[srcnum]-1)/(#t[1]))-math.floor((wallslice[srcnum]-1)/(#t[1])))%2) == 0) then
	  nucleon_ck[j] = -t[1][j_real]/nslice + (nucleon_ck[j] or 0)
	else
	  nucleon_ck[j] = t[1][j_real]/nslice + (nucleon_ck[j] or 0)
	end
      end  
    end
    printf("End wall source.\n");
  end

  return {pion5=p5, pion5_gamma4=p5_g4, pion5_ck=pion_ps_ck, pion5_gamma4_ck=pion_4_ck,
	  pion_i5=pion_i5, pion_ij=pion_ij, nucleon=nucleon, nucleon_ck=nucleon_ck, delta=delta}
end

function actmt.s4_broken_observe(a, g, mass, resid, opts, npbp)
  local s4_gauge_observables_even, s4_gauge_observables_odd = {}, {};
  local s4_gauge_observables_plaq = 0;

  a:set(g, 2) -- performs the nHYP smearing.

  -- Get smeared gauge now.
  local sg = smearGauge(g, a.smear)

  -- Get the gauge observables.
  s4_gauge_observables_plaq, s4_gauge_observables_even, s4_gauge_observables_odd = sg:s4Gauge();

  -- Next, let's get the fermionic observables.
  -- We start with pbp.
  a:set(g, 2)
  local x = getqt(a, 1)
  local y = getqt(a, 2)

  -- Repeat pbp measurements as needed.
  s4_ferm_e, s4_ferm_o = {}, {};

  for k,v in ipairs(mass) do

    s4_ferm_e[k] = {};
    s4_ferm_o[k] = {};

    for i=1,npbp do

      x:randomU1()
      a:solve({y}, x, {v}, resid, "all", opts, 0)

      -- We just pass the color vector and the gauge field
      -- onto the C level.

      s4_ferm_observables_even, s4_ferm_observables_odd = y:s4Ferm(x, sg);

      local normalization = x:norm2();
      s4_ferm_e[k][1] = v * s4_ferm_observables_even[1] / normalization + (s4_ferm_e[k][1] or 0);
      s4_ferm_o[k][1] = v *s4_ferm_observables_odd[1] / normalization + (s4_ferm_o[k][1] or 0);
      for j=2,6 do
	s4_ferm_e[k][j] = s4_ferm_observables_even[j] / normalization * 2 + (s4_ferm_e[k][j] or 0);
	s4_ferm_o[k][j] = s4_ferm_observables_odd[j] / normalization * 2 + (s4_ferm_o[k][j] or 0);
      end

    end

    for i=1,6 do
      s4_ferm_e[k][i] = s4_ferm_e[k][i] / npbp;
      s4_ferm_o[k][i] = s4_ferm_o[k][i] / npbp;
    end

  end

  return {s4_g_plaq = s4_gauge_observables_plaq, s4_g_even = s4_gauge_observables_even,
	  s4_g_odd = s4_gauge_observables_odd, s4_f_even = s4_ferm_e, s4_f_odd = s4_ferm_o}
end
