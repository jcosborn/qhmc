local function getqt(s, i)
  local q = s:quark()
  q:zero()
  return q
end

-- Use 0 for no phase, 1 for phase.
local function make_phase_term(xsign,ysign,zsign,tsign)
  return xsign + 2*ysign + 4*zsign + 8*tsign
end

-------------------------------------------------------------------------
-------------------------------------------------------------------------
function staggeredPions(s, g, wallslice, mass, resid, opts)
  local x,y,z,t,p5,ri,ps = {},{},{},{},{},{},{}
  local x2,y2 = {},{}
  local x3,y3 = {},{}
  local ps = {}
  for srcnum=1,#wallslice do

    printf("Start point source %i at t=%i.\n", srcnum, wallslice[srcnum]);
    for i=1,3 do
      printf("Start Color %i.\n", i)
      io.stdout:flush()
      x[i] = getqt(s, 3*i)
      y[i] = getqt(s, 3*i+1)
      z[i] = getqt(s, 3*i+2) -- For the rho sign change.
      x[i]:zero()
      x[i]:point({0,0,0,wallslice[srcnum]},i,1,0)
      s:solve({y[i]}, x[i], {mass}, resid, "all", opts)
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
      x2[i] = getqt(s, 2*i+9)
      y2[i] = getqt(s, 2*i+10)
      x3[i] = getqt(s, 2*i+15)
      y3[i] = getqt(s, 2*i+16)
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


-------------------------------------------------------------------------
-- Fully reproduces MILC's nl_spectrum function.
-------------------------------------------------------------------------
function staggeredPionsWall(s, g, wallslice, mass, resid, opts)
  local nslice = #wallslice;
  local t = {}

  -- These variables hold the correlators and are returned at the end of the function.
  p5,p5_g4,pion_ps_ck, pion_4_ck, pion_i5,pion_ij = {},{},{},{},{},{}
  rho_0, rho_is, rho_ij, rho_i5 = {}, {}, {}, {}
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
    even_src = getqt(s, 1);
    odd_src = getqt(s, 2);
    odd_soln = getqt(s, 3);
    o_gupta = getqt(s, 4);
    Do_gupta = getqt(s, 5);
    Dq_gupta = getqt(s, 6);
    for i=1,3 do
      even_soln[i] = getqt(s, 7+2*(i-1));
      q_gupta[i] = getqt(s, 7+2*(i-1)+1);
    end
    temp1 = getqt(s, 13);
    temp2 = getqt(s, 14);
    for i=1,3 do
      for j=1,3 do
	Dq_gupta_all[i][j] = getqt(s, 15+3*(i-1)+(j-1));
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

--[[
      even_src:zero()
      even_src:set(-0.125, i, "timeslice"..wallslice[srcnum])
      odd_src:zero()
      odd_src:set(even_src, "odd")
      even_src:zero("odd")
--]]

      -- Invert on the even source. Just specifying the "even" inverter doesn't work
      -- when it comes to reproducing MILC. We note we need each individual even solution
      -- for the nucleon.
      printf("even_src norm2: %g\n", even_src:norm2())
      s:solve({even_soln[i]}, even_src, {mass}, resid, "all", opts)

      -- Invert on the odd source.
      s:solve({odd_soln}, odd_src, {mass}, resid, "all", opts)
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

      -- NEW ADDITION! We can now grab a few rho states.
      -- First, x shifted. This lets us grab:
      -- 1. rho_0: gamma_1 gamma_4 x taste_4
      -- 2. rho_is: gamma_1 x 1
      -- These names match a MILC convention.
      -- First, symmetric shift the odd solution in the x direction,
      -- then zero out x-odd sites. 
      Do_gupta:symshift(o_gupta, sg, 1)
      temp1:set(Do_gupta);
      temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,wallslice})
      temp2:set(Do_gupta);
      Do_gupta:combine({temp1, temp2}, {0.5,0.5})

      -- Do the same thing for the even solution.
      Dq_gupta:symshift(q_gupta[i], sg, 1)
      temp1:set(Dq_gupta);
      temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,wallslice})
      temp2:set(Dq_gupta);
      Dq_gupta:combine({temp1, temp2}, {0.5,0.5})

      -- Using these states, we can compute the one-link separated states
      -- denoted rho_0 and rho_is.
      -- First, rho_0. This requires no phasing, just q Dq - o Do. 
      t[i] = q_gupta[i]:Re_dot(Dq_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	rho_0[j] = t[i][j_real]/nslice + (rho_0[j] or 0)
      end

      t[i] = o_gupta:Re_dot(Do_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	rho_0[j] = -t[i][j_real]/nslice + (rho_0[j] or 0)
      end

      -- Next, rho_is. This requires rephasing by a factor of (-1)^(x+y+z),
      -- then we take o Dq - q Do.
      Dq_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, wallslice});
      Do_gupta:rephase(make_phase_term(1,1,1,0), {0,0,0,wallslice});

      t[i] = o_gupta:Re_dot(Dq_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	rho_is[j] = t[i][j_real]/nslice + (rho_is[j] or 0)
      end

      t[i] = q_gupta[i]:Re_dot(Do_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	rho_is[j] = -t[i][j_real]/nslice + (rho_is[j] or 0)
      end

      -- First, y shifted. This lets us grab:
      -- 1. rho_ij: gamma_3 x taste_1 taste_4 taste_5
      -- 2. rho_i5: gamma_3 gamma_4 x taste_1 taste_5
      -- I couldn't find a MILC convention for these states,
      -- so I picked what I did based on looking at the taste
      -- structure of the pion states.

      -- First, symmetric shift the odd solution in the y direction,
      -- then zero out y-odd sites. 
      Do_gupta:symshift(o_gupta, sg, 2)
      temp1:set(Do_gupta);
      temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,wallslice})
      temp2:set(Do_gupta);
      Do_gupta:combine({temp1, temp2}, {0.5,0.5})

      -- Do the same thing for the even solution.
      Dq_gupta:symshift(q_gupta[i], sg, 2)
      temp1:set(Dq_gupta);
      temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,wallslice})
      temp2:set(Dq_gupta);
      Dq_gupta:combine({temp1, temp2}, {0.5,0.5})

      -- Using these states, we can compute the one-link separated states
      -- denoted rho_0 and rho_is.
      -- First, rho_0. This requires no phasing, just q Dq - o Do. 
      t[i] = q_gupta[i]:Re_dot(Dq_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	rho_ij[j] = t[i][j_real]/nslice + (rho_ij[j] or 0)
      end

      t[i] = o_gupta:Re_dot(Do_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	rho_ij[j] = -t[i][j_real]/nslice + (rho_ij[j] or 0)
      end

      -- Next, rho_is. This requires rephasing by a factor of (-1)^(x+y+z),
      -- then we take o Dq - q Do.
      Dq_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, wallslice});
      Do_gupta:rephase(make_phase_term(1,1,1,0), {0,0,0,wallslice});

      t[i] = o_gupta:Re_dot(Dq_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	rho_i5[j] = t[i][j_real]/nslice + (rho_i5[j] or 0)
      end

      t[i] = q_gupta[i]:Re_dot(Do_gupta, "timeslices");
      for j=1,#t[i] do
	j_real = (j+wallslice[srcnum]-1)%(#t[i])+1 -- Compensate for shifted wall source
	rho_i5[j] = -t[i][j_real]/nslice + (rho_i5[j] or 0)
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
  pion_i5=pion_i5, pion_ij=pion_ij, rho_0 = rho_0, rho_is=rho_is, rho_ij=rho_ij, rho_i5=rho_i5,
  nucleon=nucleon, nucleon_ck=nucleon_ck, delta=delta}
end

function staggeredS4Broken(s, sg, mass, resid, opts, npbp)
  local s4_gauge_observables_even, s4_gauge_observables_odd = {}, {};
  local s4_gauge_observables_plaq = 0;

  -- Get the gauge observables.
  local s4_gauge_observables_plaq, s4_gauge_observables_even, s4_gauge_observables_odd = sg:s4Gauge();

  -- Next, let's get the fermionic observables.
  -- We start with pbp.
  local x = getqt(s, 1)
  local y = getqt(s, 2)
  -- Repeat pbp measurements as needed.
  local s4_ferm_e, s4_ferm_o = {}, {};

  for k,v in ipairs(mass) do
    s4_ferm_e[k] = {};
    s4_ferm_o[k] = {};
    for i=1,npbp do
      x:randomU1()
      s:solve({y}, x, {v}, resid, "all", opts)
      -- We just pass the color vector and the gauge field
      -- onto the C level.
      local s4_ferm_observables_even, s4_ferm_observables_odd = y:s4Ferm(x,sg);
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
