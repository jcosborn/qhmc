-- How to measure and interpret the results of pions_wall:

-- The following code loops over four sources and prints out all outputs from the pions_wall function.
-- "act" is implicitly an asqtad fermions action object.

-- Wall pions!
for j=2,26,8 do
	pions_wall = act.f:pions_wall(G,{j},mass,1e-6, {prec = prec, restart = 5000})
	--pions_wall = act.f:pions_wall(G,{2,10,18,26},mass,1e-6, {prec = prec, restart = 5000})

	printf("Source %i\n", j);
	printf("Local Pions\n");
	printf("\tPion5\t\tPion5_4\n");

	for i = 1,#(pions_wall.pion5) do
	   printf("%i\t%.6e\t%.6e\n", i-1, pions_wall.pion5[i], pions_wall.pion5_gamma4[i]);
	end
	
	printf("Local Baryon\n");
	printf("\tNucleon\n");

	for i = 1,#(pions_wall.pion5) do
	   printf("%i\t%.6e\n", i-1, pions_wall.nucleon[i]);
	end

	printf("\nNonlocal Pions + Check\n")
	printf("\tPion5\tPion5_4\tPion_i5\tPion_ij\n");

	for i = 1,#(pions_wall.pion5) do
	   printf("%i\t%.6e\t%.6e\t%.6e\t%.6e\n", i-1, pions_wall.pion5_ck[i], pions_wall.pion5_gamma4_ck[i], pions_wall.pion_i5[i], pions_wall.pion_ij[i]);
	end
	
	printf("Nonlocal Baryons\n");
	printf("\tNucleon\tDelta\n");

	for i = 1,#(pions_wall.pion5) do
	   printf("%i\t%.6e\t%.6e\n", i-1, pions_wall.nucleon_ck[i], pions_wall.delta[i]);
	end

end

-- Alternatively, to perform measurements on multiple wall sources and automatically average over them, call something
-- along the lines of:

pions_wall = act.f:pions_wall(G,{2, 10, 18, 26},mass,1e-6, {prec = prec, restart = 5000})

-- Where the {2, 10, 18, 26} implies placing wall sources on all of these time slices and averaging over them.

-----------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------

-- How to measure and interpret the results of s4_broken_observe:

-- Arguments:
-- 1. gauge action, NOT FIELD.
-- 2. The mass you want to use to measure pbp. This supports more than one mass, but it uses the same number of repetitions and stopping conditions for all.
-- 3. Stopping condition of CG.
-- 4. Other parameters of CG. prec = 1 for single precision, 2 for double, and restart is the restart value.
-- 5. The number of stoichastic pbp measurements to perform. I just did 1, but its your call. It automatically averages over all of them.
-- Note: the gauge observables use the nHYP smeared gauge links, not the bare ones. If they don't agree with what other plaquette measurements give you, that's why!
local masses = {mass, 2*mass};
local s4_check = act.f:s4_broken_observe(G, masses, 1e-6, {prec = prec, restart = 5000}, 1)

printf("MEASplaq_ss %.12e\n", s4_check.s4_g_plaq);
printf("MEASplaq_t even %.12e odd %.12e\n", s4_check.s4_g_even[1], s4_check.s4_g_odd[1]);
printf("MEASplaq_x even %.12e odd %.12e\n", s4_check.s4_g_even[2], s4_check.s4_g_odd[2]);
printf("MEASplaq_y even %.12e odd %.12e\n", s4_check.s4_g_even[3], s4_check.s4_g_odd[3]);
printf("MEASplaq_z even %.12e odd %.12e\n", s4_check.s4_g_even[4], s4_check.s4_g_odd[4]);
printf("MEASplaq_a even %.12e odd %.12e\n", s4_check.s4_g_even[5], s4_check.s4_g_odd[5]);

for i,v in ipairs(masses) do
	printf("MEASpbp_all m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][1].r, s4_check.s4_f_even[i][1].i, s4_check.s4_f_odd[i][1].r, s4_check.s4_f_odd[i][1].i);
	printf("MEASpbp_1 m %.4e even %.12e %.12e odd %.12e %.12e\n",  v,s4_check.s4_f_even[i][2].r, s4_check.s4_f_even[i][2].i, s4_check.s4_f_odd[i][2].r, s4_check.s4_f_odd[i][2].i);
	printf("MEASpbp_t m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][3].r, s4_check.s4_f_even[i][3].i, s4_check.s4_f_odd[i][3].r, s4_check.s4_f_odd[i][3].i);
	printf("MEASpbp_x m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][4].r, s4_check.s4_f_even[i][4].i, s4_check.s4_f_odd[i][4].r, s4_check.s4_f_odd[i][4].i);
	printf("MEASpbp_y m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][5].r, s4_check.s4_f_even[i][5].i, s4_check.s4_f_odd[i][5].r, s4_check.s4_f_odd[i][5].i);
	printf("MEASpbp_z m %.4e even %.12e %.12e odd %.12e %.12e\n", v, s4_check.s4_f_even[i][6].r, s4_check.s4_f_even[i][6].i, s4_check.s4_f_odd[i][6].r, s4_check.s4_f_odd[i][6].i);
end
