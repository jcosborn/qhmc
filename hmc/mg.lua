local mgIsSet = false

function mgSetup(w)
  if mgIsSet then return end
  mgIsSet = true
  if noMG then return end
  local kappa = 0.5/(mass+4)
  --local block = {2,2,2,2}
  local block = {4,4,4,4}
  local lat = {}
  for i=1,#latsize do lat[i] = latsize[i]/block[i] end
  w:mgSet({[-1] = { nlevels=1 }})
  w:mgSet(
    {
      [-2] = { verbose=0 },
      [-1] = { verbose=1, kappa=kappa, kappanv=kappa, itmax=100, ngcr=8 },
      [0] = { lattice=lat, nvecs=24,
	      setup_res=0.4, setup_maxit=100, setup_change_fac=0.5,
	      npre=0, npost=4, scale=1, cres=0.1, itmax=100, ngcr=8 }
    }
  )
  w:mgSetup()
end
