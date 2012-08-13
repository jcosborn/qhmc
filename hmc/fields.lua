local fieldsmt = {}

function setupfields(a, p)
  local fields = {}
  fields.a = a
  fields.G = a.g:gaugeNew()
  fields.GSave = a.g:gaugeNew()
  fields.F = a.g:forceNew()
  fields.vol = p.vol
  fields.npseudo = p.npseudo
  return setmetatable(fields, {__index=fieldsmt})
end

function fieldsmt.save(f)
  f.GSave:set(f.G)
  f.F:random()
  f.a.f:refresh(f.G)
end

function fieldsmt.reverse(f)
  f.Grev = f.a.g:gaugeNew()
  f.Grev:set(f.G)
  f.F.f:scale(-1)
end

function fieldsmt.endReverse(f)
  f.G:set(f.Grev)
  f.Grev = nil
end

function fieldsmt.action(f)
  local Sgq = f.a.g:action(f.G)
  local Sgp = 0.5 * f.F:norm2() - 16*f.vol
  local Sfq = f.a.f:action(f.G)
  local S = Sgq + Sgp + Sfq
  printf("Sgq: %-12.10g  Sgp: %-12.10g  Sfq: %-12.10g\n", Sgq, Sgp, Sfq)
  return S
end

function fieldsmt.accept(f)
  local devavg,devmax = f.G:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  f.G:makeSU()
  devavg,devmax = f.G:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end

function fieldsmt.reject(f)
  f.G:set(f.GSave)
end

function fieldsmt.nfields(f)
  return 1
end

function fieldsmt.nforces(f)
  return { 1 + f.npseudo }
end

--local gt = 0
function fieldsmt.updateField(f, i, eps)
  --gt = gt + eps
  --printf("Gupdate %g %g %g\n", eps, gt, 0.5*f.F:norm2()-16*vol)
  f.G:update(f.F, eps)
  --local ss,st = f.G:plaq()
  --printf("plaq %g %g\n", 3*ss, 3*st);
end

function fieldsmt.updateMomentum(f, i, tj, teps)
  for k,j in ipairs(tj) do
    if(j==1) then
      --printf("begin GFupdate %g %g\n", eps, 0.5*f.F:norm2()-16*f.vol)
      f.a.g:updateMomentum(f.F, f.G, teps[k])
      --printf("end   GFupdate %g %g\n", eps, 0.5*f.F:norm2()-16*f.vol)
      table.remove(tj, k)
      table.remove(teps, k)
    end
  end
  if #tj > 0 then
    for k=1,#tj do tj[k] = tj[k] - 1 end
    --myprint("forces ", tj, "\n")
    --printf("begin FFupdate[%i] %g %g\n", j-1, eps, 0.5*f.F:norm2()-16*f.vol)
    f.a.f:updateMomentum(f.F, f.G, teps, tj)
    --printf("end   FFupdate[%i] %g %g\n", j-1, eps, 0.5*f.F:norm2()-16*f.vol)
  end
end
