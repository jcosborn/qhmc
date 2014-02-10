--[[
table.insert(package.searchers,1,
  function(m)
    print(m)
    print(_ENV[m])
    return function(m) return _ENV[m] end
  end
)
--]]

f = require 'lfs'
for file in f.dir(".") do
  print(file)
end
