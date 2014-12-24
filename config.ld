-- ldoc configuration file
title = 'QHMC (FUEL) API documentation'
project = 'QHMC (FUEL)'
format = 'markdown'
--convert_opt = true
add_language_extension('cdoc','c') 
merge = true

description = [[
QHMC is an application for Lattice Field Theory using the
[Lua](http://www.lua.org)
scripting language on top of the
[USQCD SciDAC libraries](http://usqcd-software.github.io).  
The main website is [here](http://jcosborn.github.io/qhmc).  

** This documentation is still very incomplete. **  
** Refer to the source and the scripts included with it for a more complete
reference. **
]]

readme = {'README','AUTHORS','NEWS'}

file = {
  'lib/base/global.cdoc',
  'lib/base/qhmc_base.c',
  'lib/base/qhmc_complex.c',
  'lib/qopqdp/qhmc_qopqdp.c',
  'lib/qopqdp/smear.c',
}

examples = {
  'test/point.lua',
  'test/norm2.lua',
  'test/subset.lua',
  'test/random.lua',
  'test/fieldIO.lua',
  'test/gauge.lua',
}
