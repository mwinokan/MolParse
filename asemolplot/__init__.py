import mout # https://github.com/mwinokan/MPyTools
import module # https://github.com/mwinokan/MPyTools

# load the correct modules for ASE/PoV-Ray
module.module('--expert','load','Boost/1.63.0-intel-2017a-Python-2.7.13')
module.module('--expert','load','zlib/1.2.8-intel-2016a')
module.module('--expert','load','libpng/1.6.24-intel-2016a')
module.module('--expert','load','libjpeg-turbo/1.5.0-intel-2016a')
module.module('--expert','load','LibTIFF/4.0.6-intel-2016a')
module.module('--expert','load','anaconda3/2019.03')

mout.out("PoV-Ray dependencies loaded.",
         printScript=True,)

# import functions
from .povplot import makePovImage
from . import styles
