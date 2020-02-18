
# load the correct modules for ASE/PoV-Ray
import os
import mout # https://github.com/mwinokan/MPyTools
os.system('source $MWAMPPATH/load_pov.sh 2> /dev/null')
mout.out("ASE & PoV-Ray dependencies loaded.",printScript=True,)

# import functions
from .povplot import makePovImage
from . import styles