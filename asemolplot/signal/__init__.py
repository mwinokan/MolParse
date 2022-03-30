
from .modify import runningAverage
from .modify import differentiate
# from .2d import peakFinder

from .io import parseDat
from .io import csv_strip

from .peaks import peakFinder

import socket
if not 'eureka' in socket.gethostname():
	from .fitwizard import fitwizard
	