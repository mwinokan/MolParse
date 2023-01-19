
from .modify import runningAverage
from .modify import differentiate
# from .2d import peakFinder

from .io import parseDat
from .io import csv_strip

from .peaks import peakFinder, closest_index, closest_value

# the fitwizard should be loaded only as needed via: from ase.signal 

# import socket
# if not 'eureka' in socket.gethostname():
# 	if not 'login' in socket.gethostname():
# 		if not 'node' in socket.gethostname():
# 			from .fitwizard import fitwizard
	