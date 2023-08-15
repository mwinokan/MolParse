
import numpy as np

def random_point_spherical(origin,radius):

	while True:

		x = np.random.uniform(-1.0,1.0)
		y = np.random.uniform(-1.0,1.0)
		z = np.random.uniform(-1.0,1.0)

		if np.linalg.norm([x,y,z]) < 1:
			return origin + radius * np.array([x,y,z])

