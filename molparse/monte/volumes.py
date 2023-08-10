
import plotly.graph_objects as go
import mgo
import numpy as np

class CompoundVolume:

	def __init__(self,reference=None):
		self.volumes = []
		self.reference = reference

	# def cull_overlapping_volumes(self):

	def add_volume(self,volume):
		self.volumes.append(volume)

	def plot(self,distance,fig=None,show=False,samples=20):
		if fig is None:
			fig = go.Figure()
		
		for vol in self.volumes:
			fig = vol.plot(distance=distance,samples=samples,fig=fig,show=False)

		if show:
			fig.show()

		return fig

	def is_inside(self,point):

		for vol in self.volumes:

			if vol.is_inside(point):
				return True

		return False

class Sphere:

	def __init__(self,centre,radius):
		
		self.centre = np.array(centre)
		self.radius = np.array(radius)

	def is_inside(self,point):
		return np.linalg.norm(point - self.centre) < self.radius

class CappedCone:

	def __init__(self,origin,target,radius,name=None):

		# parameters
		self.name = name
		self.origin = np.array(origin)
		self.target = np.array(target)
		self.radius = radius

		# derived internals
		self.vec_origin_target = self.target - self.origin
		self.dist_origin_target = np.linalg.norm(self.vec_origin_target)
		self.unit_origin_target = self.vec_origin_target / self.dist_origin_target
		self.theta = np.arctan(self.radius/self.dist_origin_target)
		self.cos_theta = np.cos(self.theta)
		self.start = self.dist_origin_target * self.cos_theta * self.cos_theta

		# child volumes
		self.cap_radius = self.dist_origin_target * np.sin(self.theta)
		self.cap_sphere = Sphere(target,self.cap_radius)

		# plotting variables
		self.vec_origin_point = None
	
	def is_inside(self,point):

		self._is_inside_cap = self.cap_sphere.is_inside(point)

		if self._is_inside_cap:
			return True

		self.vec_origin_point = point - self.origin
		self.dist_origin_point = np.linalg.norm(self.vec_origin_point)
		self.unit_origin_point = self.vec_origin_point / self.dist_origin_point

		if self.dist_origin_point < self.start:
			return False

		self.origin_point_dot_origin_target = np.dot(self.unit_origin_point,self.unit_origin_target)

		if self.origin_point_dot_origin_target > self.cos_theta:
			# print(f'{self.origin_point_dot_origin_target=}')
			# print(f'{self.cos_theta=}')
			return True

		return False

	def test_random_point(self,n=1,distance=5):

		fig = self.plot(distance)

		for i in range(n):

			# point = random_point(self.origin,distance)
			point = np.array([-2.127,2.62,-17.7])

			# print(f'{point=}')

			is_inside = self.is_inside(point)

			# print(f'{is_inside=}')
			# print(f'{self._is_inside_cap=}')

			if self.vec_origin_point is not None:
				trace = mgo.vector_trace(self.origin,self.vec_origin_point,name=f'{is_inside=}')
				fig.add_trace(trace)

		fig.show()

	def plot(self,distance,samples=20,fig=None,show=False,verbosity=1):

		if fig is None:
			fig = go.Figure()

		trace = mgo.cone_trace(self.origin, self.target, self.radius, distance, samples=samples, name=self.name, start_at_target=True, plot_caps=True, verbosity=verbosity-1)
		fig.add_trace(trace)

		if show:
			fig.show()

		return fig

def random_point(origin,radius):

	while True:

		x = np.random.uniform(-1.0,1.0)
		y = np.random.uniform(-1.0,1.0)
		z = np.random.uniform(-1.0,1.0)

		if np.linalg.norm([x,y,z]) < 1:
			return origin + radius * np.array([x,y,z])
