import plotly.graph_objects as go
import mgo
import numpy as np
from .rand import random_point_spherical


class CompoundVolume:
    def __init__(self, reference=None):
        self.volumes = []
        self.reference = reference

    # def cull_overlapping_volumes(self):

    def add_volume(self, volume):
        volume.index = self.num_volumes
        self.volumes.append(volume)

    def plot(self, distance, fig=None, show=False, samples=20):
        if fig is None:
            fig = go.Figure()

        for vol in self.volumes:
            fig = vol.plot(distance=distance, samples=samples, fig=fig, show=False)

        if show:
            fig.show()

        return fig

    def is_inside(self, point):

        for vol in self.volumes:

            if vol.is_inside(point):
                return True
        return False

    def is_within(self, point, within):
        for vol in self.volumes:
            if vol.is_within(point, within):
                return True
        return False

    @property
    def num_volumes(self):
        return len(self.volumes)

    def simplify(self):

        del_list = []

        for this_volume in self.volumes:
            for other_volume in self.volumes:

                if this_volume.index == other_volume.index:
                    continue

                if this_volume in del_list:
                    continue

                if other_volume in del_list:
                    continue

                print(this_volume, this_volume.index, other_volume, other_volume.index)

                if this_volume.always_inside(other_volume):
                    del_list.append(this_volume)

        print(del_list)

        for volume in reversed(self.volumes):
            if volume in del_list:
                del volume


class Sphere:
    def __init__(self, centre, radius, index=None):
        self.index = index
        self.centre = np.array(centre)
        self.radius = np.array(radius)

    def is_inside(self, point):
        return np.linalg.norm(point - self.centre) < self.radius


import mout


class CappedCone:
    def __init__(self, origin, target, radius, name=None, index=None):

        # parameters
        self.name = name
        self.origin = np.array(origin)
        self.target = np.array(target)
        self.radius = radius
        self.index = index

        # derived internals
        self.vec_origin_target = self.target - self.origin
        self.dist_origin_target = np.linalg.norm(self.vec_origin_target)
        self.unit_origin_target = self.vec_origin_target / self.dist_origin_target
        self.theta = np.arctan(self.radius / self.dist_origin_target)
        self.cos_theta = np.cos(self.theta)
        self.start = self.dist_origin_target * self.cos_theta * self.cos_theta

        # child volumes
        self.cap_radius = self.dist_origin_target * np.sin(self.theta)
        self.cap_sphere = Sphere(target, self.cap_radius)

        # plotting variables
        self.vec_origin_point = None

        self._expanded_volume = None

    def is_inside(self, point):

        self._is_inside_cap = self.cap_sphere.is_inside(point)

        if self._is_inside_cap:
            return True

        self.vec_origin_point = point - self.origin
        self.dist_origin_point = np.linalg.norm(self.vec_origin_point)

        if self.dist_origin_point == 0:
            return False

        self.unit_origin_point = self.vec_origin_point / self.dist_origin_point

        if self.dist_origin_point < self.start:
            return False

        self.origin_point_dot_origin_target = np.dot(
            self.unit_origin_point, self.unit_origin_target
        )

        if self.origin_point_dot_origin_target > self.cos_theta:
            # print(f'{self.origin_point_dot_origin_target=}')
            # print(f'{self.cos_theta=}')
            return True

        return False

    def is_within(self, point, within):

        if self._expanded_volume is None or self._expanded_volume._expansion != within:
            new_origin = self.origin - self.unit_origin_target * within / np.sin(
                self.theta
            )

            self._expanded_volume = CappedCone(
                origin=new_origin,
                target=self.target,
                radius=self.radius + within,
                name=self.name + " expanded",
            )
            self._expanded_volume._expansion = within

        return self._expanded_volume.is_inside(point)

    def test_random_point(self, n=1, distance=5):

        fig = self.plot(distance)

        for i in range(n):

            # point = random_point(self.origin,distance)
            point = np.array([-2.127, 2.62, -17.7])

            # print(f'{point=}')

            is_inside = self.is_inside(point)

            # print(f'{is_inside=}')
            # print(f'{self._is_inside_cap=}')

            if self.vec_origin_point is not None:
                trace = mgo.vector_trace(
                    self.origin, self.vec_origin_point, name=f"{is_inside=}"
                )
                fig.add_trace(trace)

        fig.show()

    def plot(self, distance, samples=20, fig=None, show=False, verbosity=1):

        if fig is None:
            fig = go.Figure()

        trace = mgo.cone_trace(
            self.origin,
            self.target,
            self.radius,
            distance,
            samples=samples,
            name=self.name,
            start_at_target=True,
            plot_caps=True,
            verbosity=verbosity - 1,
        )
        fig.add_trace(trace)

        if show:
            fig.show()

        return fig

    @mout.debug_log
    def always_inside(self, other):

        if isinstance(other, CappedCone):

            assert np.linalg.norm(self.origin - other.origin) == 0

            if (
                other.dist_origin_target - other.cap_radius
                < self.dist_origin_target - self.cap_radius
            ):
                print("case1")
                return False

            # if not other.is_inside(self.target):
            # 	print(f'{self}.target is not inside {other}')
            # 	return False

            # print(f'{self}.target is inside {other}')

            if other.theta < self.theta:
                print("case2")
                return False

            phi = np.arccos(np.dot(self.unit_origin_target, other.unit_origin_target))

            if other.theta < phi + self.theta:
                print("case3")
                return False

            return True

        else:
            raise Exception(
                f"Unsupported: CappedCone.always_inside(self,type({other}))"
            )

    def __repr__(self):
        if self.name:
            return self.name
        else:
            return f"CappedCone index={self.index}"
