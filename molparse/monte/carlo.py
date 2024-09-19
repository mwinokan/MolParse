from .rand import random_point_spherical
import numpy as np


def mc_spherical(
    volume,
    radius: float,
    proximity: float,
    origin: np.ndarray = None,
    samples: int = 4000,
    boolean: bool = False,
    inverse: bool = False,
    points: bool = False,
):
    """
    volume: CompoundVolume object
    radius: radius of testing sphere
    proximity: value to pass to volume.is_within
    samples: number of MC samples
    boolean: return True/False if any valid point in testing sphere
    inverse: points within the volume are invalid
    """

    num_in = 0
    num_out = 0
    points_in = []
    points_out = []

    if origin is None:
        origin = volume.origin

    for i in range(samples):
        point = random_point_spherical(origin, radius)

        if volume.is_within(point, proximity):

            # point is in volume

            if inverse:
                num_out += 1
                points_out.append(point)
            else:
                num_in += 1
                if boolean:
                    return True
                points_in.append(point)
        else:

            # point is not in volume

            if inverse:
                num_in += 1
                if boolean:
                    return True
                points_in.append(point)
            else:
                num_out += 1
                points_out.append(point)

    # fig = volume.plot(radius)

    # import mgo

    # print(f'{point=}')
    # print(f'{points_in=}')
    # print(f'{points_out=}')

    # fig.add_trace(mgo.point_trace(point))

    # fig.show()

    # exit()

    if boolean:
        return False

    sphere_volume = np.pi * np.power(radius, 3)
    fraction = num_out / samples
    valid_volume = fraction * sphere_volume

    if points:
        return fraction, valid_volume, points_in, points_out
    else:
        return fraction, valid_volume
