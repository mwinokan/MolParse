def apply_rototranslation(positions, source_center, target_center, rotation_matrix):
    from numpy import dot

    return dot(positions - source_center, rotation_matrix.T) + target_center


def compute_rototranslation_matrix(source_coords, target_coords):
    """Compute the rototranslation matrix required to align one protein to another"""
    # rotation_matrix_from_points(pself.T, ptarget.T)

    raise NotImplementedError

    from numpy import dot, mean, matmul, array, double
    from ase.build.rotate import rotation_matrix_from_points

    source_center = mean(source_coords, axis=0)
    source_coords -= source_center

    target_center = mean(target_coords, axis=0)
    target_coords -= target_center

    # get the rotation matrix
    rotation = rotation_matrix_from_points(source_coords.T, target_coords.T)

    # matrix equivalent?

    # dot(positions - source_center, rotation_matrix.T) + target_center

    transformation_1 = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [*(-source_center), 1],
    ]

    transformation_2 = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [*(target_center), 1],
    ]

    transformation_1 = array(transformation_1, dtype=double)
    transformation_2 = array(transformation_2, dtype=double)
    rotation = array(rotation, dtype=double)
    # return matmul(transformation_2.T, transformation_1.T)

    # final = transformation_2 * rotation * transformation_1

    rotation = [
        [*rotation[0], 0],
        [*rotation[1], 0],
        [*rotation[2], 0],
        [0, 0, 0, 1],
    ]

    return transformation_1

    # print(rotation)
    # print(transformation_1)

    step_1 = matmul(rotation.T, transformation_1.T)
    step_2 = matmul(transformation_2.T, step_1)

    return step_2


def apply_transformation(positions, matrix):
    from numpy import matmul, mean

    positions = add_unity_column(positions)
    # print(mean(positions, axis=0))
    positions = matmul(matrix.T, positions.T)
    # print(mean(positions.T, axis=0))
    return positions[:-1, :].T
    # return dot(positions, matrix.T)


def add_unity_column(positions):
    from numpy import c_, ones  # , mean

    # print(positions)
    positions = c_[positions, ones(positions.shape[0])]
    # print(positions)
    # print(positions[:,:-1])
    return positions


# def remove_last_column(positions):
