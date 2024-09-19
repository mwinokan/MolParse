# to-do:

"""

	improve performance by creating a smaller subsystem just of varying atoms, interpolating those and writing them to disk. Then use SED or other bash tools to insert the modifications into the original pdb

"""


def interpolate(
    input1,
    input2,
    frames=10,
    verbosity=2,
    smooth=False,
    indices=None,
    frame_padding=0,
    grid=False,
    ratios=None,
):
    import mout
    from ase import Atoms as aseatoms
    from .system import System

    if verbosity > 0:
        mout.headerOut("Interpolating...")
    if verbosity > 1:
        mout.varOut("input1", str(input1))
        mout.varOut("input2", str(input2))
        mout.varOut("frames", frames)
        mout.varOut("frame_padding", frame_padding)

    if isinstance(input1, aseatoms):
        mout.out("Input1 = ase.Atoms")
        mout.errorOut("Not supported yet", fatal=True)
    elif isinstance(input1, System):
        mout.out("Input1 = amp.system")
    else:
        mout.out(str(type(input1)))
        mout.errorOut("Not supported yet", fatal=True)

    if isinstance(input2, aseatoms):
        mout.out("Input2 = ase.Atoms")
        mout.errorOut("Not supported yet", fatal=True)
    elif isinstance(input2, System):
        mout.out("Input2 = amp.system")
    else:
        mout.out(str(type(input2)))
        mout.errorOut("Not supported yet", fatal=True)

    # these copy calls take a very long time
    system = input1.copy(alt=False)

    if not grid:

        system_array = []
        iterange = range(-frame_padding, frames + frame_padding)

        print(iterange)

        if indices is None and ratios is not None:
            mout.warningOut("Ratios not used as no indices passed!")
        elif indices is not None:
            assert len(ratios) == len(indices)

        # mout.progress(0,frames+2*frame_padding)

        for i in iterange:

            system.name = "Interpolation Frame " + str(i)

            if indices is not None:
                for j, index in enumerate(indices):
                    atom = system.atoms[index]

                    if ratios is not None:
                        atom.position = custom_interpolate(
                            input1.atoms[index].np_pos,
                            input2.atoms[index].np_pos,
                            ratios[j][i],
                        )
                    elif smooth:
                        atom.position = smooth_interpolate(
                            input1.atoms[index].np_pos,
                            input2.atoms[index].np_pos,
                            frames,
                            i,
                        )
                    else:
                        atom.position = simple_interpolate(
                            input1.atoms[index].np_pos,
                            input2.atoms[index].np_pos,
                            frames,
                            i,
                        )

            else:

                for index, atom in enumerate(system.atoms):
                    if smooth:
                        atom.position = smooth_interpolate(
                            input1.atoms[index].np_pos,
                            input2.atoms[index].np_pos,
                            frames,
                            i,
                        )
                    else:
                        atom.position = simple_interpolate(
                            input1.atoms[index].np_pos,
                            input2.atoms[index].np_pos,
                            frames,
                            i,
                        )

            # these copy calls take up a very long time
            system_array.append(system.copy(alt=False))
            mout.progress(i + frame_padding, frames + 2 * frame_padding)

        mout.progress(frames + 2 * frame_padding, frames + 2 * frame_padding)
        return system_array

    else:

        assert indices is not None
        assert len(indices) == 2
        assert frames <= 10
        assert not smooth
        assert ratios is None

        """

        AA AB AC
        BA BB BC
        CA CB CC

        """

        import numpy as np

        system_array = np.empty(
            (2 * frame_padding + frames, 2 * frame_padding + frames), dtype=object
        )
        # names_array = np.empty((2*frame_padding+frames,2*frame_padding+frames),dtype=str)

        iterange = range(-frame_padding, frames + frame_padding)

        system_array[frame_padding][frame_padding] = input1.copy(alt=False)
        system_array[-1 - frame_padding][-1 - frame_padding] = input2.copy(alt=False)

        corner1 = input1.copy(alt=False)
        corner2 = input1.copy(alt=False)

        corner1.atoms[indices[0]].position = input2.atoms[indices[0]].position
        corner2.atoms[indices[1]].position = input2.atoms[indices[1]].position

        # naming scheme
        from string import ascii_uppercase

        names = []
        for i in range(frame_padding):
            names.append(ascii_uppercase[-frame_padding + i])
        for i in range(frames):
            names.append(f"{i}")
        for i in range(frame_padding):
            names.append(ascii_uppercase[i])

        for i in iterange:
            for j in iterange:

                if system_array[i][j] is None:
                    system = input1.copy()

                    atom1 = system.atoms[indices[0]]
                    atom1.position = simple_interpolate(
                        input1.atoms[indices[0]].np_pos,
                        corner1.atoms[indices[0]].np_pos,
                        frames,
                        i,
                    )
                    atom2 = system.atoms[indices[1]]
                    atom2.position = simple_interpolate(
                        input1.atoms[indices[1]].np_pos,
                        corner2.atoms[indices[1]].np_pos,
                        frames,
                        j,
                    )

                    system_array[i][j] = system

                system_array[i][j].name = f"{names[i]}{names[j]}"

        return system_array


def custom_interpolate(start, end, ratio):
    return start + ratio * (end - start)


def simple_interpolate(start, end, frames, i):
    return start + i * (end - start) / (frames - 1)


def smooth_interpolate(start, end, frames, i):
    import math

    angle = simple_interpolate(0.0, math.pi, frames, i)
    return start + 0.5 * (1 - math.cos(angle)) * (end - start)


def auto_rotate(atoms):
    from .system import System

    if isinstance(atoms, System):
        sys = atoms
        atoms = auto_rotate(sys.ase_atoms)
        sys.set_coordinates(atoms)
        return sys

    atoms = atoms.copy()

    positions = atoms.get_positions()

    positions = [
        [p[0] for p in positions],
        [p[1] for p in positions],
        [p[2] for p in positions],
    ]

    import numpy as np
    from numpy.linalg import svd

    # fit a plane to the atomic positions
    points = np.reshape(
        positions, (np.shape(positions)[0], -1)
    )  # Collapse trialing dimensions
    assert (
        points.shape[0] <= points.shape[1]
    ), "There are only {} points in {} dimensions.".format(
        points.shape[1], points.shape[0]
    )
    central_point = points.mean(axis=1)
    x = points - central_point[:, np.newaxis]
    M = np.dot(x, x.T)  # Could also use np.cov(x) here.
    normal_vector = svd(M)[0][:, -1]

    # rotate the normal onto the Z-axis
    atoms.rotate(normal_vector, "z")

    # rotate the position vector of the first atom onto the x-axis
    positions = atoms.get_positions()
    index0_vector = np.array(
        [np.dot(positions[0], [1, 0, 0]), np.dot(positions[0], [0, 1, 0]), 0]
    )
    index1_vector = np.array(
        [np.dot(positions[1], [1, 0, 0]), np.dot(positions[1], [0, 1, 0]), 0]
    )
    # index0_vector = positions[0]
    # index1_vector = positions[1]

    if index0_vector[0] < 0:
        index0_vector = [index0_vector[0], index0_vector[1], 0]
    # print(index0_vector)
    atoms.rotate(index1_vector - index0_vector, "x")
    # print(atoms.get_positions()[0])

    return atoms
