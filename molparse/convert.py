def pdb2traj(input, output, verbosity=1, printScript=False, tagging=True):
    import mout
    import mcol
    from .io import read, write

    if not input.endswith(".pdb"):
        mout.errorOut("Input is not a PDB!")
        return None

    if verbosity > 0:
        mout.out(
            "Converting "
            + mcol.file
            + input
            + mcol.clear
            + " to "
            + mcol.file
            + output
            + mcol.clear
            + " ... ",
            printScript=printScript,
            end="",
        )
        if verbosity > 1:
            mout.out("")

        in_traj = read(input, index=":", verbosity=verbosity - 1, tagging=tagging)

        write(output, in_traj, verbosity=verbosity - 1)

        if verbosity == 1:
            mout.out("Done.")

        return in_traj


def gro2traj(input, output, verbosity=1, printScript=False, tagging=True):
    import mout
    import mcol
    from .io import read, write

    if not input.endswith(".gro"):
        mout.errorOut("Input is not a .gro file!")
        return None

    if verbosity > 0:
        mout.out(
            "Converting "
            + mcol.file
            + input
            + mcol.clear
            + " to "
            + mcol.file
            + output
            + mcol.clear
            + " ... ",
            printScript=printScript,
            end="",
        )
        if verbosity > 1:
            mout.out("")

    if tagging:
        # initialise arrays
        taglist = []

        if verbosity > 1:
            mout.out("parsing GRO for tags ... ", printScript=printScript, end="")

        line_number = 0
        max_line_number = 1000

        # Parse the GRO to get the tags
        with open(input, "r") as input_gro:
            for line in input_gro:
                if line_number == 1:
                    max_line_number = int(line.strip().split()[0]) + 2
                if line_number == max_line_number:
                    break
                if line_number > 1:
                    taglist.append(
                        int(
                            "".join(
                                filter(lambda i: i.isdigit(), line.strip().split()[1])
                            )
                        )
                    )
                line_number = line_number + 1

        if verbosity > 1:
            mout.out("Done.")

    in_traj = read(input, index=":", verbosity=verbosity - 1)

    if tagging:
        # set the tags
        for atoms in in_traj:
            for index, tag in enumerate(taglist):
                atoms[index].tag = tag

    write(output, in_traj, verbosity=verbosity - 1)

    if verbosity == 1:
        mout.out("Done.")

    return in_traj


def xyz2traj(input, output, verbosity=1, printScript=False):
    import mout
    import mcol
    from .io import read, write

    if not input.endswith(".xyz"):
        mout.errorOut("Input is not a .xyz file!")
        return None

    if verbosity > 0:
        mout.out(
            "Converting "
            + mcol.file
            + input
            + mcol.clear
            + " to "
            + mcol.file
            + output
            + mcol.clear
            + " ... ",
            printScript=printScript,
            end="",
        )
        if verbosity > 1:
            mout.out("")

    in_traj = read(input, index=":", verbosity=verbosity - 1)

    write(output, in_traj, verbosity=verbosity - 1)

    if verbosity == 1:
        mout.out("Done.")

    return in_traj


def fileLength(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
