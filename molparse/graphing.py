"""

  To-Do's
    * If show = True close all other figures?
    * graphVelocity
    * graphTemperature
    * graphGyration

"""


def graphEnergy(
    trajectory,
    perAtom=True,
    filename=None,
    show=True,
    verbosity=2,
    kJpermol=False,
    xlab=None,
    timestep=None,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
):
    import mcol
    import mout
    import mplot
    from ase import units

    if verbosity > 0:
        mout.out(
            "graphing " + mcol.varName + "Energy" + mcol.clear + " ... ",
            printScript=True,
            end="",
        )  # user output

    xdata = []
    ekins = []
    epots = []
    etots = []

    if not kJpermol:
        if perAtom:
            ylab = "Energy [eV/atom]"
        else:
            ylab = "Energy [eV]"
    else:
        ylab = "Energy [kJ/mol]"

    for n, atoms in enumerate(trajectory):
        if timestep is None:
            if xlab is None:
                xlab = "MD Steps"
            xdata.append(n)
        else:
            xlab = "Time [fs]"
            xdata.append(n * timestep / units.fs)
        epot = atoms.get_potential_energy()
        ekin = atoms.get_kinetic_energy()

        if not kJpermol:
            if perAtom:
                epot /= len(atoms)
                ekin /= len(atoms)
                ylab = "Energy [eV/atom]"
            else:
                ylab = "Energy [eV]"
        else:
            epot *= units.eV / units.kJ * units.mol
            ekin *= units.eV / units.kJ * units.mol
            ylab = "Energy [kJ/mol]"

        epots.append(epot)
        ekins.append(ekin)
        etots.append(epot + ekin)

    mplot.graph2D(
        xdata,
        [epots, ekins, etots],
        ytitles=["Potential", "Kinetic", "Total"],
        show=show,
        xlab=xlab,
        ylab=ylab,
        filename=filename,
        verbosity=verbosity - 1,
        xmin=xmin,
        xmax=xmax,
        ymin=ymin,
        ymax=ymax,
    )

    if verbosity > 0:
        mout.out("Done.")  # user output


def graphForces(
    trajectory,
    filename=None,
    max=True,
    show=True,
    verbosity=2,
    xlab="Step",
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
    yLog=False,
):
    import mcol
    import mout
    import mplot
    import numpy as np

    if verbosity > 0:
        mout.out(
            "graphing " + mcol.varName + "Forces" + mcol.clear + " ... ",
            printScript=True,
            end="",
        )  # user output

    xdata = []
    fmax = []
    favg = []

    if max:
        ylab = "Force"
    else:
        ylab = "Average Force"

    for n, atoms in enumerate(trajectory):
        xdata.append(n)

        forces = atoms.get_forces()

        this_fmax = np.max(forces)
        this_favg = np.average(forces)

        fmax.append(this_fmax)
        favg.append(this_favg)

    if max:
        ydata = [fmax, favg]
        ytitles = ["Maximum", "Average"]
    else:
        ydata = favg
        ytitles = "Average"

    mplot.graph2D(
        xdata,
        ydata,
        ytitles=ytitles,
        show=show,
        xlab=xlab,
        ylab=ylab,
        filename=filename,
        verbosity=verbosity - 1,
        xmin=xmin,
        xmax=xmax,
        ymin=ymin,
        ymax=ymax,
        ySci=True,
        yLog=yLog,
    )

    if verbosity > 0:
        mout.out("Done.")  # user output


def graphDisplacement(
    trajectory, show=True, filename=None, relative=True, verbosity=2, timestep=None
):
    """
    Root Mean Square Displacement

    """

    import mcol
    import mout
    import mplot
    import numpy as np
    from ase import units

    if verbosity > 0:
        mout.out(
            "graphing " + mcol.varName + "RMSD" + mcol.clear + " ... ",
            printScript=True,
            end="",
        )  # user output

    xdata = []
    rmsd = []

    for n, atoms in enumerate(trajectory):

        if timestep is None:
            xlab = "MD Steps"
            xdata.append(n)
        else:
            xlab = "Time [fs]"
            xdata.append(n * timestep / units.fs)

        positions = atoms.get_positions()

        if relative:
            if n == 0:
                reference = positions.copy()
            positions -= reference

        rmsd.append(np.sqrt(np.mean(positions**2)))

    mplot.graph2D(
        xdata,
        rmsd,
        show=show,
        xlab="MD Steps",
        ylab="RMS Displacement",
        filename=filename,
        verbosity=verbosity - 1,
    )

    if verbosity > 0:
        mout.out("Done.")  # user output


def graphBondLength(
    trajectory,
    indices,
    torsion_indices=None,
    printScript=False,
    show=True,
    filename=None,
    fitMin=None,
    fitMax=None,
    verbosity=2,
    timestep=None,
    title=None,
    fitOrder=None,
    yUnit="Angstroms",
    dataFile=None,
    ymin=0,
    ymax=None,
    xmin=None,
    xmax=None,
    noplot=False,
    write_data=False,
    write_ang=False,
    write_torsion=False,
    write_exang=False,
    debug_plot=False,
    return_data=False,
):
    from .analysis import bondLengthStats

    """
      Graph the bond lengths (displacement) between atoms.

    """

    # import mcol
    import mout
    import mplot
    import numpy as np

    many = any(isinstance(el, list) for el in indices)

    ydata = []
    labels = []

    if timestep is None:
        xlab = "MD Steps"
        xUnit = "MD Step"
    else:
        xlab = "Time [fs]"
        xUnit = "femtoseconds"

    if many:

        for i, pair in enumerate(indices):

            if verbosity > 2:
                print("")
            val, err, label, xdata, this_data = bondLengthStats(
                trajectory,
                pair,
                printScript=printScript,
                verbosity=verbosity - 2,
                timestep=timestep,
                yUnit=yUnit,
                returnData=True,
            )

            labels.append(label)
            ydata.append(this_data)

        if fitOrder is None:
            if not noplot:
                mplot.graph2D(
                    xdata,
                    ydata,
                    ytitles=labels,
                    show=show,
                    xlab=xlab,
                    ylab="Distance [Angstrom]",
                    filename=filename,
                    title=title,
                    verbosity=verbosity - 1,
                    ymin=ymin,
                    ymax=ymax,
                    xmin=xmin,
                    xmax=xmax,
                )
        else:
            if verbosity > 1:
                print("")
            if len(xdata) < fitOrder + 1:
                mout.warningOut("Not enough data-points for fitting.")
                fitOrder = None
            if title is not None:
                dataFile.write("# ")
                dataFile.write(title)
                dataFile.write("\n")
            val, err, fit_func = mplot.fit(
                xdata,
                ydata,
                rank=fitOrder,
                verbosity=verbosity - 1,
                title=title,
                fitMin=fitMin,
                fitMax=fitMax,
                yUnit=yUnit,
                xUnit=xUnit,
                dataFile=dataFile,
            )
            text = mplot.getCoeffStr(val, err, 1, yUnit=yUnit, xUnit=xUnit)
            if not noplot:
                mplot.graph2D(
                    xdata,
                    ydata,
                    fitFunc=fit_func,
                    ytitles=labels,
                    show=show,
                    xlab=xlab,
                    ylab="Distance [Angstrom]",
                    filename=filename,
                    title=title,
                    verbosity=verbosity,
                    subtitle=text,
                    ymin=ymin,
                    ymax=ymax,
                    xmin=xmin,
                    xmax=xmax,
                )

    else:

        if verbosity > 2:
            print("")
        val, err, label, xdata, ydata = bondLengthStats(
            trajectory,
            indices,
            printScript=printScript,
            verbosity=verbosity - 2,
            timestep=timestep,
            yUnit=yUnit,
            returnData=True,
        )

        if fitOrder is None:
            if not noplot:
                mplot.graph2D(
                    xdata,
                    ydata,
                    show=show,
                    xlab=xlab,
                    ylab="Distance [Angstrom]",
                    filename=filename,
                    title=label,
                    verbosity=verbosity - 1,
                    ymin=ymin,
                    ymax=ymax,
                    xmin=xmin,
                    xmax=xmax,
                )
        else:
            if verbosity > 1:
                print("")
            if len(xdata) < fitOrder + 1:
                mout.warningOut("Not enough data-points for fitting.")
                fitOrder = None
            if title is not None:
                dataFile.write("#")
                dataFile.write(title)
                dataFile.write("\n")
            val, err, fit_func = mplot.fit(
                xdata,
                ydata,
                rank=fitOrder,
                verbosity=verbosity - 1,
                fitMin=fitMin,
                fitMax=fitMax,
                title=title,
                yUnit=yUnit,
                xUnit=xUnit,
                dataFile=dataFile,
            )
            text = mplot.getCoeffStr(val, err, 1, yUnit=yUnit, xUnit=xUnit)
            if not noplot:
                mplot.graph2D(
                    xdata,
                    ydata,
                    fitFunc=fit_func,
                    show=show,
                    xlab=xlab,
                    ylab="Distance [Angstrom]",
                    filename=filename,
                    title=label,
                    verbosity=verbosity,
                    subtitle=text,
                    ymin=ymin,
                    ymax=ymax,
                    xmin=xmin,
                    xmax=xmax,
                )

    if write_data:

        import os

        # base=os.path.basename(filename)
        # data_dump = open(os.path.splitext(base)[0]+".dat",'w')
        data_dump = open(filename.replace(".png", ".dat"), "w")

        if many:

            data_dump.write("# x " + str(labels))

            if write_torsion:
                data_dump.write(" opening_angle opening_torsion")
                if write_exang:
                    data_dump.write(" psi phi")
            elif write_ang:
                data_dump.write(" opening_angle")

            data_dump.write("\n")
            data_dump.close()

            # generate the torsion data
            if write_torsion:
                angdata = []
                torsdata = []
                exangdata1 = []
                exangdata2 = []
                for image in trajectory:
                    positions = image.get_positions()
                    symbols = image.get_chemical_symbols()

                    # DC:O2 -> DC:N4
                    # print(1,symbols[torsion_indices[1]],symbols[torsion_indices[0]])
                    vec1 = positions[torsion_indices[1]] - positions[torsion_indices[0]]

                    # DG:N2 -> DG:O6
                    # print(2,symbols[torsion_indices[5]],symbols[torsion_indices[4]])
                    vec2 = positions[torsion_indices[5]] - positions[torsion_indices[4]]

                    # DC:N1 -> DC:N3
                    # print(3,symbols[torsion_indices[3]],symbols[torsion_indices[2]])
                    vec1p = (
                        positions[torsion_indices[3]] - positions[torsion_indices[2]]
                    )

                    # DG:C4 -> DG:N1
                    # print(4,symbols[torsion_indices[7]],symbols[torsion_indices[6]])
                    vec2p = (
                        positions[torsion_indices[7]] - positions[torsion_indices[6]]
                    )

                    vec1pp = np.cross(vec1, vec1p)
                    vec2pp = np.cross(vec2, vec2p)

                    unit_vec1 = vec1 / np.linalg.norm(vec1)
                    unit_vec2 = vec2 / np.linalg.norm(vec2)
                    # unit_vec1p = vec1p / np.linalg.norm(vec1p)
                    # unit_vec2p = vec2p / np.linalg.norm(vec2p)
                    unit_vec1pp = vec1pp / np.linalg.norm(vec1pp)
                    unit_vec2pp = vec2pp / np.linalg.norm(vec2pp)

                    unit_vec1pp_dot_vec2pp = np.dot(unit_vec1pp, unit_vec2pp)
                    tau = np.arccos(unit_vec1pp_dot_vec2pp)
                    torsdata.append(np.pi - tau)

                    unit_vec1_dot_vec2 = np.dot(unit_vec1, unit_vec2)
                    theta = np.arccos(unit_vec1_dot_vec2)

                    # we potentially need to flip the angles
                    vec1_cross_vec2 = np.cross(vec1, vec2)
                    flip = np.dot(vec1_cross_vec2, vec2pp) > 0
                    if flip:
                        theta = -theta
                    angdata.append(theta)

                    if write_exang:
                        # calculate the projected opening angles
                        vec1_dot_vec2 = np.dot(vec1, vec2)
                        vec1ppp = vec1p - np.dot(vec1p, vec1) * vec1
                        vec2ppp = vec2p - np.dot(vec2p, vec2) * vec2
                        psi = np.arctan(np.dot(vec2, vec1ppp) / vec1_dot_vec2)
                        phi = np.arctan(np.dot(vec1, vec2ppp) / vec1_dot_vec2)

                        # Because it is based on atan, it does not need flipping
                        exangdata1.append(psi)
                        exangdata2.append(phi)

                if debug_plot:
                    import matplotlib.pyplot as plt

                    fig, ax = plt.subplots()
                    plt.plot(xdata, ydata[0], label="bot")
                    # plt.plot(xdata,ydata[1],label="mid")
                    plt.plot(xdata, ydata[2], label="top")
                    plt.plot(xdata, angdata, label="ang")
                    plt.plot(xdata, exangdata1, label="psi")
                    plt.plot(xdata, exangdata2, label="phi")
                    plt.plot(xdata, torsdata, label="tors")
                    ax.set_xlim(0, xdata[100])
                    plt.legend(loc="best")
                    plt.show()
                    plt.close()
                    exit()

            # generate the angle data
            elif write_ang:
                angdata = []
                for image in trajectory:
                    positions = image.get_positions()
                    vec1 = positions[indices[0][-1]] - positions[indices[-1][0]]
                    vec2 = positions[indices[0][0]] - positions[indices[-1][-1]]
                    unit_vec1 = vec1 / np.linalg.norm(vec1)
                    unit_vec2 = vec2 / np.linalg.norm(vec2)
                    dot_product = np.dot(unit_vec1, unit_vec2)
                    angle = np.arccos(dot_product)
                    if angle > np.pi / 2:
                        angle = angle - np.pi
                    angdata.append(angle)

            # data_dump = open(os.path.splitext(base)[0]+".dat",'a')
            data_dump = open(filename.replace(".png", ".dat"), "a")

            for index, x in enumerate(xdata):
                # data_dump.write(str(xdata[index])+" ")
                data_dump.write(f"{xdata[index]:.5f} ")
                for data in ydata:
                    # data_dump.write(str(data[index])+" ")
                    data_dump.write(f"{data[index]:.5f} ")
                if write_torsion:
                    # data_dump.write(str(angdata[index])+" ")
                    data_dump.write(f"{angdata[index]:.5f} ")
                    data_dump.write(f"{torsdata[index]:.5f} ")
                    if write_exang:
                        data_dump.write(f"{exangdata1[index]:.5f} ")
                        data_dump.write(f"{exangdata2[index]:.5f} ")
                elif write_ang:
                    # data_dump.write(str(angdata[index])+" ")
                    data_dump.write(f"{angdata[index]:.5f} ")
                data_dump.write("\n")

        else:

            data_dump.write("# x " + label + "\n")
            data_dump.close()

            # data_dump = open(os.path.splitext(base)[0]+".dat",'a')
            data_dump = open(filename.replace(".png", ".dat"), "a")

            for index, x in enumerate(xdata):
                data_dump.write(str(xdata[index]) + " ")
                data_dump.write(str(ydata[index]) + " ")
                data_dump.write("\n")

        data_dump.close()

        mout.out("Data dumped to " + filename.replace(".png", ".dat"))

    if fitOrder is not None:
        return val, err, fit_func
    elif return_data:
        if torsion_indices is not None:
            return xdata, ydata, angdata, torsdata
        return xdata, ydata
    else:
        return None, None, None


def graphBondVibSpec(
    trajectory,
    indices,
    printScript=False,
    show=True,
    filename=None,
    verbosity=2,
    timestep=None,
    title=None,
    ymin=None,
    ymax=None,
    power_spec=False,
    power_segment=None,
    wave_number=False,
    xlab="Frequency [Hz]",
    write_data=False,
):
    # import mcol
    import mout
    import mplot
    from .analysis import getBondLabel

    """
      Graph the fourier transfort of bond lengths (displacement) between atoms.

    """

    # if (verbosity > 0):
    #   mout.out("graphing "+mcol.varName+
    #            title+
    #            mcol.clear+" ... ",
    #            printScript=printScript,
    #            end='') # user output

    many = any(isinstance(el, list) for el in indices)

    if many:
        mout.errorOut("Unsupported.", fatal=True)

    if not many:

        xdata = []
        ydata = []

        bond_label = getBondLabel(trajectory, indices)

        index1 = indices[0]
        index2 = indices[1]

        for n, atoms in enumerate(trajectory):
            # if timestep is None:
            xdata.append(n)

            dist = atoms.get_distance(index1, index2)
            ydata.append(dist)

        # print(len(xdata))
        # print(len(ydata))

        # import scipy
        import scipy.fftpack

        Ydata = scipy.fftpack.fft(ydata)[: len(ydata) // 2]

        if timestep is not None:
            Xdata = scipy.fftpack.fftfreq(n, timestep)[: len(ydata) // 2]
        else:
            Xdata = xdata[: len(ydata) // 2]

        if power_spec:
            if timestep is None:
                mout.errorOut("Timestep required for power spectrum", fatal=True)
            import scipy.signal

            if power_segment is None:
                power_segment = len(ydata) // 10
            f, Pxx_den = scipy.signal.welch(
                ydata, fs=1 / timestep, nperseg=power_segment
            )

            Xdata = f
            Ydata = Pxx_den

        if wave_number:
            temp_x = []
            light_speed = 299792458
            for x in Xdata:
                temp_x.append(x / light_speed / 100.0)
            Xdata = temp_x
            xlab = "Wave Number [/cm]"

        # ydata = np.fft.fft(ydata)

        # if ymin is None:
        #   ymin=min(ydata)
        # if ymax is None:
        #   ymax=max(ydata)

        mplot.graph2D(
            Xdata,
            Ydata,
            ytitles=bond_label,
            show=show,
            filename=filename,
            verbosity=verbosity - 1,
            yLog=True,
            ymin=ymin,
            ymax=ymax,
            xlab=xlab,
        )

        if write_data:

            import os

            base = os.path.basename(filename)
            data_dump = open(os.path.splitext(base)[0] + ".dat", "w")

            if many:

                # NOT UPDATED

                data_dump.write("# x " + str(labels) + "\n")
                data_dump.close()

                data_dump = open(os.path.splitext(base)[0] + ".dat", "a")

                for index, x in enumerate(xdata):
                    data_dump.write(str(xdata[index]) + " ")
                    for data in ydata:
                        data_dump.write(str(data[index]) + " ")
                    data_dump.write("\n")

            else:

                # bond_label = getBondLabel(trajectory,indices)

                if wave_number:
                    data_dump.write("# wave_number [/cm], " + bond_label + "\n")
                else:
                    data_dump.write("# frequency [Hz], " + bond_label + "\n")
                data_dump.close()

                data_dump = open(os.path.splitext(base)[0] + ".dat", "a")

                for index, x in enumerate(Xdata):
                    data_dump.write(str(Xdata[index]) + " ")
                    data_dump.write(str(Ydata[index]) + " ")
                    data_dump.write("\n")

            data_dump.close()


# just a wrapper for mplot.show()
def showFigs(verbosity=1):
    mplot.show(verbosity=verbosity)


def graphBondAngle(
    trajectory,
    indices,
    printScript=False,
    show=True,
    filename=None,
    fitMin=None,
    fitMax=None,
    verbosity=2,
    timestep=None,
    title=None,
    fitOrder=None,
    yUnit="Degrees",
    dataFile=None,
    ymin=0,
    ymax=180,
    xmin=None,
    xmax=None,
):
    import mout
    import mplot
    from .analysis import bondAngleStats

    """
      Graph the bond angle (acute) between atoms.

    """

    many = any(isinstance(el, list) for el in indices)

    ydata = []
    labels = []

    if timestep is None:
        xlab = "MD Steps"
        xUnit = "MD Step"
    else:
        xlab = "Time [fs]"
        xUnit = "picosecond"

    if many:

        for i, triple in enumerate(indices):

            if verbosity > 2:
                print("")
            val, err, label, xdata, this_data = bondAngleStats(
                trajectory,
                triple,
                printScript=printScript,
                verbosity=verbosity - 2,
                timestep=timestep,
                yUnit=yUnit,
                returnData=True,
            )

            labels.append(label)
            ydata.append(this_data)

        if fitOrder is None:
            mplot.graph2D(
                xdata,
                ydata,
                ytitles=labels,
                show=show,
                xlab=xlab,
                ylab="Distance [Angstrom]",
                filename=filename,
                title=title,
                verbosity=verbosity - 1,
            )
        else:
            if verbosity > 1:
                print("")
            if len(xdata) < fitOrder + 1:
                mout.warningOut("Not enough data-points for fitting.")
                fitOrder = None
            if title is not None:
                dataFile.write("#")
                dataFile.write(title)
                dataFile.write("\n")
            val, err, fit_func = mplot.fit(
                xdata,
                ydata,
                rank=fitOrder,
                verbosity=verbosity - 1,
                title=title,
                fitMin=fitMin,
                fitMax=fitMax,
                yUnit=yUnit,
                xUnit=xUnit,
                dataFile=dataFile,
            )
            text = mplot.getCoeffStr(val, err, 1, yUnit=yUnit, xUnit=xUnit)
            mplot.graph2D(
                xdata,
                ydata,
                fitFunc=fit_func,
                ytitles=labels,
                show=show,
                xlab=xlab,
                ylab="Angle [Degrees]",
                filename=filename,
                title=title,
                verbosity=verbosity,
                subtitle=text,
                ymin=ymin,
                ymax=ymax,
                xmin=xmin,
                xmax=xmax,
            )

    else:

        if verbosity > 2:
            print("")
        val, err, label, xdata, ydata = bondAngleStats(
            trajectory,
            indices,
            printScript=printScript,
            verbosity=verbosity - 2,
            timestep=timestep,
            yUnit=yUnit,
            returnData=True,
        )

        if fitOrder is None:
            mplot.graph2D(
                xdata,
                ydata,
                show=show,
                xlab=xlab,
                ylab="Distance [Angstrom]",
                filename=filename,
                title=label,
                verbosity=verbosity - 1,
            )
        else:
            if verbosity > 1:
                print("")
            if len(xdata) < fitOrder + 1:
                mout.warningOut("Not enough data-points for fitting.")
                fitOrder = None
            if title is not None:
                dataFile.write("# ")
                dataFile.write(title)
                dataFile.write("\n")
            val, err, fit_func = mplot.fit(
                xdata,
                ydata,
                rank=fitOrder,
                verbosity=verbosity - 1,
                fitMin=fitMin,
                fitMax=fitMax,
                title=title,
                yUnit=yUnit,
                xUnit=xUnit,
                dataFile=dataFile,
            )
            text = mplot.getCoeffStr(val, err, 1, yUnit=yUnit, xUnit=xUnit)
            mplot.graph2D(
                xdata,
                ydata,
                fitFunc=fit_func,
                show=show,
                xlab=xlab,
                ylab="Distance [Angstrom]",
                filename=filename,
                title=label,
                verbosity=verbosity,
                subtitle=text,
                ymin=ymin,
                ymax=ymax,
                xmin=xmin,
                xmax=xmax,
            )

    if fitOrder is not None:
        return val, err, fit_func
    else:
        return None, None, None
