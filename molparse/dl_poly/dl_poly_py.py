def SPE(
    control="CONTROL",
    config="CONFIG",
    field="FIELD",
    workdir="DL_POLY",
    num_procs=8,
    log="_amp.dl_poly.log",
    verbosity=2,
    output="OUTPUT",
):
    import os
    import mout
    import mcol

    if verbosity > 0:
        mout.headerOut(
            "Calling " + mcol.func + "DL_POLY" + mcol.clear + mcol.bold + "..."
        )

    if verbosity > 1:
        mout.varOut("IO: Workdir", workdir, valCol=mcol.file)
        mout.varOut("IO: Field", field, valCol=mcol.file)
        mout.varOut("IO: Config", config, valCol=mcol.file)
        mout.varOut("IO: Control", control, valCol=mcol.file)
        mout.varOut("IO: Output", output, valCol=mcol.file)

    if verbosity < 2:
        mout.redirectPrint(log)

    from dlpoly import DLPoly

    dl_poly_root = os.environ["DL_POLY"]
    dlp_exec = dl_poly_root + "/execute/DLPOLY.Z"

    DL_POLY = DLPoly(
        control=control, config=config, field=field, workdir=workdir, output=output
    )

    DL_POLY.run(executable=dlp_exec, numProcs=num_procs)

    DL_POLY.load_statis(workdir + "/STATIS")

    if verbosity < 2:
        mout.enablePrint()

    # Parse DL_POLY's OUTPUT
    num_warn = int(
        os.popen("grep 'warning - ' " + workdir + "/" + output + " | wc -l").read()
    )
    num_errs = int(
        os.popen("grep 'error - ' " + workdir + "/" + output + " | wc -l").read()
    )

    if num_warn > 0:
        mout.warningOut(
            "DL_POLY issued " + str(num_warn) + " warnings.",
            code="amp.dl_poly.dl_poly_py.SPE[1]",
        )
    if num_errs > 0:
        mout.errorOut(
            "DL_POLY issued " + str(num_errs) + " errors.",
            fatal=True,
            code="amp.dl_poly.dl_poly_py.SPE[2]",
        )

    e_tot = totalEnergy(DL_POLY, verbosity=verbosity - 1)

    return e_tot


def totalEnergy(dlpoly, verbosity=1):
    import mout
    import mcol

    if dlpoly.statis is None:
        mout.errorOut("STATIS file not loaded", code="amp.dl_poly.totalEnergy[1]")
        return None

    e_tot = dlpoly.statis.data[:, 3][0]

    if verbosity > 0:
        mout.varOut("Total Energy", e_tot, unit="kcal/mol", valCol=mcol.result)

    return e_tot
