def plot3d(
    atoms,
    extra=[],
    bonds=[],
    alpha=1.0,
    velocity=False,
    v_scale=1.0,
    fig=None,
    flat=False,
    show=True,
    transform=None,
    title=None,
    plot_atoms=True,
    features=None,
    extra_labels=None,
):
    """Render the atoms with plotly graph objects.
    extra can contain pairs of coordinates to be shown as vectors."""

    import plotly.graph_objects as go
    from ase.data import vdw_radii, atomic_numbers
    from ase.data.colors import jmol_colors

    from .group import AtomGroup

    import mout

    group = AtomGroup.from_any("Compound AtomGroup", atoms)
    atoms = group.atoms

    species = [a.symbol for a in atoms]
    species = list(set(species))

    if fig is None:
        fig = go.Figure()

    if features:

        from .rdkit import Feature

        if isinstance(features[0], Feature):
            features = [
                dict(
                    family=f.family,
                    position=f.position,
                    indices=None,
                    x=f.x,
                    y=f.y,
                    z=f.z,
                )
                for f in features
            ]

        import plotly.express as px

        px_fig = px.scatter_3d(
            features,
            x="x",
            y="y",
            z="z",
            color="family",
            size_max=100,
            size=[100] * len(features),
            opacity=0.5,
        )

        for trace in px_fig.data:
            fig.add_trace(trace)

    if bonds:
        x, y, z = [], [], []

        for a, b in bonds:
            x.extend([a[0], b[0], None])
            y.extend([a[1], b[1], None])
            z.extend([a[2], b[2], None])

        if transform and flat:
            x, y = zip(*transform(list(zip(x, y))))

        if flat:
            trace = go.Scatter(
                x=x, y=y, mode="lines", name="bonds", line=dict(color="black", width=16)
            )
        else:
            trace = go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                name="bonds",
                line=dict(color="black", width=16),
            )

        fig.add_trace(trace)

    if not flat and velocity:

        x, y, z = [], [], []

        for atom in atoms:

            if atom.velocity is None:
                mout.warningOut(f"Atom {atom.name} (#{atom.number}) has no velocity")
                continue

            a = atom.np_pos
            b = v_scale * atom.np_vel + a

            x.extend([a[0], b[0], None])
            y.extend([a[1], b[1], None])
            z.extend([a[2], b[2], None])

        trace = go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="lines",
            name="velocity",
            line=dict(color="red", width=16),
        )

        fig.add_trace(trace)

    if plot_atoms:
        for s in species:

            atom_subset = [a for a in atoms if a.symbol == s]

            data = {}

            positions = [a.np_pos for a in atom_subset]
            x, y, z = zip(*positions)

            if transform and flat:
                x, y = zip(*transform(list(zip(x, y))))

            data["x"] = x
            data["y"] = y
            data["z"] = z

            atomic_number = atomic_numbers[s]
            size = vdw_radii[atomic_number] * 15
            color = jmol_colors[atomic_number]
            color = (color[0] * 256, color[1] * 256, color[2] * 256, alpha)

            data["index"] = [a.index for a in atom_subset]
            data["residue"] = [a.residue for a in atom_subset]
            data["res_number"] = [a.res_number for a in atom_subset]
            data["name"] = [a.name for a in atom_subset]

            customdata = [
                f"name={a.name}<br>chain={a.chain}<br>index={a.index}<br>number={a.number}<br>residue={a.residue}<br>res_index={a.res_index}<br>res_number={a.res_number}<br>x={a.x:.3f}<br>y={a.y:.3f}<br>z={a.z:.3f}"
                for a in atom_subset
            ]

            if flat:
                trace = go.Scatter(
                    x=x,
                    y=y,
                    mode="markers",
                    name=s,
                    marker=dict(
                        size=size,
                        color=f"rgba{color}",
                        line=dict(color="black", width=2),
                    ),
                    customdata=customdata,
                    hovertemplate="%{customdata}<extra></extra>",
                )
            else:
                trace = go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode="markers",
                    name=s,
                    marker=dict(
                        size=size,
                        color=f"rgba{color}",
                        line=dict(color="black", width=2),
                    ),
                    customdata=customdata,
                    hovertemplate="%{customdata}<extra></extra>",
                )

            fig.add_trace(trace)

    if not flat:

        if extra_labels:

            for (a, b), l in zip(extra, extra_labels):
                trace = go.Scatter3d(
                    x=[a[0], b[0]], y=[a[1], b[1]], z=[a[2], b[2]], name=l
                )

                fig.add_trace(trace)

        else:
            for i, (a, b) in enumerate(extra):
                trace = go.Scatter3d(
                    x=[a[0], b[0]], y=[a[1], b[1]], z=[a[2], b[2]], name=f"extra[{i}]"
                )

                fig.add_trace(trace)

    fig.update_layout(scene_aspectmode="data")

    if title:
        fig.update_layout(title=title)

    if show:
        fig.show()

    return fig
