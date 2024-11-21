# implementation of SuCOS from https://github.com/susanhleung/SuCOS https://doi.org/10.26434/chemrxiv.8100203.v1

import os
from numpy import clip

from rdkit import RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem.rdmolops import CombineMols
import mrich
from pandas import DataFrame

# feature setup

feature_factory = AllChem.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
)

feature_map_params = {}
for k in feature_factory.GetFeatureFamilies():
    feature_params = FeatMaps.FeatMapParams()
    feature_map_params[k] = feature_params

keep = (
    "Donor",
    "Acceptor",
    "NegIonizable",
    "PosIonizable",
    "ZnBinder",
    "Aromatic",
    "Hydrophobe",
    "LumpedHydrophobe",
)


def feature_map_score(
    inspiration,
    derivative,
    score_mode=FeatMaps.FeatMapScoreMode.All,
    debug=False,
    draw=False,
    return_data: bool = False,
):

    inspiration_features = [
        f
        for f in feature_factory.GetFeaturesForMol(inspiration)
        if f.GetFamily() in keep
    ]
    derivative_features = [
        f
        for f in feature_factory.GetFeaturesForMol(derivative)
        if f.GetFamily() in keep
    ]

    if draw:
        from .draw import draw_mol

        mrich.h3("inspiration")
        draw_mol(inspiration, feats=inspiration_features)
        mrich.h3("derivative")
        draw_mol(derivative, feats=derivative_features)

    inspiration_feature_map = FeatMaps.FeatMap(
        feats=inspiration_features,
        weights=[1] * len(inspiration_features),
        params=feature_map_params,
    )
    derivative_feature_map = FeatMaps.FeatMap(
        feats=derivative_features,
        weights=[1] * len(derivative_features),
        params=feature_map_params,
    )

    inspiration_feature_map.scoreMode = score_mode

    feature_score = inspiration_feature_map.ScoreFeats(derivative_features)

    score_vector = [0] * inspiration_feature_map.GetNumFeatures()

    inspiration_feature_map.ScoreFeats(derivative_features, mapScoreVect=score_vector)

    if return_data or debug:
        inspiration_data = []
        for feature, score in zip(inspiration_features, score_vector):
            inspiration_data.append(
                dict(
                    family=feature.GetFamily(),
                    type=feature.GetType(),
                    atom_ids=str(feature.GetAtomIds()),
                    score=score,
                )
            )

    if debug:

        from rich.table import Table

        mrich.h3("Debug")

        if draw:
            from .draw import draw_flat

            display(draw_flat(inspiration, indices=True))
            display(draw_flat(derivative, indices=True))

        table = Table(title="inspiration features")
        table.add_column("family")
        table.add_column("type")
        table.add_column("atom_ids")
        table.add_column("recapitulation score")

        for d in inspiration_data:
            table.add_row(
                d["family"], d["type"], str(d["atom_ids"]), f"{d['score']:.3f}"
            )
        mrich.print(table)

        table = Table(title="derivative features")
        table.add_column("family")
        table.add_column("type")
        table.add_column("atom_ids")
        # table.add_column("recapitulation score")

        for feature in derivative_features:
            # mrich.debug(feature.GetFamily(), feature.GetType(), feature.GetAtomIds(), score)
            table.add_row(
                feature.GetFamily(), feature.GetType(), str(feature.GetAtomIds())
            )
        mrich.print(table)

        # mrich.print(dir(feature.GetPos()))

    feature_score /= min(
        inspiration_feature_map.GetNumFeatures(), len(derivative_features)
    )

    if return_data:
        return feature_score, inspiration_data

    return feature_score


def multi_feature_map_score(
    inspirations,
    derivative,
    draw=False,
    debug=False,
    score_mode=FeatMaps.FeatMapScoreMode.All,
):

    inspiration_feature_lists = [
        [f for f in feature_factory.GetFeaturesForMol(x) if f.GetFamily() in keep]
        for x in inspirations
    ]
    derivative_features = [
        f
        for f in feature_factory.GetFeaturesForMol(derivative)
        if f.GetFamily() in keep
    ]

    combined_inspiration_features = sum(inspiration_feature_lists, [])

    if draw:
        from .draw import draw_mol

        for inspiration, features in zip(inspirations, inspiration_feature_lists):
            draw_mol(inspiration, feats=features)

        draw_mol(derivative, feats=combined_inspiration_features)

        draw_mol(derivative, feats=derivative_features)

    derivative_feature_map = FeatMaps.FeatMap(
        feats=derivative_features,
        weights=[1] * len(derivative_features),
        params=feature_map_params,
    )

    # loop over inspirations

    score_dict = {}

    for i, inspiration in enumerate(inspirations):

        score_dict[i] = {}

        if debug:
            mrich.title(inspiration)

        inspiration_features = [
            f
            for f in feature_factory.GetFeaturesForMol(inspiration)
            if f.GetFamily() in keep
        ]
        inspiration_feature_map = FeatMaps.FeatMap(
            feats=inspiration_features,
            weights=[1] * len(inspiration_features),
            params=feature_map_params,
        )

        score_vector = [0] * derivative_feature_map.GetNumFeatures()

        feature_score = derivative_feature_map.ScoreFeats(
            inspiration_features, mapScoreVect=score_vector
        )

        for j, (feature, score) in enumerate(zip(derivative_features, score_vector)):
            if debug:
                mrich.var(feature.GetFamily(), score)

            score_dict[i][j] = score

    i_list = list(score_dict.keys())

    combined_score_vector = []

    for j, feature in enumerate(derivative_features):

        scores = [score_dict[i][j] for i in i_list]

        best_score = max(scores)

        combined_score_vector.append(best_score)

        # print(feature.GetFamily(), scores, best_score)

    feature_score = sum(combined_score_vector) / len(derivative_features)

    return feature_score


def SuCOS_score(
    inspiration,
    derivative,
    recapitulation_threshold: float = 0.6,
    return_all=False,
    print_scores=False,
    debug: bool = False,
    draw: bool = False,
    **kwargs,
):

    if isinstance(inspiration, list):
        if len(inspiration) == 1:
            multi = False
            inspiration = inspiration[0]
        else:
            multi = True
    else:
        multi = False

    if multi:
        feature_score = multi_feature_map_score(
            inspiration, derivative, debug=debug, draw=draw, **kwargs
        )
        inspiration_data = []
    else:
        feature_score, inspiration_data = feature_map_score(
            inspiration, derivative, debug=debug, draw=draw, return_data=True, **kwargs
        )

    if debug or print_scores or return_all:
        df = DataFrame(inspiration_data)
        df = df.groupby("atom_ids").max("score")
        recapitulation_count = len(df[df["score"] > recapitulation_threshold])
        recapitulation_fraction = recapitulation_count / len(df)

    feature_score = clip(feature_score, 0, 1)

    if multi:

        mol = inspiration.pop()

        while inspiration:
            mol = CombineMols(mol, inspiration.pop())

        protrude_dist = rdShapeHelpers.ShapeProtrudeDist(
            derivative, mol, allowReordering=False
        )

    else:
        protrude_dist = rdShapeHelpers.ShapeProtrudeDist(
            inspiration, derivative, allowReordering=False
        )

    protrude_dist = clip(protrude_dist, 0, 1)
    volume_score = 1 - protrude_dist

    SuCOS_score = (feature_score + volume_score) * 0.5

    if debug or print_scores:
        if debug:
            mrich.h3("scores")
        mrich.var("feature_score", feature_score)
        mrich.var("volume_score", volume_score)
        mrich.var("average_score", SuCOS_score)
        mrich.var("recapitulation_count", recapitulation_count)

    if return_all:
        return dict(
            average_score=SuCOS_score,
            feature_score=feature_score,
            volume_score=volume_score,
            recapitulation_count=recapitulation_count,
            recapitulation_fraction=recapitulation_fraction,
        )
    else:
        return SuCOS_score
