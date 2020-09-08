import glob
import pickle
import pandas as pd
import numpy as np
from typing import Tuple, List
from configuration import config as cfg
from utilities.featurizer import VscreenMLAggregateFeatures
from p_tqdm import p_map

predictors = [
    "fa_atr",
    "fa_rep",
    "fa_sol",
    "fa_elec",
    "hbond_bb_sc",
    "hbond_sc",
    "Interface_Energy",
    "Total_BSA",
    "Interface_HB",
    "Total_packstats",
    "Interface_unsat",
    "Total_pose_exposed_SASA",
    "interface_hydrophobic_sasa",
    "interface_polar_sasa",
    "6.6",
    "7.6",
    "8.6",
    "16.6",
    "6.7",
    "7.7",
    "8.7",
    "16.7",
    "6.8",
    "7.8",
    "8.8",
    "16.8",
    "6.9",
    "7.9",
    "8.9",
    "16.9",
    "6.15",
    "7.15",
    "8.15",
    "16.15",
    "6.16",
    "7.16",
    "8.16",
    "16.16",
    "6.17",
    "7.17",
    "8.17",
    "16.17",
    "6.35",
    "7.35",
    "8.35",
    "16.35",
    "6.53",
    "7.53",
    "8.53",
    "16.53",
    "SIDE_flex_ALPHA",
    "SIDE_flex_BETA",
    "SIDE_flex_OTHER",
    "BACK_flex_ALPHA",
    "BACK_flex_BETA",
    "BACK_flex_OTHER",
    "TotalElec",
    "TotalHbond",
    "Hphobe_contact",
    "Pi_Pi",
    "T-stack",
    "Cation-pi",
    "Salt_Bridge",
    "diff_entropy",
    "pienergy",
    "fsp3",
    "polarsurfacearea",
    "vdwsa",
]


def get_score_single_complex(path_tag_tuple: Tuple) -> float:
    """Get score for single complex.

    Parameters
    ----------
    path_tag_tuple: Tuple
        Tuple of path where complexes are stores and the tags for each complex
        complex files must be formatted as mini_complex_{tag}.pdb. Rosetta
        parameter file for small molecules must be formatted as {tag}.params

    Returns
    -------
    float
        vscreenml score

    """
    folder, filetag = path_tag_tuple
    data = (
        np.array(VscreenMLAggregateFeatures(folder, filetag).vscreenml_featurize())
        .reshape((1, -1))
        .astype(float)
    )
    with open(cfg.VSCREENML + "/model/XGB_CLASSIFIER_alldata.pickle.dat", "rb") as f:
        model = pickle.load(f)

    df = pd.DataFrame(data=data, columns=predictors)
    return model.predict_proba(df)[:, 1][0]


def get_vs_score(dir: str, ncpu: int = 4) -> List:
    """Get scores for a list of compounds in a given directory

    Parameters
    ----------
    dir: str
        Directory containing Rosetta minimized complexes and
        their corresponding .params files.
        The files must be follow the following convention
        mini_complex_{tag}.pdb, {tag}.params

    ncpu: int
        number of parallel processes to spin up

    Returns
    -------
    List[float]
        List of vscreen

    """

    list_of_files = glob.glob(dir + "/*.params")
    list_of_dir_tags = [(dir, i.split("/")[-1].split(".")[0]) for i in list_of_files]
    list_of_tags = [i[1] for i in list_of_dir_tags]
    result = p_map(get_score_single_complex, list_of_dir_tags, num_cpus=4)
    return pd.DataFrame(list(zip(list_of_tags, result)), columns=["complex", "score"])


if __name__ == "__main__":
    dir = "/Users/yusufadeshina/vscreenml/data"
    print(get_vs_score(dir))
