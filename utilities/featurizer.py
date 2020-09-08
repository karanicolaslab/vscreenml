#! usr/bin/python
import os
import re
import subprocess
import pandas as pd
from configuration import config as cfg
from rfscore_like_feature import AtomCountFeaturizer
from typing import Tuple, List

__author__ = "Yusuf Adeshina"
__email__ = "dryusufadeshina@gmail.com"

# TODO: Change string path to PoxisPath objects


def parse_binana_file(binana_file: str) -> Tuple[float]:
    """Parse the output file of binana.

    Parameters
    ----------
    binana_file: str
        binana output file

    Returns
    -------
    Tuple[float]
        Tuple of the extracted features

    """
    with open(binana_file, "r") as f:
        content = f.read()

    data = {}

    for block in content.split("\n\n"):
        block = block.lstrip()
        block = block.rstrip()

        if re.search("Atom-type pair counts within 2.5 angstroms:", block):
            for line in block.split("\n"):
                matchObj = re.search(r"Atom|.*----", line)
                if not matchObj:
                    k, v, p = line.split("|")
                    k = k.strip()
                    v = v.strip()
                    p = p.strip()

                    k_v = str(k) + "_" + "2.5A" + "_" + str(v)
                    data[k_v] = float(p)

        if re.search("Atom-type pair counts within 4.0 angstroms:", block):
            for line in block.split("\n"):
                matchObj = re.search(r"Atom|.*----", line)
                if not matchObj:
                    k, v, p = line.split("|")
                    k = k.strip()
                    v = v.strip()
                    p = p.strip()
                    k_v = str(k) + "_" + "4.0A" + "_" + str(v)
                    data[k_v] = float(p)

        if re.search("Summed electrostatic energy by atom-type pair", block):
            TotalElec = 0.0
            for line in block.split("\n"):
                matchObj = re.search(r"Summed|Atom|.*----", line)
                if not matchObj:
                    k, v, p = line.split("|")
                    k = k.strip()
                    v = v.strip()
                    p = p.strip()
                    k_v = str(k) + "_" + "Elec" + "_" + str(v)
                    TotalElec += float(p)
            data["TotalElec"] = TotalElec

        if re.search("Number of rotatable bonds in the ligand:", block):
            for line in block.split("\n"):

                v, p = line.split(":")
                v = v.strip()
                p = p.strip()
                k_v = "RotBond"

                data[k_v] = float(p)

        if re.search("Active-site flexibility:", block):
            for line in block.split("\n"):
                matchObj = re.search(r"Active-site|Sidechain/Backbone|.*----", line)
                if not matchObj:
                    k, v, p = line.split("|")
                    k = k.strip()
                    v = v.strip()
                    p = p.strip()

                    k_v = str(k) + "_" + "flex" + "_" + str(v)
                    data[k_v] = float(p)

        if re.search("Hydrogen bonds:", block):
            TotalHbond = 0.0
            for line in block.split("\n"):
                matchObj = re.search(r"Hydrogen bonds|Location of Donor|.*----", line)
                if not matchObj:
                    r, k, v, p = line.split("|")
                    r = r.strip()
                    k = k.strip()
                    v = v.strip()
                    p = p.strip()

                    k_v = str(r) + "_" + str(k) + "_" + str(v) + "Hbond"
                    TotalHbond += float(p)
            data["TotalHbond"] = TotalHbond

        if re.search("Hydrophobic contacts", block):
            TotalHphob = 0.0
            for line in block.split("\n"):
                matchObj = re.search(
                    r"Hydrophobic contacts|Sidechain/Backbone|.*----", line
                )
                if not matchObj:
                    k, v, p = line.split("|")
                    k = k.strip()
                    v = v.strip()
                    p = p.strip()

                    k_v = str(k) + "_" + "Hphob" + "_" + str(v)
                    TotalHphob += float(p)
            data["Hphob"] = TotalHphob

        if re.search("pi-pi stacking interactions:", block):
            TotalPipi = 0.0
            for line in block.split("\n"):
                matchObj = re.search(r"pi-pi|Secondary Structure|.*----", line)
                if not matchObj:
                    v, p = line.split("|")
                    v = v.strip()
                    p = p.strip()
                    k_v = str(v) + "_" + "Pipi"
                    TotalPipi += float(p)
            data["Pipi"] = TotalPipi

        if re.search("T-stacking (face-to-edge) interactions:", block):
            TotalTstac = 0.0
            for line in block.split("\n"):
                matchObj = re.search(r"T-stacking|Secondary Structure|.*----", line)
                if not matchObj:
                    v, p = line.split("|")
                    v = v.strip()
                    p = p.strip()

                    k_v = str(v) + "_" + "Tstac"
                    TotalTstac += float(p)
            data["Tstac"] = TotalTstac

        if re.search("Cation-pi interactions:", block):
            TotalCapi = 0.0
            for line in block.split("\n"):
                matchObj = re.search(r"Cation-pi|Secondary Structure|.*----", line)
                if not matchObj:
                    k, v, p = line.split("|")
                    k = k.strip()
                    v = v.strip()
                    p = p.strip()
                    k_v = str(k) + "_" + "Capi" + "_" + str(v)
                    TotalCapi += float(p)
            data["Capi"] = TotalCapi

        if re.search("Salt Bridges:", block):
            TotalSalt = 0.0
            for line in block.split("\n"):
                matchObj = re.search(r"Salt Bridges:|Secondary Structure|.*----", line)
                if not matchObj:
                    v, p = line.split("|")
                    v = v.strip()
                    p = p.strip()
                    k_v = str(v) + "_" + "Salt"
                    TotalSalt += float(p)
            data["Salt"] = TotalSalt

    return (
        data.get("SIDECHAIN_flex_ALPHA", 0.0),
        data.get("SIDECHAIN_flex_BETA", 0.0),
        data.get("SIDECHAIN_flex_OTHER", 0.0),
        data.get("BACKBONE_flex_ALPHA", 0.0),
        data.get("BACKBONE_flex_BETA", 0.0),
        data.get("BACKBONE_flex_OTHER", 0.0),
        data.get("TotalElec", 0.0),
        data.get("TotalHbond", 0.0),
        data.get("Hphob", 0.0),
        data.get("Pipi", 0.0),
        data.get("Tstac", 0.0),
        data.get("Capi", 0.0),
        data.get("Salt", 0.0),
    )


class VscreenMLFeaturizer(object):
    "This featurizes protein-ligand complexes for vscreenml scoring"

    def __init__(self, folder: str, filetag: str) -> None:
        """Initialize class.

        Parameters
        ----------
        folder: str
            string of the path of the working directory

        filetag: str
            file tag. The input pdbs should be formatted
            as mini_complex_{tagfile}.pdb

        Returns
        -------
        None

        """
        self.filetag = filetag
        self.folder = folder
        self.feature_list = []

    def prepare_protein_and_ligand_for_scoring(self) -> None:
        """"Extract ligands and protein structures."""

        p = subprocess.Popen(
            'grep -E "ATOM|ZN|MG|MN" '
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + ".pdb > "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_unb.pdb && grep HETATM "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + ".pdb|grep -v 'ZN'|grep -v 'MG'|grep -v 'MN' > "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_lig.pdb",
            shell=True,
        )

        p.communicate()

        # reformatting the ligands
        p = subprocess.Popen(
            cfg.BABEL
            + "/bin/obabel "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_lig.pdb -O "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_lig1.sdf",  # >/dev/null 2>&1",
            shell=True,
        )

        p = subprocess.Popen(
            cfg.CHEMAXON
            + '/bin/structurecheck -c "valenceerror" -m fix -f sdf -o '
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_lig.sdf "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_lig1.sdf",
            shell=True,
        )
        p.communicate()

    def rosetta_minimization(self) -> Tuple[str]:
        #
        """Score Rosetta pre-minimized complex

        Returns
        -------
        Tuple[str]
            Tuple of Rosetta interface scores and other non-canonical rosetta
            features

        """
        try:
            p_rosetta = subprocess.Popen(
                cfg.ROSETTA_DIR
                + "/main/source/bin/minimize_ppi."
                + cfg.ROSETTA_BUILD
                + " -s "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + ".pdb -extra_res_fa "
                + self.folder
                + "/"
                + self.filetag
                + ".params -database "
                + cfg.ROSETTA_DIR
                + "/main/database -ignore_zero_occupancy false -score_only \
                |grep 'Interface_Scores:mini' |awk '{print $4,$5,$6,$7,$8}' ",
                shell=True,
                stdout=subprocess.PIPE,
            )
            rosetta = p_rosetta.communicate()[0].split()

        except OSError:
            rosetta = ["-23.6933", "719.7050", "1.0000", "0.6410", "4.0000"]

        return rosetta[0], rosetta[1], rosetta[2], rosetta[3], rosetta[4]

    def calculate_sasa(self) -> Tuple[str]:
        """Calculate theta-lig features.

        Returns
        -------
        Tuple of theta-lig features

        """
        try:
            p_sasa = subprocess.Popen(
                cfg.ROSETTA_DIR
                + "/main/source/bin/ragul_get_ligand_sasa."
                + cfg.ROSETTA_BUILD
                + " -database "
                + cfg.ROSETTA_DIR
                + "/main/database -input_bound_pdb "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + ".pdb -input_unbound_pdb "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + "_unb.pdb -input_ligand_pdb "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + "_lig.pdb -extra_res_fa "
                + self.folder
                + "/"
                + self.filetag
                + ".params |grep 'Scores:' |awk '{print $2,$3,$4}'",
                shell=True,
                stdout=subprocess.PIPE,
            )
            sasa = p_sasa.communicate()[0].split()

        except OSError:
            # use typical value (median) from the training set to impute missing value
            sasa = ["0.3304", "0.7138", "0.2867"]

        return sasa[0], sasa[1], sasa[2]

    def calculate_component_score(self) -> Tuple[str]:
        """Calculate Rosetta componenet scores

        Returns
        -------
        Tuple[str]
            Tuple of components scores

        """
        try:
            p_compo = subprocess.Popen(
                cfg.ROSETTA_DIR
                + "/main/source/bin/vscrenml_interface_ddg."
                + cfg.ROSETTA_BUILD
                + " -database "
                + cfg.ROSETTA_DIR
                + "/main/database -s "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + ".pdb -extra_res_fa "
                + self.folder
                + "/"
                + self.filetag
                + ".params -ignore_zero_occupancy false  \
                 |grep 'Component_score:' |awk '{print $4,$5,$6,$10,$11,$12}'",
                shell=True,
                stdout=subprocess.PIPE,
            )
            compo = p_compo.communicate()[0].split()

        except OSError:
            # use typical value (median) from the training set to impute missing value
            compo = ["-28.0043", "6.45592", "7.77609", "-4.8185", "0.0000", "0.0000"]

        return compo[0], compo[1], compo[2], compo[3], compo[4], compo[5]

    def generate_dot_pdbqt_files(self) -> None:
        """Generate .pdbqt filetags for cfg.BINANA features"""

        p = subprocess.Popen(
            cfg.MGLTOOLS
            + "/bin/pythonsh "
            + cfg.MGLTOOLS
            + "/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_unb.pdb -A hydrogens -o "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_unb.pdbqt >/dev/null 2>&1",
            shell=True,
        )
        p.wait()

        p = subprocess.Popen(
            cfg.MGLTOOLS
            + "/bin/pythonsh "
            + cfg.MGLTOOLS
            + "/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_lig.pdb -A hydrogens -o "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_lig.pdbqt >/dev/null 2>&1",
            shell=True,
        )
        p.wait()

    def generate_binana_features(self) -> Tuple[str]:
        """Calculate binana features.

        Returns
        -------
        Tuple of binana features

        """
        try:
            p = subprocess.Popen(
                "python "
                + cfg.BINANA
                + "/binana.py  -receptor "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + "_unb.pdbqt -ligand "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + "_lig.pdbqt > "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + "_binana.out",
                shell=True,
            )
            p.communicate()

            binana = parse_binana_file(
                self.folder + "/mini_complex_" + self.filetag + "_binana.out"
            )

        except (OSError, FileNotFoundError):
            # use typical value (median) from the training set to impute missing value
            binana = [
                "5",
                "22",
                "17",
                "0",
                "4",
                "5",
                "-56280.936",
                "1",
                "21",
                "0",
                "0",
                "0",
                "0",
            ]
        return (
            binana[0],
            binana[1],
            binana[2],
            binana[3],
            binana[4],
            binana[5],
            binana[6],
            binana[7],
            binana[8],
            binana[9],
            binana[10],
            binana[11],
            binana[12],
        )

    def calculate_szybki(self) -> str:
        # generating cfg.OMEGA conformers for entropy calculations
        """Calculate Szybki features.

        Returns
        -------
        Tuple of Szybki features

        """
        p = subprocess.Popen(
            cfg.OMEGA
            + " -in "
            + self.folder
            + "/mini_complex_"
            + self.filetag
            + "_lig1.sdf -out "
            + self.folder
            + "/"
            + self.filetag
            + ".oeb.gz -strictstereo false -strictatomtyping false -warts true -maxconfs 1000 -rms 0.1  >/dev/null 2>&1",
            shell=True,
        )

        try:
            # p = subprocess.Popen(
            #    cfg.SZYBKI
            #    + " -entropy AN  -prefix "
            #    + self.folder
            #    + "/bou_"
            #    + self.filetag
            #    + " -complex "
            #    + self.folder
            #    + "/mini_complex_"
            #    + self.filetag
            #    + ".pdb -strict false >/dev/null 2>&1",
            #    shell=True,
            #    stdout=subprocess.PIPE,
            # )

            # versions other than 1.9.0.3 follows this convention
            p = subprocess.Popen(
                cfg.SZYBKI
                + " -entropy AN  -prefix "
                + self.folder
                + "/bou_"
                + self.filetag
                + " -p "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + "_unb.pdb -strict false"
                + " "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + "_lig.sdf ",  # >/dev/null 2>&1",
                shell=True,
                stdout=subprocess.PIPE,
            )
            p.communicate()

            p = subprocess.Popen(
                cfg.SZYBKI
                + " -entropy AN  -prefix "
                + self.folder
                + "/unb_"
                + self.filetag
                + " -sheffield -in "
                + self.folder
                + "/"
                + self.filetag
                + ".oeb.gz -strict false >/dev/null 2>&1",
                shell=True,
                stdout=subprocess.PIPE,
            )
            p.communicate()

            p_bound = subprocess.Popen(
                "grep 'At T=298.15 K,' "
                + self.folder
                + "/bou_"
                + self.filetag
                + ".log|awk '{print $8}'|tail -1",
                shell=True,
                stdout=subprocess.PIPE,
            )
            bound = p_bound.communicate()[0].split()[0]
            p_unbound = subprocess.Popen(
                "grep 'At T=298.15 K,' "
                + self.folder
                + "/unb_"
                + self.filetag
                + ".log|awk '{print $8}'|tail -1",
                shell=True,
                stdout=subprocess.PIPE,
            )
            unbound = p_unbound.communicate()[0].split()[0]
            entropy = float(bound) - float(unbound)
            entropy = "0.0"
        except OSError:
            # use typical value (median) from the training set to impute missing value
            entropy = "-1.4"
        return entropy

    def calculate_chemaxon_features(self) -> Tuple[str]:
        """Calculate Chemaxon ligand descriptors features.

        Returns
        -------
        Tuple[str]
            Tuple of ligand descriptors

        """

        try:
            p_chemax = subprocess.Popen(
                cfg.CXCALC
                + " -N i logp pienergy acceptorcount donorcount fsp3 polarsurfacearea vdwsa "
                + self.folder
                + "/mini_complex_"
                + self.filetag
                + "_lig1.sdf |tail -1",
                shell=True,
                stdout=subprocess.PIPE,
            )
            chemax = p_chemax.communicate()[0].split()
        except OSError:
            # use typical value (median) from the training set to impute missing value
            chemax = ["FAIL", "39.46", "FAIL", "FAIL", "0.28", "79.90", "468.78"]

        return chemax[1], chemax[4], chemax[5], chemax[6]

    def calculate_rfscore_like_features(self) -> List[str]:
        filename = self.folder + "/mini_complex_" + self.filetag + ".pdb"
        feat_dict = AtomCountFeaturizer(filename).calc_feature()
        feat_df = pd.DataFrame.from_dict(feat_dict, orient="index")
        return feat_df.T.values.tolist()[0]

    def clean_up(self) -> None:
        """Clean up by deleting files."""
        os.remove(self.folder + "/" + self.filetag + ".params")
        os.remove(self.folder + "/" + self.filetag + "_0001.pdb")
        os.remove(self.folder + "/mini_complex_" + self.filetag + "_lig1.sdf")
        os.remove(self.folder + "/mini_complex_" + self.filetag + "_binana.out")
        os.remove(self.folder + "/mini_complex_" + self.filetag + "_unb.pdbqt")
        os.remove(self.folder + "/mini_complex_" + self.filetag + "_lig.pdbqt")
        os.remove(self.folder + "/mini_complex_" + self.filetag + "_unb.pdb")
        os.remove(self.folder + "/mini_complex_" + self.filetag + "_lig.pdb")
        os.remove(self.folder + "/" + self.filetag + ".oeb.gz")
        os.remove(self.folder + "/bou_" + self.filetag + ".log")
        os.remove(self.folder + "/bou_" + self.filetag + "_out.oeb")
        os.remove(self.folder + "/bou_" + self.filetag + ".status")
        os.remove(self.folder + "/bou_" + self.filetag + ".param")
        os.remove(self.folder + "/unb_" + self.filetag + ".log")
        os.remove(self.folder + "/unb_" + self.filetag + "_out.oeb")
        os.remove(self.folder + "/unb_" + self.filetag + ".status")
        # os.remove(self.folder + "/unb_" + self.filetag + ".param")


class VscreenMLAggregateFeatures(object):
    "This featurizes protein-ligand complexes for vscreenml scoring"

    def __init__(self, folder: str, filetag: str, clean_up: bool = False) -> None:
        """Initialize class.

        Parameters
        ----------
        folder: str
            string of the path of the working directory

        filetag: str
            file tag. The input pdbs should be formatted
            as mini_complex_{tagfile}.pdb

        Returns
        -------
        None

        """
        self.filetag = filetag
        self.folder = folder
        self.feature_list: List = []
        self.clean_up = clean_up

    def vscreenml_featurize(self) -> List[str]:
        """Get overall vscreenml features.

        Returns
        -------
        List[str]
            Tuple of the extracted features

        """
        features = VscreenMLFeaturizer(self.folder, self.filetag)

        features.prepare_protein_and_ligand_for_scoring()
        self.feature_list.extend(features.calculate_component_score())
        self.feature_list.extend(features.rosetta_minimization())
        self.feature_list.extend(features.calculate_sasa())
        self.feature_list.extend(features.calculate_rfscore_like_features())
        features.generate_dot_pdbqt_files()
        self.feature_list.extend(features.generate_binana_features())
        self.feature_list.append(features.calculate_szybki())
        self.feature_list.extend(features.calculate_chemaxon_features())
        if self.clean_up:
            features.clean_up()
        return self.feature_list


if __name__ == "__main__":
    folder = "/Users/yusufadeshina/vscreenml/data"
    filetag = "Z2488751527_29"
    feat = VscreenMLFeaturizer(folder, filetag, False)
    feat.generate_binana_features()
