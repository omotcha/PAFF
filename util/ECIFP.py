"""
platform: any
env: any
name: ECIFP.py
ECIF-plus utils
"""
import numpy as np
import pandas as pd
import os
from rdkit import Chem
from scipy.spatial.distance import cdist
from itertools import product
from config import *
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator


class ECIFP:
    _ds_version = 2020

    # ECIF::atom_types
    _ProteinAtoms = ['C;4;1;3;0;0', 'C;4;2;1;1;1', 'C;4;2;2;0;0', 'C;4;2;2;0;1',
                     'C;4;3;0;0;0', 'C;4;3;0;1;1', 'C;4;3;1;0;0', 'C;4;3;1;0;1',
                     'C;5;3;0;0;0', 'C;6;3;0;0;0', 'N;3;1;2;0;0', 'N;3;2;0;1;1',
                     'N;3;2;1;0;0', 'N;3;2;1;1;1', 'N;3;3;0;0;1', 'N;4;1;2;0;0',
                     'N;4;1;3;0;0', 'N;4;2;1;0;0', 'O;2;1;0;0;0', 'O;2;1;1;0;0',
                     'S;2;1;1;0;0', 'S;2;2;0;0;0', 'OTH']

    _LigandAtoms_2020 = ['Br;1;1;0;0;0', 'C;3;3;0;1;1', 'C;4;1;1;0;0', 'C;4;1;2;0;0',
                         'C;4;1;3;0;0', 'C;4;2;0;0;0', 'C;4;2;1;0;0', 'C;4;2;1;0;1',
                         'C;4;2;1;1;1', 'C;4;2;2;0;0', 'C;4;2;2;0;1', 'C;4;3;0;0;0',
                         'C;4;3;0;0;1', 'C;4;3;0;1;1', 'C;4;3;1;0;0', 'C;4;3;1;0;1',
                         'C;4;4;0;0;0', 'C;4;4;0;0;1', 'C;5;3;0;0;0', 'C;5;3;0;1;1',
                         'C;6;3;0;0;0', 'Cl;1;1;0;0;0', 'F;1;1;0;0;0', 'I;1;1;0;0;0',
                         'N;3;1;0;0;0', 'N;3;1;1;0;0', 'N;3;1;2;0;0', 'N;3;2;0;0;0',
                         'N;3;2;0;0;1', 'N;3;2;0;1;1', 'N;3;2;1;0;0', 'N;3;2;1;0;1',
                         'N;3;2;1;1;1', 'N;3;3;0;0;0', 'N;3;3;0;0;1', 'N;3;3;0;1;1',
                         'N;4;1;2;0;0', 'N;4;1;3;0;0', 'N;4;2;1;0;0', 'N;4;2;2;0;0',
                         'N;4;2;2;0;1', 'N;4;3;0;0;0', 'N;4;3;0;0;1', 'N;4;3;1;0;0',
                         'N;4;3;1;0;1', 'N;4;4;0;0;0', 'N;4;4;0;0;1', 'N;5;2;0;0;0',
                         'N;5;3;0;0;0', 'N;5;3;0;1;1', 'O;2;1;0;0;0', 'O;2;1;1;0;0',
                         'O;2;2;0;0;0', 'O;2;2;0;0;1', 'O;2;2;0;1;1', 'P;5;4;0;0;0',
                         'P;6;4;0;0;0', 'P;6;4;0;0;1', 'P;7;4;0;0;0', 'S;2;1;0;0;0',
                         'S;2;1;1;0;0', 'S;2;2;0;0;0', 'S;2;2;0;0;1', 'S;2;2;0;1;1',
                         'S;3;3;0;0;0', 'S;3;3;0;0;1', 'S;4;3;0;0;0', 'S;6;4;0;0;0',
                         'S;6;4;0;0;1', 'S;7;4;0;0;0', 'OTH']

    # ELEMENTS:atom_types
    _ProteinELEMENTS = ["C", "N", "O", "S"]
    _LigandELEMENTS = ["Br", "C", "Cl", "F", "I", "N", "O", "P", "S"]

    # {ECIF|ELEMENTS}::LD
    _LigandDescriptors = ['MaxEStateIndex', 'MinEStateIndex', 'MaxAbsEStateIndex', 'MinAbsEStateIndex',
                          'qed', 'MolWt', 'HeavyAtomMolWt', 'ExactMolWt', 'NumValenceElectrons',
                          'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'BalabanJ',
                          'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n',
                          'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'HallKierAlpha', 'Kappa1',
                          'Kappa2', 'Kappa3', 'LabuteASA', 'PEOE_VSA14', 'SMR_VSA1', 'SMR_VSA10',
                          'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7',
                          'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12',
                          'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6',
                          'SlogP_VSA7', 'SlogP_VSA8', 'TPSA', 'EState_VSA1', 'EState_VSA10',
                          'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5',
                          'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'VSA_EState1',
                          'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5',
                          'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'FractionCSP3',
                          'HeavyAtomCount', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles',
                          'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles',
                          'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 'NumHDonors',
                          'NumHeteroatoms', 'NumRotatableBonds', 'NumSaturatedCarbocycles',
                          'NumSaturatedHeterocycles', 'NumSaturatedRings', 'RingCount', 'MolLogP',
                          'MolMR', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_N',
                          'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO',
                          'fr_C_S', 'fr_HOCCN', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O',
                          'fr_Ndealkylation1', 'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde',
                          'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide',
                          'fr_amidine', 'fr_aniline', 'fr_aryl_methyl', 'fr_azo', 'fr_barbitur',
                          'fr_benzene', 'fr_bicyclic', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester',
                          'fr_ether', 'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone',
                          'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone',
                          'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine',
                          'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitroso', 'fr_oxazole',
                          'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol', 'fr_phenol_noOrthoHbond',
                          'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_pyridine', 'fr_quatN',
                          'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole',
                          'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_urea']

    _DescCalc = MolecularDescriptorCalculator(_LigandDescriptors)

    def _get_possible_ecif(self):
        return [i[0]+"-"+i[1] for i in product(
            self._ProteinAtoms, self.__getattribute__('_LigandAtoms_{}'.format(self._ds_version)))]

    def _get_possible_element(self):
        return [i[0]+"-"+i[1] for i in product(self._ProteinELEMENTS, self._LigandELEMENTS)]

    def _get_atom_type(self, atom):
        """

        :param atom: RDKit supported atom type
        :return:
          an ECIF::atom_types
        #     Atom symbol;
        #     Explicit valence;
        #     Attached heavy atoms;
        #     Attached hydrogens;
        #     Aromaticity;
        #     Ring membership
        """
        symbol = atom.GetSymbol()
        # add escape element symbol OTH that stands for other types
        if symbol not in self._LigandELEMENTS:
            symbol = "OTH"
        AtomType = [symbol,
                    str(atom.GetExplicitValence()),
                    str(len([x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() != "H"])),
                    str(len([x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() == "H"])),
                    str(int(atom.GetIsAromatic())),
                    str(int(atom.IsInRing())),
                    ]

        return ";".join(AtomType)

    # helpers for loading ligands and proteins

    def _load_ligand(self, sdf):
        """
        This function takes an SDF for a ligand as input and returns it as a pandas DataFrame with its atom types labeled according to ECIF
        :param sdf: ligand file
        :return: a pandas DataFrame for the ligand with ECIF::atom_types
        """
        # This function takes an SDF for a ligand as input and returns it as a pandas DataFrame

        m = Chem.MolFromMolFile(sdf, sanitize=False)
        m.UpdatePropertyCache(strict=False)

        ECIF_atoms = []

        for atom in m.GetAtoms():
            if atom.GetSymbol() != "H":  # Include only non-hydrogen atoms
                entry = [int(atom.GetIdx())]
                entry.append(self._get_atom_type(atom))
                pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                entry.append(float("{0:.4f}".format(pos.x)))
                entry.append(float("{0:.4f}".format(pos.y)))
                entry.append(float("{0:.4f}".format(pos.z)))
                ECIF_atoms.append(entry)

        df = pd.DataFrame(ECIF_atoms)
        df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]
        df["ECIF_ATOM_TYPE"].fillna("OTH")
        return df

    def _load_protein(self, pdb):
        """
        This function takes a PDB for a protein as input and returns it as a pandas DataFrame with its atom types labeled according to ECIF
        :param pdb: protein file
        :return: a pandas DataFrame for the protein with ECIF::atom_types
        """
        Atom_Keys = pd.read_csv("../ECIF_PDB_Atom_Keys.csv", sep=",")
        ECIF_atoms = []

        f = open(pdb)
        for i in f:
            if i[:4] == "ATOM":
                # Include only non-hydrogen atoms
                if (len(i[12:16].replace(" ", "")) < 4 and i[12:16].replace(" ", "")[0] != "H") or (
                        len(i[12:16].replace(" ", "")) == 4 and i[12:16].replace(" ", "")[1] != "H" and
                        i[12:16].replace(" ", "")[0] != "H"):
                    ECIF_atoms.append([int(i[6:11]),
                                       i[17:20] + "-" + i[12:16].replace(" ", ""),
                                       float(i[30:38]),
                                       float(i[38:46]),
                                       float(i[46:54])
                                       ])

        f.close()

        df = pd.DataFrame(ECIF_atoms, columns=["ATOM_INDEX", "PDB_ATOM", "X", "Y", "Z"])
        df = df.merge(Atom_Keys, left_on='PDB_ATOM', right_on='PDB_ATOM')[
            ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]].sort_values(by="ATOM_INDEX").reset_index(drop=True)
        df["ECIF_ATOM_TYPE"].fillna("OTH")
        return df

    def _get_pl_pairs(self, protein_f, ligand_f, distance_cutoff=6.0):
        """
        This function returns the protein-ligand atom-type pairs for a given distance cutoff
        :param protein_f: pdb file name with dir
        :param ligand_f:  sdf file name with dir
        :param distance_cutoff:
        :return:
        """

        Target = self._load_protein(protein_f)
        Ligand = self._load_ligand(ligand_f)

        for i in ["X", "Y", "Z"]:
            Target = Target[Target[i] < float(Ligand[i].max()) + distance_cutoff]
            Target = Target[Target[i] > float(Ligand[i].min()) - distance_cutoff]

        # Get all possible pairs
        Pairs = list(product(Target["ECIF_ATOM_TYPE"], Ligand["ECIF_ATOM_TYPE"]))
        Pairs = [x[0] + "-" + x[1] for x in Pairs]
        Pairs = pd.DataFrame(Pairs, columns=["ECIF_PAIR"])
        Distances = cdist(Target[["X", "Y", "Z"]], Ligand[["X", "Y", "Z"]], metric="euclidean")
        Distances = Distances.reshape(Distances.shape[0] * Distances.shape[1], 1)
        Distances = pd.DataFrame(Distances, columns=["DISTANCE"])

        Pairs = pd.concat([Pairs, Distances], axis=1)
        Pairs = Pairs[Pairs["DISTANCE"] <= distance_cutoff].reset_index(drop=True)
        # Pairs from ELEMENTS could be easily obtained froms pairs from ECIF
        Pairs["ELEMENTS_PAIR"] = [x.split("-")[0].split(";")[0] + "-" + x.split("-")[1].split(";")[0] for x in
                                  Pairs["ECIF_PAIR"]]
        return Pairs

    # callables
    def get_ligand_features_by_file(self, ligand_f):
        """
        calculate ligand descriptors using RDKit
        :param ligand_f: sdf file name with dir
        :return:
        """
        ligand = Chem.MolFromMolFile(ligand_f, sanitize=False)
        ligand.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(ligand)
        return self._DescCalc.CalcDescriptors(ligand)

    def get_ecif(self, protein_f, ligand_f, distance_cutoff=6.0):
        """
        get the fingerprint-like array (ECIF) for a protein-ligand pair
        :param protein_f: pdb file name with dir
        :param ligand_f: sdf file name with dir
        :param distance_cutoff:
        :return:
        """
        Pairs = self._get_pl_pairs(protein_f, ligand_f, distance_cutoff=distance_cutoff)
        ECIF = [list(Pairs["ECIF_PAIR"]).count(x) for x in self._get_possible_ecif()]
        return ECIF

    def get_elements(self, protein_f, ligand_f, distance_cutoff=6.0):
        """
        get the fingerprint-like array (ELEMENTS) for a protein-ligand pair
        :param protein_f: pdb file name with dir
        :param ligand_f: sdf file name with dir
        :param distance_cutoff:
        :return:
        """
        Pairs = self._get_pl_pairs(protein_f, ligand_f, distance_cutoff=distance_cutoff)
        ELEMENTS = [list(Pairs["ELEMENTS_PAIR"]).count(x) for x in self._get_possible_element()]
        return ELEMENTS

    def get_possible_ecif(self):
        return self._get_possible_ecif()

    def get_possible_elements(self):
        return self._get_possible_element()

    def get_ligand_descriptors(self):
        return self._LigandDescriptors

    def __init__(self, version=2020):
        self._ds_version = version


def my_test():
    """
    test
    :return:
    """
    ligand_file = os.path.join(tmpdata_dir, '4glr_ligand.sdf')
    protein_file = os.path.join(tmpdata_dir, '4glr_protein.pdb')
    ecifp_helper = ECIFP(2020)
    print(len(ecifp_helper.get_ecif(protein_file, ligand_file, distance_cutoff=6.0)))


if __name__ == '__main__':
    my_test()
