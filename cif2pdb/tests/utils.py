import re

from cif2pdb.convert import _convert_cif_to_pdb
from click.testing import CliRunner

runner = CliRunner()


# Utility functions

def dict2str(dict_):
    """
    Create string from Python dictionary.

    Parameters
    ----------
    dict_ : dict (str : str)

    Returns
    -------
    str
    """
    return " ".join([f'{k} {v}' for k, v in dict_.items()])


def generate_pdb_file(params, inpath, outpath):
    """
    Generate testing  PDB file.

    Parameters
    ----------
    params : str
        parameters of conversion
    inpath : str
        input path to CIF file
    outpath : str
        output path for PDB file
    """
    response = runner.invoke(_convert_cif_to_pdb, f"{params} -i {inpath} -o {outpath}")
    assert response.exit_code == 0


def split_pdb_atom_lines(lines):
    """
    Extract PDB atom properties
    for a list of strings and
    return a list of dictionaries.

    Parameters
    ----------
    lines : list of str
        list of PDB atom strings
    Returns
    -------
    list of dict
        list of PDB atom properties
    """
    return [{"atype": line[:6].strip(),
             "index": line[6:11].strip(),
             "atom": line[12:16].strip(),
             "resid": line[17:20].strip(),
             "chain": line[21:22].strip(),
             "resseq": int(line[22:26].strip()),
             "icode": line[26:27].strip(),
             "pos_x": line[30:38].strip(),
             "pos_y": line[38:46].strip(),
             "pos_z": line[46:54].strip(),
             "occ": line[54:60].strip(),
             "tfactor": line[60:66].strip(),
             "symbol": line[76:78].strip(),
             "charge": line[78:79].strip()}
            for line in lines]


def compare_pdb_files(file_1, file_2):
    """
    Compare two PDB files.

    Parameters
    ----------
    file_1 : str
        path to the first file
    file_2 : str
        path to the second file
    """
    # Load files
    with open(file_1, 'r') as f:
        atoms_1 = f.readlines()
    with open(file_2, 'r') as f:
        atoms_2 = f.readlines()

    # Consider only ATOM or HETATM lines
    r = re.compile("^(ATOM|HETATM)")
    atoms_1 = list(filter(r.match, atoms_1))
    atoms_2 = list(filter(r.match, atoms_2))

    assert len(atoms_1) == len(atoms_2)

    # Fetch PDB properties from each line
    atoms_1 = split_pdb_atom_lines(atoms_1)
    atoms_2 = split_pdb_atom_lines(atoms_2)

    # Sort by residue number and atom type
    # (sometimes order of atoms in two PDB files
    # is different within a given residue)
    atoms_1 = sorted(atoms_1, key=lambda d: (d['resseq'], d['atom']))
    atoms_2 = sorted(atoms_2, key=lambda d: (d['resseq'], d['atom']))

    # Compare all atoms in generated and expected files
    for i in range(len(atoms_1)):
        atom_1, atom_2 = atoms_1[i], atoms_2[i]
        for k in ('atype', 'atom', 'resid',
                  'icode', 'pos_x', 'pos_y', 'pos_z',
                  'occ', 'tfactor', 'symbol'):
            assert atom_1[k] == atom_2[k], f' ERROR for key "{k}", ' \
                                           f'id_1 {atom_1["index"]}: ' \
                                           f'id_2 {atom_2["index"]}: ' \
                                           f'{atom_1[k]} different from {atom_2[k]}'

