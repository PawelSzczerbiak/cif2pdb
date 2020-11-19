# =====================================================
# Script calculates number of residues in PDB file
# =====================================================

import sys
import os

from os.path import join
from os.path import normpath


AA_SEQ_DICT = {"ALA": "A",
               "ARG": "R",
               "ASN": "N",
               "ASP": "D",
               "ASX": "B",
               "CYS": "C",
               "GLU": "E",
               "GLN": "Q",
               "GLX": "Z",
               "GLY": "G",
               "HIS": "H",
               "ILE": "I",
               "LEU": "L",
               "LYS": "K",
               "MET": "M",
               "PHE": "F",
               "PRO": "P",
               "SER": "S",
               "THR": "T",
               "TRP": "W",
               "TYR": "Y",
               "VAL": "V",
               "SEC": "U",
               "PYL": "O",
               "UNK": "X",
               }


def fetch_residues_from_pdb_file(path):
    """
    Fetch atom residues, their names and corresponding indices
    from PDB file (ATOM flag).

    Parameters
    ----------
    path : str
        path to PDB file

    Returns
    -------
    residues : list of str
        e.g. "2", "3" etc.
    names : list of int
        e.g. "ALA", "ASN" etc.
    indices : list of int
        e.g. 0, 1 et.c
    """
    residues, names, indices = [], [], []
    with open(path, 'r') as f:
        data = f.readlines()
        residues, names, indices = fetch_residues_from_pdb_list(data)
    return residues, names, indices


def fetch_residues_from_pdb_list(list_atoms):
    """
    Fetch atom residues, their names and corresponding indices
    from list of PDB atoms.

    Parameters
    ----------
    list_atoms : list of str
        list of PDB atoms

    Returns
    -------
    residues : list of str
        e.g. "2", "3" etc.
    names : list of int
        e.g. "ALA", "ASN" etc.
    indices : list of int
        e.g. 0, 1 et.c
    """
    residues, names, indices = [], [], []
    for i, line in enumerate(list_atoms):
        if line.startswith("ATOM"):
            residues.append(line[22:26].strip())
            names.append(line[17:20].strip())
            indices.append(i)
    return residues, names, indices


def transform_ranges(ranges):
    """
    Transform sequence of ranges into a list
    of sorted integeres.

    Parameteres
    -----------
    ranges : str
        ranges separated by colon (,) e.g.:
            100,1-5,2,300,200-202

    Returns
        list of integers
        In the example above we would get
        [1, 2, 3, 4, 5, 100, 200, 201, 202, 300]
    -------
    """
    ids = set()
    ranges = ranges.split(',')
    for range_ in ranges:
        el = range_.split('-')
        if len(el) == 1:
            ids.add(int(el[0]))
        else:
            # take care of negative left limit
            a, b = '-'.join(el[:-1]), el[-1]
            [ids.add(i) for i in range(int(a), int(b) + 1)]
    ids = sorted(list(ids))
    return ids


def get_atoms_for_residues(list_atoms, list_residues):
    """
    Get atoms from list of PDB atoms (ATOM flag)
    that belong to specfied residues.

    Parameters
    ----------
    list_atoms : list of str
        list of PDB atoms

    list_residues : list of int
        list of residues e.g.:
            [1, 5, 7, 200]

    Returns
    -------
    atoms : list of str
    """
    atoms = []
    residues, _, indices = fetch_residues_from_pdb_list(list_atoms)
    for i, residue in enumerate(residues):
        if int(residue) in list_residues:
            atoms.append(list_atoms[indices[i]])
    return atoms


def get_number_of_residues(residues):
    """
    Return number of residues.

    Parameters
    ----------
    residues : list
        can be list or strings e.g.:
            ["1", "1", "2", "3", "3"]
        or sequenece of chars e.g.:
            "aaacccdddaaa"
    Returns
    -------
    int
        Number of elements changes in the list. In the example
        above we would get 3 and 4 respectively.
    """
    prev = residues[0]
    number = 1
    for i in range(len(residues) - 1):
        cont = residues[i + 1]
        if prev != cont:
            number += 1
            prev = cont
    return number


# TODO: generalize to include residue ranges and PDB list
def get_sequence_for_pdb_file(path):
    """
    Return sequence from a PDB file.

    Parameters
    ----------
    path : str
        path to PDB file
    Returns
    -------
    seq : str
        Sequence of aminoacids e.g. AAKBCA etc.
    """
    residues, names, _ = fetch_residues_from_pdb_file(path)
    seq = AA_SEQ_DICT[names[0]]
    prev = residues[0]
    for i in range(len(residues) - 1):
        cont = residues[i + 1]
        if prev != cont:
            seq += AA_SEQ_DICT[names[i + 1]]
            prev = cont
    return seq


# TODO: rewrite using Click
if __name__ == "__main__":

    assert len(sys.argv) == 6, "{0} ERROR: wrong number of arguments".format(sys.argv[1:])

    input_path = sys.argv[1]  # full input path to the PDB files
    output_path = sys.argv[2]  # output path
    output_file = sys.argv[3]  # output file
    limits = sys.argv[4]  # residue number limits (comma separated)
    indicator = sys.argv[5]  # which input path index to use as a structure name indicator
    # counting from the end e.g. for path/to/file.pdb we'd get:
    # "file.pdb" for -1
    # "to"       for -2
    # "path"     for -3

    limits = list(map(lambda x: int(x), limits.split(",")))

    residues, _, _ = fetch_residues_from_pdb_file(input_path)
    number = get_number_of_residues(residues)

    assert len(limits) >= 1, "{0} ERROR: wrong number of limits".format(limits)

    indicator = int(indicator)
    indicator_name = normpath(input_path).split(os.sep)[indicator]
    line = f"{indicator_name} {number}\n"

    if number < limits[0]:
        with open(join(output_path,
                       f"{output_file}__less_{limits[0]}"), 'a') as f:
            f.write(line)
    elif number >= limits[-1]:
        with open(join(output_path,
                       f"{output_file}__geq_{limits[-1]}"), 'a') as f:
            f.write(line)
    else:
        for i in range(len(limits[1:])):
            if limits[i] <= number < limits[i + 1]:
                with open(join(output_path,
                               f"{output_file}__geq_{limits[i]}_less_{limits[i + 1]}"), 'a') as f:
                    f.write(line)
