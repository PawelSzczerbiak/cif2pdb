import os
import click

from os.path import join
from os.path import normpath

import cif2pdb.convert


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


def fetch_residues_from_pdb_file(path, atoms_only=False):
    """
    Fetch atom residues, their names and corresponding indices
    (line numbers) from PDB file (ATOM/HETATM flag).

    Parameters
    ----------
    path : str
        path to PDB file
    atoms_only : bool
        fetch only atoms (ATOM)

    Returns
    -------
    residues : list of str
        e.g. "2", "3" etc.
    names : list of int
        e.g. "ALA", "ASN" etc.
    indices : list of int
        e.g. 0, 1 etc.
    """
    with open(path, 'r') as f:
        data = f.readlines()
        return fetch_residues_from_pdb_list(data, atoms_only)


def fetch_residues_from_pdb_list(list_atoms, atoms_only=False):
    """
    Fetch atom residues, their names and corresponding indices
    (line numbers) from list of PDB atoms (ATOM/HETATM flag).

    Parameters
    ----------
    list_atoms : list of str
        list of PDB atoms
    atoms_only : bool
        fetch only atoms (ATOM)

    Returns
    -------
    residues : list ocif2pdb.convertf str
        e.g. "2", "3" etc.
    names : list of int
        e.g. "ALA", "ASN" etc.
    indices : list of int
        e.g. 0, 1 et.c
    """
    residues, names, indices = [], [], []
    for i, line in enumerate(list_atoms):
        if line.startswith(cif2pdb.convert.ATOM_ID) or \
                (not atoms_only and line.startswith(cif2pdb.convert.HETATM_ID)):
            residues.append(line[22:26].strip())
            names.append(line[17:20].strip())
            indices.append(i)
    return residues, names, indices


def transform_ranges(ranges):
    """
    Transform sequence of ranges into a list
    of sorted integers.

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
    if not ranges:
        return []

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
    Get atoms from list of PDB atoms (ATOM/HETATM flag)
    that belong to specified residues.

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


def get_number_of_residues_from_pdb_file(path):
    """
    Return number of residues from PDB file.

    Parameters
    ----------
    path : str
        path to PDB file

    Returns
    -------
    int
        Number of residues.
    """
    residues, _, _ = fetch_residues_from_pdb_file(path)
    number = get_number_of_residues(residues)
    return number


def get_number_of_residues(residues):
    """
    Return number of residues from residue list.

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
    if not residues:
        return 0

    prev = residues[0]
    number = 1
    for i in range(len(residues) - 1):
        cont = residues[i + 1]
        if prev != cont:
            number += 1
            prev = cont
    return number


def fetch_seqres_from_pdb_file(path):
    """
    Fetch SEQRES records from PDB file.

    Parameters
    ----------
    path : str
        path to PDB file

    Returns
    -------
    seqres : list of str
        e.g. "SEQRES   1 A  133  GLU ALA GLU ALA HIS ..."
    """
    with open(path, 'r') as f:
        data = f.readlines()
    seqres = [line for line in data
              if line.startswith("SEQRES")]
    return seqres


# TODO: generalize to include chains, residue ranges and PDB list
def get_sequence_from_pdb_file(path, seqres=True):
    """
    Return sequence from a PDB file.
    NOTE: only 3-letter aminoacid codes are supported!

    Parameters
    ----------
    path : str
        path to PDB file
    seqres : bool
        get sequence based on SEQRES information (if available)

    Returns
    -------
    seq : str
        Sequence of aminoacids e.g. AAKBCA etc.
    """

    if seqres:
        seq = []
        records = fetch_seqres_from_pdb_file(path)
        if records:
            for record in records:
                seq.extend([AA_SEQ_DICT[el] for el in record[19:].split()])
            seq = "".join(seq)
        else:
            # SEQRES information not found
            seqres = False
    if not seqres:
        residues, names, _ = fetch_residues_from_pdb_file(path,
                                                          atoms_only=True)
        seq = AA_SEQ_DICT[names[0]]
        prev = residues[0]
        for i in range(len(residues) - 1):
            cont = residues[i + 1]
            if prev != cont:
                seq += AA_SEQ_DICT[names[i + 1]]
                prev = cont
    return seq


@click.command()
@click.option('--input_path', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Full input path to the .pdb file.')
@click.option('--outdir', '-o', required=True, default=None,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Directory where the output file is (or will be) located.')
@click.option('--limits', '-l', required=False, default=None, type=str,
              help='Residue number limits (comma separated). '
                   'Examples: 100 or 40,100,200 or 10,100.')
@click.option('--file_root', '-f', required=True, default=None, type=str,
              help='Output file root.')
@click.option('--name_id', '-n', required=False, default=-1, type=int,
              help='Integer specifying which input path index to use '
                   'as a structure name indicator counting from the end. '
                   'E.g. for path/to/file.pdb we would get (default: -1): '
                   '"file.pdb" for -1, '
                   '"to" for -2, '
                   '"path" for -3.')
def _split_files_based_on_length(input_path, outdir, limits, file_root, name_id):
    """
    The script loads PDB file, calculates number of residues
    and saves results (format: 'name_id\\nlength') to file
    indicated by 'file_root' and 'limits' (if specified).
    When paralleled, It might be useful to split PDB files
    according to their length.

    Example usage (suppose file.pdb has 40 residues):

    python cif2pdb/resid.py  -i input_path/id_1/file.pdb
    -o output_path -l 40 -f split_len -n -2

    would produce file: 'output_path/split_len__geq_40'
    with one line appended: 'id_1 40\\n'.
    """
    residues, _, _ = fetch_residues_from_pdb_file(input_path)
    number = get_number_of_residues(residues)

    indicator_name = normpath(input_path).split(os.sep)[name_id]
    line = f"{indicator_name} {number}\n"

    if limits:
        limits = list(map(lambda x: int(x), limits.split(",")))

        if number < limits[0]:
            with open(join(outdir,
                           f"{file_root}__less_{limits[0]}"), 'a') as f:
                f.write(line)
        elif number >= limits[-1]:
            with open(join(outdir,
                           f"{file_root}__geq_{limits[-1]}"), 'a') as f:
                f.write(line)
        else:
            for i in range(len(limits[1:])):
                if limits[i] <= number < limits[i + 1]:
                    with open(join(outdir,
                                   f"{file_root}__geq_{limits[i]}"
                                   f"_less_{limits[i + 1]}"), 'a') as f:
                        f.write(line)
    else:
        with open(join(outdir, f"{file_root}"), 'a') as f:
            f.write(line)


if __name__ == "__main__":
    _split_files_based_on_length()
