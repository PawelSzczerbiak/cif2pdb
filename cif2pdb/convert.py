# =====================================================
# The script extracts atoms from a CIF file, converts
# them into the PDB format and save in a .pdb file.
#
# Run:
#
#    python convert.py identifier ranges chain_type row_type indir outdir
#
# identifier - nameCHAINrest e.g. 1vcrAB or 1h8pA02
#              note: name contains exactly 4 chars;
#              chain(s) can contain any number of chars
#              rest may be e.g. domainID for CATH files
# ranges     - sequence of residue ranges e.g. 2,10-30
# chain_type - 1: single chain  e.g. 1vcrA
#              2: multiple chains e.g. 1vcrABC
#              3: single chain + non-zero rest
#                 e.g. 1h8pA02 (chain = A, rest = 02)
# row_type   - A: include atoms (ATOM)
#            - H: include heteroatoms (HETATM)
#            - AH: include both atoms & heteroatoms
# indir      - input dir where to look for the .cif
#              file e.g. indir/1vcr.cif
# outdir     - output dir where the .pdb file will
#              be saved e.g. outdir/1vcrA0.pdb
#
# Cif files may be download e.g. from:
# http://files.rcsb.org/download/ABCD.cif.gz
# =====================================================

import sys
from os.path import join

from cif2pdb.resid import transform_ranges
from cif2pdb.resid import get_atoms_for_residues


# Globals (used in CIF file)

LOOP_ID = "loop_"
LOOP_EL_ID = "_atom_site."
ATOM_ID = "ATOM"
HETATM_ID = "HETATM"

# Keys for dict that stores indices for specific fields in CIF file
KEY_RECORD = "_atom_site.group_PDB"
KEY_SERIAL = "_atom_site.id"
KEY_ATOM = "_atom_site.label_atom_id"
KEY_ALTLOC = "_atom_site.label_alt_id"
KEY_RES = "_atom_site.label_comp_id"
KEY_CHAIN = "_atom_site.auth_asym_id"  # AUTH
KEY_RESSEQ = "_atom_site.auth_seq_id"  # AUTH --- remove !!!!
KEY_ICODE = "_atom_site.pdbx_PDB_ins_code"
KEY_POS_X = "_atom_site.Cartn_x"
KEY_POS_Y = "_atom_site.Cartn_y"
KEY_POS_Z = "_atom_site.Cartn_z"
KEY_OCC = "_atom_site.occupancy"
KEY_TFACTOR = "_atom_site.B_iso_or_equiv"
KEY_SYMBOL = "_atom_site.type_symbol"
KEY_CHARGE = "_atom_site.pdbx_formal_charge"
KEY_MODEL_NUM = "_atom_site.pdbx_PDB_model_num"


def fetch_atoms_from_cif(name, chain, chain_type, row_type, indir):
    """
    Fetch atoms from cif file for speific chain.

    Parameters
    ----------
    name : str
        4-letter protein identification
    chain : str
        chain identification
    chain_type : str
        if 1: interpret chain as single chain
        if 2: interpret chain as multiple chains
    row_type : str
        if 'A': fetch ATOM rows
        if 'H': fetch HETATM rows
        if 'AH': fetch ATOM and HETATM rows

    Returns
    -------
    list of str
        lines with atoms details in CIF format
    dict (str : int)
        key: field name
        value: index where to look for specific field
    """

    assert chain, "{0} ERROR: chain not provided".format(identifier)
    assert chain_type in ["1", "2", "3"], \
        "{0} ERROR: chain type different from 1 or 2".format(chain_type)
    assert row_type in ["A", "H", "AH"], \
        "{0} ERROR: Row type different from A, H or AH".format(row_type)

    chains = []
    if chain_type in ["1", "3"]:
        chains = [chain]
    elif chain_type == "2":
        chains = list(chain)

    atoms = []
    with open(join(indir, "%s.cif" % name), 'r') as f:
        data = f.readlines()
        is_atom_loop = False
        is_atom_flag = False
        val_modelnum_first = None

        for i, line in enumerate(data):
            if not is_atom_loop:
                # if atom loop reached
                # create dictionary
                if line.startswith(LOOP_ID) and \
                        data[i + 1].startswith(LOOP_EL_ID):
                    fields = {}
                    num = 0
                    is_atom_loop = True
                continue

            # if atom flag not yet reached
            # fill the dictionary
            if not is_atom_flag:
                key = line.strip("\n").strip()
                fields[key] = num
                num += 1
                if data[i + 1].startswith(ATOM_ID) or \
                        data[i + 1].startswith(HETATM_ID):
                    is_atom_flag = True
                    assert KEY_CHAIN in fields, "{0}{1} ERROR: " \
                                                "chain key not found.".format(name, chain)
                    # index where to look for chain ID
                    idx_chain = fields[KEY_CHAIN]
                    # index where to look for model num (shall be 1)
                    idx_modelnum = fields[KEY_MODEL_NUM]
                    # value for the first encountered model num
                    val_modelnum_first = data[i + 1].split()[idx_modelnum]
                continue

            # go through all atoms and fetch atoms
            # belonging to the chain(s)
            line_start_atom = line.startswith(ATOM_ID)
            line_start_hetm = line.startswith(HETATM_ID)
            line_start_atom_hetm = line_start_atom or line_start_hetm
            if line_start_atom_hetm:
                if (row_type == "AH" and line_start_atom_hetm) or \
                        (row_type == "A" and line_start_atom) or \
                        (row_type == "H" and line_start_hetm):
                    line_splitted = line.split()
                    if line_splitted[idx_chain] in chains \
                            and line_splitted[idx_modelnum] == val_modelnum_first:
                        atoms.append(line.strip())
            else:
                break

        return atoms, fields


def create_pdb_atoms_from_cif(cif_atoms, cif_fields, identifier):
    """
    Transform cif atoms into pdb atoms.

    Parameters
    ----------
    cif_atoms : list of str
        lines with atoms details in CIF format
    cif_fields : dict (str : int)
        key: field name
        value: index where to look for specific field
    identifier : str
        nameCHAIN - used for logging

    Returns
    -------
    list of str
        lines with atoms details in PDB format
    """
    pdb_atoms = []

    # Indices
    idx_record = cif_fields[KEY_RECORD]
    # idx_serial = cif_fields[KEY_SERIAL]
    idx_atom = cif_fields[KEY_ATOM]
    idx_altloc = cif_fields[KEY_ALTLOC]
    idx_res = cif_fields[KEY_RES]
    idx_chain = cif_fields[KEY_CHAIN]
    idx_resseq = cif_fields[KEY_RESSEQ]
    idx_icode = cif_fields[KEY_ICODE]
    idx_pos_x = cif_fields[KEY_POS_X]
    idx_pos_y = cif_fields[KEY_POS_Y]
    idx_pos_z = cif_fields[KEY_POS_Z]
    idx_occ = cif_fields[KEY_OCC]
    idx_tfactor = cif_fields[KEY_TFACTOR]
    idx_symbol = cif_fields[KEY_SYMBOL]
    idx_charge = cif_fields[KEY_CHARGE]

    # Below we save only atoms belonging to the first
    # encountered conformation (if exists)
    first_conf = None  # first conformation
    is_first_conf = False  # whether first conformation is encountered

    for i, atom in enumerate(cif_atoms):
        elements = atom.split()
        assert len(elements) == len(cif_fields), \
            "{0} ERROR: wrong number of fields for atom at position {1}" \
                .format(identifier, i + 1)

        # Preprocessing

        if elements[idx_altloc] == '.':
            elements[idx_altloc] = ""
        else:
            if not is_first_conf:
                is_first_conf = True
                first_conf = elements[idx_altloc]
            else:
                if elements[idx_altloc] != first_conf:
                    continue

        if elements[idx_icode] == '?':
            elements[idx_icode] = ""
        if elements[idx_charge] == '?':
            elements[idx_charge] = ""
        elements[idx_atom] = elements[idx_atom] \
            .replace('\'', "").replace('\"', "")

        # Create and save line

        line = f"{elements[idx_record]:<6}" \
               f"{str(i + 1)[:5]:>5}" \
               f" " \
               f"{elements[idx_atom]:^4}" \
               f" " \
               f"{elements[idx_res]:>3}" \
               f" " \
               f"{elements[idx_chain][0]:>1}" \
               f"{elements[idx_resseq][-4:]:>4}" \
               f"{elements[idx_icode]:>1}" \
               f"   " \
               f"{elements[idx_pos_x][:7]:>8}" \
               f"{elements[idx_pos_y][:7]:>8}" \
               f"{elements[idx_pos_z][:7]:>8}" \
               f"{elements[idx_occ][:5]:>6}" \
               f"{elements[idx_tfactor][:5]:>6}" \
               f"          " \
               f"{elements[idx_symbol]:>2}" \
               f"{elements[idx_charge]:>1}" \
               f"\n"

        pdb_atoms.append(line)

    return pdb_atoms


if __name__ == "__main__":

    assert len(sys.argv) >= 7, \
        "{0} ERROR: wrong number of arguments".format(sys.argv[1:])

    identifier = sys.argv[1]  # nameCHAINrest
    ranges = sys.argv[2]  # sequence of ranges
    chain_type = sys.argv[3]  # chain type
    row_type = sys.argv[4]  # atoms / hetatms
    indir = sys.argv[5]  # full input path
    outdir = sys.argv[6]  # full output path
    if len(sys.argv) == 8:
        filename = sys.argv[7]  # alternative filename
    else:
        filename = identifier

    name = identifier[0:4]
    if chain_type == "3":
        chain = identifier[4]
    else:
        chain = identifier[4:]

    cif_atoms, cif_fields = fetch_atoms_from_cif(name, chain, chain_type, row_type, indir)

    assert cif_atoms, "{0} ERROR: nothing found for chain(s) {1}" \
        .format(identifier, chain)

    pdb_atoms = create_pdb_atoms_from_cif(cif_atoms, cif_fields, identifier)

    if ranges != "None":
        residues_to_fetch = transform_ranges(ranges)
        pdb_atoms = get_atoms_for_residues(pdb_atoms, residues_to_fetch)

    with open(join(outdir, "%s.pdb" % (filename)), "w") as f:
        for line in pdb_atoms:
            f.write(line)
        f.write("TER\n")
        f.write("END\n")

    print("{0} OK".format(filename))
