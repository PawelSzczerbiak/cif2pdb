import click

from os.path import join

import cif2pdb.resid


# Globals (used in mmCIF file)

LOOP_ID = "loop_"
LOOP_EL_ID = "_atom_site."
ATOM_ID = "ATOM"
HETATM_ID = "HETATM"

# Dictionary keys that specify indices for fields in the CIF file
KEY_RECORD = "_atom_site.group_PDB"  # record name
KEY_SERIAL = "_atom_site.id"  # atom serial number
KEY_ATOM = "_atom_site.label_atom_id"  # atom name
KEY_ALTLOC = "_atom_site.label_alt_id"  # alternate location indicator
KEY_RES = "_atom_site.label_comp_id"  # residue name
KEY_CHAIN = "_atom_site.auth_asym_id"  # strand ID / chain ID (AUTH)
KEY_RESSEQ = "_atom_site.label_seq_id"  # residue sequence number
KEY_ICODE = "_atom_site.pdbx_PDB_ins_code"  # code for insertion of residues
KEY_POS_X = "_atom_site.Cartn_x"  # orthogonal coordinate for X in [A]
KEY_POS_Y = "_atom_site.Cartn_y"  # orthogonal coordinate for Y in [A]
KEY_POS_Z = "_atom_site.Cartn_z"  # orthogonal coordinate for X in [A]
KEY_OCC = "_atom_site.occupancy"  # occupancy
KEY_TFACTOR = "_atom_site.B_iso_or_equiv"  # temperature factor
KEY_SYMBOL = "_atom_site.type_symbol"  # element symbol, right-justified
KEY_CHARGE = "_atom_site.pdbx_formal_charge"  # charge on the atom
KEY_MODEL_NUM = "_atom_site.pdbx_PDB_model_num"  # model number


def _fetch_atoms_from_cif(name, chain, chain_type, row_type, indir):
    """
    Fetch atoms from mmCIF file for specific chain.

    Parameters
    ----------
    name : str
        4-letter protein identification
    chain : str
        chain identification
    chain_type : str
        if 1 or 3: interpret chain as single chain
        if 2: interpret chain as multiple 1-letter chains
    row_type : str
        if 'A': fetch ATOM rows
        if 'H': fetch HETATM rows
        if 'AH': fetch ATOM and HETATM rows

    Returns
    -------
    list of str
        lines with atoms details in mmCIF format
    dict (str : int)
        key: field name
        value: index where to look for specific field
    """

    assert chain, "{0} ERROR: chain not provided".format(name)
    assert chain_type in ["1", "2", "3"], \
        "{0} ERROR: chain type different from 1, 2 or 3".format(chain_type)
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


def _create_pdb_atoms_from_cif(cif_atoms, cif_fields, identifier):
    """
    Transform mmCIF atoms into pdb atoms.

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
               f"{elements[idx_pos_x][:8]:>8}" \
               f"{elements[idx_pos_y][:8]:>8}" \
               f"{elements[idx_pos_z][:8]:>8}" \
               f"{elements[idx_occ][:6]:>6}" \
               f"{elements[idx_tfactor][:6]:>6}" \
               f"          " \
               f"{elements[idx_symbol]:>2}" \
               f"{elements[idx_charge]:>1}" \
               f"\n"

        pdb_atoms.append(line)

    return pdb_atoms


# TODO Add possibility to fetch cif file if not present in the input directory

@click.command()
@click.option('--identifier', '-d', required=True, default=None, type=str,
              help='mmCIF file identifier in the form: nameCHAINrest. '
                   'Examples: 1vbzA, 1vcrAB, 1h8pA02. '
                   'Notes: "name" (PDB ID) contains exactly 4 chars; '
                   '"CHAIN" can contain any number of chars; '
                   '"rest" might be e.g. domainID for CATH files or be empty.')
@click.option('--ranges', '-s', required=False, default=None, type=str,
              help='Sequence of residue ranges separated by comma or dash. '
                   'Examples: 1,2,3,4 or 1-10 or 2,5,10-30.')
@click.option('--chain_type', '-c', required=False, default='1', type=str,
              help='Type of chain(s) specified in the identifier. '
                   '1: single chain e.g. A in 1vcrA or A1 in 1vcrA1, '
                   '2: multiple 1-letter chains e.g. A, B, C in 1vcrABC, '
                   '3: single 1-letter chain + non-zero rest e.g. A in 1h8pA02.')
@click.option('--row_type', '-r', required=False, default='AH', type=str,
              help='Type of rows to be fetched. '
                   'A: include atoms (ATOM), '
                   'H: include heteroatoms (HETATM), '
                   'AH: include both atoms abd heteroatoms.')
@click.option('--indir', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Input directory where the name.cif file is located. '
                   'Note: "name" is the 4-letter PDB ID specified in the identifier.')
@click.option('--outdir', '-o', required=True, default=None,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output directory where the .pdb file will saved.')
@click.option('--outfile', '-f', required=False, default=None, type=str,
              help='Alternative name of the output .pdb file.')
def _convert_cif_to_pdb(identifier, ranges, chain_type, row_type, indir, outdir, outfile):
    """
    The script extracts atoms from mmCIF file, converts
    them into the PDB format and save in a .pdb file.
    mmCIF files may be download e.g. from:
    http://files.rcsb.org/download/abcd.cif.gz

    Example usage (note: 'abcdA.cif' shall be present in 'input_directory'):

    python cif2pdb/convert.py -d abcdA -i input_directory/ -o output_directory/
    """

    name = identifier[0:4]
    if chain_type == "3":
        chain = identifier[4]
    else:
        chain = identifier[4:]

    cif_atoms, cif_fields = _fetch_atoms_from_cif(name, chain, chain_type, row_type, indir)

    assert cif_atoms, "{0} ERROR: nothing found for chain(s) {1}" \
        .format(identifier, chain)

    pdb_atoms = _create_pdb_atoms_from_cif(cif_atoms, cif_fields, identifier)

    if ranges:
        residues_to_fetch = cif2pdb.resid.transform_ranges(ranges)
        pdb_atoms = cif2pdb.resid.get_atoms_for_residues(pdb_atoms, residues_to_fetch)

    # Save .pdb file

    if not outfile:
        outfile = identifier

    with open(join(outdir, "%s.pdb" % (outfile)), "w") as f:
        for line in pdb_atoms:
            f.write(line)
        f.write("TER\n")
        f.write("END\n")

    print("{0} OK".format(outfile))


if __name__ == "__main__":
    _convert_cif_to_pdb()
