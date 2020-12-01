import os

from os.path import join
from click.testing import CliRunner

from cif2pdb.tests.utils import dict2str
from cif2pdb.tests.utils import generate_pdb_file
from cif2pdb.tests.utils import compare_pdb_files
from cif2pdb.convert import _convert_cif_to_pdb

# ================================================================
# Testing command-line behaviour of the cif2pdb/convert.py script
# ================================================================

# Different databases encode their PDB IDs in a different way
# especially, regarding chain and residue range (below).
# CATH:
# - chain: AUTHOR-defined
# - residue number: PDB-defined
# PDB90 etc.:
# - chain: AUTHOR-defined
# - residue number: all residues
# SCOP:
# - chain: AUTHOR-defined
# - residue number: AUTHOR-defined

# Default input and output paths
INPATH = join(os.getcwd(), "data/cif")
OUTPATH = join(os.getcwd(), "data/pdb_generated")
EXPPATH = join(os.getcwd(), "data/pdb_expected")

runner = CliRunner()


def test_help():
    response = runner.invoke(_convert_cif_to_pdb, ["--help"])
    assert response.exit_code == 0
    assert "The script extracts atoms from a CIF file" in response.output


def test_convert_cath_2fphX01():
    # AUTHORS chain diff from PDB
    # PDB residues overlap with AUTHORS
    params = { '-d': '2fphX01',
               '-s': '1-77',
               '-c': '3',
               '-r': 'A'}
    generate_pdb_file(dict2str(params), INPATH, OUTPATH)
    compare_pdb_files(join(OUTPATH, '2fphX01.pdb'),
                      join(EXPPATH, '2fphX_1-77.pdb'))


def test_convert_cath_d3bJ00():
    # AUTHORS chain diff from PDB
    # some PDB residues unique
    params = { '-d': '3d3bJ00',
               '-s': '3-87',
               '-c': '3',
               '-r': 'A'}
    generate_pdb_file(dict2str(params), INPATH, OUTPATH)
    compare_pdb_files(join(OUTPATH, '3d3bJ00.pdb'),
                      join(EXPPATH, '3d3bJ_3-87.pdb'))


def test_convert_pdb_2fjhL():
    # AUTHORS chain diff from PDB
    # all residues
    params = { '-d': '2fjhL',
               '-c': '1',
               '-r': 'A'}
    generate_pdb_file(dict2str(params), INPATH, OUTPATH)
    compare_pdb_files(join(OUTPATH, '2fjhL.pdb'),
                      join(EXPPATH, '2fjhL.pdb'))


def test_convert_pdb_5uyl32():
    # multi-letter AUTHORS chain
    # all residues
    params = {'-d': '5uyl32',
              '-c': '1',
              '-r': 'A'}
    generate_pdb_file(dict2str(params), INPATH, OUTPATH)
    compare_pdb_files(join(OUTPATH, '5uyl32.pdb'),
                      join(EXPPATH, '5uyl32.pdb'))

# TODO
# SCOP
# - many chains (-c 2)
# - -r H
# - -r AH
# some problematic...
