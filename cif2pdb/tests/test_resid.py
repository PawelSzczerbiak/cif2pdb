import os
import glob
import pytest

from os.path import join
from click.testing import CliRunner

from cif2pdb.tests.utils import dict2str
from cif2pdb.tests.utils import get_expected_chains
from cif2pdb.tests.utils import generate_resid_file
from cif2pdb.tests.utils import compare_resid_files
from cif2pdb.resid import fetch_residues_from_pdb_file
from cif2pdb.resid import fetch_residues_from_pdb_list
from cif2pdb.resid import transform_ranges
from cif2pdb.resid import get_atoms_for_residues
from cif2pdb.resid import get_number_of_residues
from cif2pdb.resid import get_number_of_residues_from_pdb_file
from cif2pdb.resid import get_sequence_from_pdb_file
from cif2pdb.resid import _split_files_based_on_length

# Globals

FAKEPATH = join(os.getcwd(), "data/fake")
PDBPPATH = join(os.getcwd(), "data/pdb_expected")
EXPPATH = join(os.getcwd(), "data/resid_expected")
OUTPATH = join(os.getcwd(), "data/resid_generated")

FAKE_ATOM_LIST = \
    ['ATOM      1  N   MET 3   1     241.114 202.652 206.225  1.00 84.60      CA   N  \n',
     'ATOM      2  CA  MET 3   1     240.396 201.506 206.847  1.00 90.87      CA   C  \n',
     'ATOM      3  C   MET 3   1     241.214 200.870 207.974  1.00 75.98      CA   C  \n',
     'ATOM      4  O   MET 3   1     242.445 200.800 207.907  1.00 80.82      CA   O  \n',
     'TER \n',
     'HETATM    9  O   HOH 3   2     240.518 200.424 209.017  1.00 71.07      C    N  \n',
     'HETATM   10  O   HOH 3   2     241.163 199.814 210.179  1.00 56.20      C    C  \n',
     'HETATM   11  O   HOH 3   2     241.501 198.342 209.943  1.00 66.97      C    C  \n',
     'TER \n',
     'ATOM     18  N   ARG 3   3     242.789 198.014 209.975  1.00 44.72      CA   N  \n',
     'ATOM     19  CA  ARG 3   3     243.193 196.632 209.774  1.00 51.14      CA   C  \n',
     'ATOM     20  C   ARG 3   3     242.976 195.831 211.060  1.00 43.92      CA   C  \n',
     'ATOM     21  O   ARG 3   3     243.117 196.354 212.168  1.00 50.71      CA   O  \n',
     'ATOM     22  CB  ARG 3   3     244.654 196.546 209.322  1.00 32.77      CA   C  \n']

runner = CliRunner()


@pytest.fixture(scope="session", autouse=True)
def clean_generated_files():
    print("\nRemoving old generated files...")
    for f in glob.glob(join(OUTPATH, '*')):
        os.remove(f)
    assert glob.glob(join(OUTPATH, '*')) == []


# =========================================================
# Testing function behaviour of the cif2pdb/resid.py script
# =========================================================

def test_fetch_residues_from_pdb_file_not_exists():
    with pytest.raises(FileNotFoundError):
        fetch_residues_from_pdb_file("some_wrong_path")


def test_fetch_residues_from_pdb_file():
    # Testing on fake PDB file with a few ATOM and HETATM
    fakefile = join(FAKEPATH, "simp.pdb")
    expected = (['277', '277', '277', '277', '277', '277', '277', '277', '277', '277', '277',
                 '278', '278', '278', '278', '278', '278', '278', '278', '278', '278', '278',
                 '279', '279', '279', '279', '279', '279', '280', '280', '280', '280', '280',
                 '280', '280', '280', '281', '281', '281', '281', '281', '281', '281', '282',
                 '282', '282', '282', '282', '282', '282', '282', '282', '282', '282', '285',
                 '286', '287', '288', '289', '290', '291', '292'],
                ['ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG',
                 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG',
                 'SER', 'SER', 'SER', 'SER', 'SER', 'SER', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE',
                 'ILE', 'ILE', 'ILE', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'ARG',
                 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'HOH',
                 'HOH', 'HOH', 'HOH', 'HOH', 'HOH', 'HOH', 'HOH'],
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
                 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62])
    assert fetch_residues_from_pdb_file(fakefile) == expected, f" ERROR for file {fakefile}"


def test_fetch_residues_from_pdb_list_empty():
    assert fetch_residues_from_pdb_list([]) == ([], [], [])


def test_fetch_residues_from_pdb_list():
    input = FAKE_ATOM_LIST
    expected = (['1', '1', '1', '1', '2', '2', '2', '3', '3', '3', '3', '3'],
                ['MET', 'MET', 'MET', 'MET', 'HOH', 'HOH', 'HOH',
                 'ARG', 'ARG', 'ARG', 'ARG', 'ARG'],
                [0, 1, 2, 3, 5, 6, 7, 9, 10, 11, 12, 13])
    assert fetch_residues_from_pdb_list(input) == expected


def test_transform_ranges():
    inputs = ["", "100,1-5,2,300,200-202", "20,1-5",
              "2,1-5", "2,30,1-5", "2-5,40", "2,3,4",
              "2", "1-5,4-7", "200,1-5,4-7", "1-5,4-7,200",
              "1-5,2,4-7", "1-5,20,7-7"]
    expected = [[], [1, 2, 3, 4, 5, 100, 200, 201, 202, 300],
                [1, 2, 3, 4, 5, 20], [1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 30],
                [2, 3, 4, 5, 40], [2, 3, 4], [2], [1, 2, 3, 4, 5, 6, 7],
                [1, 2, 3, 4, 5, 6, 7, 200], [1, 2, 3, 4, 5, 6, 7, 200],
                [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 7, 20]]
    for inp, exp in zip(inputs, expected):
        assert transform_ranges(inp) == exp, f" ERROR transforming {inp}"


def test_get_atoms_for_residues_empty_both():
    assert get_atoms_for_residues([], []) == []


def test_get_atoms_for_residues_empty_first():
    assert get_atoms_for_residues([], [1, 2, 3]) == []


def test_get_atoms_for_residues_empty_second():
    assert get_atoms_for_residues(FAKE_ATOM_LIST, []) == []


def test_get_atoms_for_residues():
    input_list = FAKE_ATOM_LIST
    input_ranges = [2, 3]
    expected = list(FAKE_ATOM_LIST[i] for i in [5, 6, 7, 9, 10, 11, 12, 13])

    assert get_atoms_for_residues(input_list, input_ranges) == expected


def test_get_number_of_residues():
    inputs = ["", ["A"], ["1", "1", "2", "3", "3"],
              ["A", "B", "C", "B", "B", "A", "D", "D"],
              "aaacccdddaaa"]
    expected = [0, 1, 3, 6, 4]
    for inp, exp in zip(inputs, expected):
        assert get_number_of_residues(inp) == exp, f" ERROR transforming {inp}"


def test_get_sequence_from_pdb_file_not_exists():
    with pytest.raises(FileNotFoundError):
        get_sequence_from_pdb_file("some_wrong_path")


def test_get_sequence_from_pdb_file():
    pdb_files = glob.glob(join(PDBPPATH, "*.pdb"))
    assert pdb_files, f" ERROR directory with test data cannot be empty"
    for file in pdb_files:
        chains = get_expected_chains(file)
        # we expect only one-chain structures
        joined_chains = "".join(chains.values())
        assert get_sequence_from_pdb_file(file) == joined_chains, \
            f" ERROR for file {file}"


def test_get_sequence_from_pdb_file_seqres():
    pdb_file = join(FAKEPATH, "1AAQA.pdb")
    assert pdb_file, f" ERROR directory with test data cannot be empty"
    # Note: ATOM data are not complete!
    expected = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGG" \
               "IGGFIKVRQYDQIIIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
    assert get_sequence_from_pdb_file(pdb_file) == expected, \
        f" ERROR for file {pdb_file}"


def test_get_sequence_from_pdb_file_seqres_false():
    pdb_file = join(FAKEPATH, "1AAQA.pdb")
    assert pdb_file, f" ERROR directory with test data cannot be empty"
    # Note: ATOM data are not complete, but use that inseted of SEQRES
    expected = "PQITL"
    assert get_sequence_from_pdb_file(pdb_file, seqres=False) == expected, \
        f" ERROR for file {pdb_file}"


def test_get_number_of_residues_from_pdb_file_not_exists():
    with pytest.raises(FileNotFoundError):
        get_sequence_from_pdb_file("some_wrong_path")


def test_get_number_of_residues_from_pdb_file():
    pdb_files = glob.glob(join(PDBPPATH, "*.pdb"))
    assert pdb_files, f" ERROR directory with test data cannot be empty"
    for file in pdb_files:
        chains = get_expected_chains(file)
        # we expect only one-chain structures
        joined_chains = "".join(chains.values())
        assert get_number_of_residues_from_pdb_file(file) == len(joined_chains), \
            f" ERROR for file {file}"


# =============================================================
# Testing command line behaviour of the cif2pdb/resid.py script
# =============================================================

def test_help():
    response = runner.invoke(_split_files_based_on_length, ["--help"])
    assert response.exit_code == 0
    assert "The script loads PDB file, calculates number of residues" in response.output


def test_split_many_files():
    # Testing on some selected files only!
    files = ["2fjhL.pdb", "2fphX_1-77.pdb", "3d3bJ_3-87.pdb", "5uyl32.pdb"]
    file_root = "resid_many"
    params = {'-l': '47,85', '-f': file_root}
    # Generate/append 'resid' files
    for file in files:
        INPATH = join(PDBPPATH, file)
        generate_resid_file(dict2str(params), INPATH, OUTPATH)
    # Compare 'resid' files
    for expected in glob.glob(join(EXPPATH, f"{file_root}*")):
        generated = join(OUTPATH, os.path.basename(expected))
        compare_resid_files(generated, expected)


def test_split_one_file():
    # Testing -n flag on one file
    files = ["5uyl32.pdb"]
    file_root = "resid_one"
    params = {'-l': '47,85', '-f': file_root, '-n': -2}
    # Generate/append 'resid' files
    for file in files:
        INPATH = join(PDBPPATH, file)
        generate_resid_file(dict2str(params), INPATH, OUTPATH)
    # Compare 'resid' files
    for expected in glob.glob(join(EXPPATH, f"{file_root}*")):
        generated = join(OUTPATH, os.path.basename(expected))
        compare_resid_files(generated, expected)


def test_do_not_split():
    # Testing on some selected files only!
    files = ["2fjhL.pdb", "2fphX_1-77.pdb", "3d3bJ_3-87.pdb", "5uyl32.pdb"]
    file_root = "resid_not_split"
    params = {'-f': file_root}
    # Generate/append 'resid' files
    for file in files:
        INPATH = join(PDBPPATH, file)
        generate_resid_file(dict2str(params), INPATH, OUTPATH)
    # Compare 'resid' files
    for expected in glob.glob(join(EXPPATH, f"{file_root}*")):
        generated = join(OUTPATH, os.path.basename(expected))
        compare_resid_files(generated, expected)

