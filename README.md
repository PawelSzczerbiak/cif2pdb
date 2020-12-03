# cif2pdb

## Introduction
The package:
- converts `mmCIF` file into `PDB` structure file format (extracting `ATOM` / `HETATM` information)
- calculates sequence length (number of residues) from a `PDB` file and saves results based on specified residue limits
- extracts list of residues, sequence, etc. from a `PDB` file

## Installation

### Requirements

- `Click`
- `PyTest`
- `Biopython`

### Conda

In main `cif2pdb` directory run e.g.

```
conda create --name cif2pdb
conda activate cif2pdb
python setup.py install
```

### Virtualenv

In main `cif2pdb` directory run e.g.

```
virtualenv --python=/usr/bin/python3.6 cif2pdb_env
source cif2pdb_env/bin/activate
pip install .
```

## Examples

### `mmCIF` to `PDB` convertion

In order to fetch atoms for `L` chain residues in the range `10-30` from `mmCIF` file located in `cif2pdb/tests/data/cif/2fjh.cif` and save results to `cif2pdb/tests/data/pdb_generated/result.pdb` run:

`python cif2pdb/convert.py -d 2fjhL -s 10-30 -c 1 -r A -i cif2pdb/tests/data/cif/ -o cif2pdb/tests/data/pdb_generated/ -f result`

Detailed help:  

`python cif2pdb/convert.py --help`

### Splitting based on sequence length

Suppose we have a few `PDB` files located in `cif2pdb/tests/data/pdb_expected/` and we want to specify which of them are shorter than 40 residues. To this purpose, on each file (below for `5upl32.pdb`) we can run:  

`python cif2pdb/resid.py -i cif2pdb/tests/data/pdb_expected/5uyl32.pdb -o cif2pdb/tests/data/resid_generated/ -l 40 -f resid_len`

which creates two files in `cif2pdb/tests/data/resid_generated/` directory: 
- `resid_len__less_40`
- `resid_len__geq_40` 

with the lines of the form `{filename}.pdb length\n`.

The above command might be easily paralelled e.g. by using `xargs`. 

We can also specify the `-n` flag that replaces the default `{filename}.pdb` by the name of some previous directories in the path. It is particularily usefull when we have the following tree:
```
SOME/PATH/
  STRUCT_1/
    model_1.pdb
    model_2.pdb
    model_3.pdb
  STRUCT_2/
    model_1.pdb
    model_2.pdb
```
and we want to split the `model_1.pdb` files based on their sequence length (shorter / longer or equal than 100 residues) and differ between structures based on the preceeding directory name (`STRUCT_1` etc.). It is enough to run:

`python cif2pdb/resid.py -i SOME/PATH/STRUCT_1/model_1.pdb -o SOME/PATH/split_results/ -l 100 -f resid_len -n -2`

If we don't specify the `-l` flag we can produce one output file with all structures and their corresponding sequence lengths. 

Detailed help:  

`python cif2pdb/resid.py --help`
