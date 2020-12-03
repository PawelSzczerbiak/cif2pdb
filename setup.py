from setuptools import setup, find_packages

setup(
    name="cif2pdb",
    packages=find_packages(),
    install_requires=['pytest',
                      'click'],
    python_requires='>=3.6',
    author="Paweł Szczerbiak",
    author_email="pawel.szczerbiak@uj.edu.pl",
    description="Convert ATOM/HETATM information from "
                "mmCIF file into PDB format",
    license='BSD-3-Clause',
    url="https://github.com/PawelSzczerbiak/cif2pdb",
)
