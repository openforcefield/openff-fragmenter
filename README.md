# Fragmenter

[![Test Status](https://github.com/openforcefield/fragmenter/actions/workflows/ci.yaml/badge.svg?branch=main)](https://github.com/openforcefield/fragmenter/actions/workflows/ci.yaml)
[![Documentation Status](https://readthedocs.org/projects/fragmenter/badge/?version=latest)](https://fragmenter.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/openforcefield/openff-fragmenter/branch/main/graph/badge.svg)](https://codecov.io/gh/openforcefield/openff-fragmenter/branch/main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Software DOI](https://img.shields.io/badge/Code%20DOI-zenodo.127185286-blue)](https://zenodo.org/badge/latestdoi/127185286)
[![Paper DOI](https://img.shields.io/badge/Paper%20DOI-10.1101%2F2020.08.27.270934-blue)](https://doi.org/10.1101/2020.08.27.270934)

A package for fragmenting molecules for quantum mechanics torsion scans.

More information about using this package and its features can be found in the [documentation](
https://fragmenter.readthedocs.io/en/latest/).

**Warning:** This code is currently experimental and under active development. If you are using this code,
be aware that it is not guaranteed to provide the correct results, and the API can change without notice.

## Installation

The package and its dependencies can be installed using the `conda` package manager:

```shell
conda install -c conda-forge openff-fragmenter
```

## Getting Started

*We recommend viewing the getting started example in a Jupyter notebook. [This full example can be found here](
https://github.com/openforcefield/fragmenter/blob/master/examples/fragment-molecules.ipynb)*. 

Here will will show how a drug-like molecule can be fragmented using this framework, and how those fragments can 
be easily visualised using its built-in helper utilities.

To begin with we load in the molecule to be fragmented. Here we load Cobimetinib directly using its SMILES 
representation using the [Open Force Field toolkit](https://github.com/openforcefield/openff-toolkit):

```python
from openff.toolkit.topology import Molecule

parent_molecule = Molecule.from_smiles(
    "OC1(CN(C1)C(=O)C1=C(NC2=C(F)C=C(I)C=C2)C(F)=C(F)C=C1)[C@@H]1CCCCN1"
)
```

Next we create the fragmentation engine which will perform the actual fragmentation. Here we will use the recommended 
`WBOFragmenter` with default options:

```python
from openff.fragmenter.fragment import WBOFragmenter

frag_engine = WBOFragmenter()
# Export the engine's settings directly to JSON
frag_engine.json()
```

Use the engine to fragment the molecule:

```python
result = frag_engine.fragment(parent_molecule)
# Export the result directly to JSON
result.json()
```

Any generated fragments will be returned in a ``FragmentationResult`` object. We can loop over each of the generated 
fragments and print both the SMILES representation of the fragment as well as the map indices of the bond that the
fragment was built around:

```python
for fragment in result.fragments:
    print(f"{fragment.bond_indices}: {fragment.smiles}")
```

Finally, we can visualize the produced fragments:

```python
from openff.fragmenter.depiction import depict_fragmentation_result

depict_fragmentation_result(result=result, output_file="example_fragments.html")
```

### Copyright

Copyright (c) 2018, Chaya D. Stern

#### Acknowledgements

Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
