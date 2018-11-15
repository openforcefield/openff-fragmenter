""" Tests generating files for qm torsion scan """

import unittest
import json
from fragmenter.tests.utils import get_fn, has_openeye
import fragmenter.torsions as torsions
from fragmenter import utils, chemi
from cmiles import to_canonical_smiles_oe
import warnings


import copy
import pytest
import qcfractal.interface as portal
from qcfractal import testing
from qcfractal.testing import dask_server_fixture

@testing.using_rdkit
@testing.using_geometric
@testing.using_torsiondrive
def test_torsiondrive_run(dask_server_fixture):

    # Cannot use this fixture without these services. Also cannot use `mark` and `fixture` decorators
    pytest.importorskip("torsiondrive")
    pytest.importorskip("geometric")
    pytest.importorskip("rdkit")

    client = portal.FractalClient(dask_server_fixture.get_address())

    # Add a HOOH
    hooh = {
        'symbols': ['H', 'O', 'O', 'H'],
        'geometry': [
             1.84719633,  1.47046223,  0.80987166,
             1.3126021,  -0.13023157, -0.0513322,
            -1.31320906,  0.13130216, -0.05020593,
            -1.83756335, -1.48745318,  0.80161212
        ],
        'name': 'HOOH',
        'connectivity': [[0, 1, 1], [1, 2, 1], [2, 3, 1]],
    }
    mol_ret = client.add_molecules({"hooh": hooh})

    # Geometric options
    instance_options = {
        "torsiondrive_meta": {
            "dihedrals": [[0, 1, 2, 3]],
            "grid_spacing": [90]
        },
        "optimization_meta": {
            "program": "geometric",
            "coordsys": "tric",
        },
        "qc_meta": {
            "driver": "gradient",
            "method": "UFF",
            "basis": "",
            "options": "none",
            "program": "rdkit",
        },
    }

    ret = client.add_service("torsiondrive", [mol_ret["hooh"]], instance_options)
    dask_server_fixture.await_services()
    assert len(dask_server_fixture.list_current_tasks()) == 0





