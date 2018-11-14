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
from qcfractal.testing import dask_server_fixture, recursive_dict_merge

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
    hooh = portal.data.get_molecule("hooh.json")
    mol_ret = client.add_molecules({"hooh": hooh})
    default_grid_spacing = 90

    # Geometric options
    torsiondrive_options = {
        "torsiondrive_meta": {
            "dihedrals": [[0, 1, 2, 3]],
            "grid_spacing": [default_grid_spacing]
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

    instance_options = copy.deepcopy(torsiondrive_options)
    instance_options["torsiondrive_meta"]["grid_spacing"] = [grid_spacing]

    # instance_options = {**instance_options, **keyword_augments}
    recursive_dict_merge(instance_options, keyword_augments)
    ret = client.add_service("torsiondrive", [mol_ret["hooh"]], instance_options)
    dask_server_fixture.await_services()
    assert len(dask_server_fixture.list_current_tasks()) == 0





