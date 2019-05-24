""" Test launch fragmenter """

import pytest
import qcfractal.interface as portal
import fragmenter
from fragmenter import workflow_api, chemi
from fragmenter.tests.utils import get_fn, has_crank, has_openeye
import json
import copy
from qcfractal import testing
from qcfractal.testing import fractal_compute_server
qcfractal = pytest.importorskip("qcfractal")

@testing.using_rdkit
@testing.using_geometric
@testing.using_torsiondrive
def test_workflow(fractal_compute_server):
    """Fragmenter regression test"""

    pytest.importorskip("torsiondrive")
    pytest.importorskip("geometric")
    pytest.importorskip("rdkit")

    client = portal.FractalClient(fractal_compute_server)
    smiles = 'CCCCC'
    workflow_id = 'example'
    workflow_json = get_fn('workflows.json')

    workflow = workflow_api.WorkFlow(client=client, workflow_json=workflow_json, workflow_id=workflow_id)
    workflow.workflow(molecules_smiles=smiles, molecule_titles=['pentane'])

    assert len(workflow.qcfractal_jobs) == 1
    key = list(workflow.qcfractal_jobs.keys())[0]
    assert len(workflow.qcfractal_jobs[key]) == 3
    assert len(workflow.qcfractal_jobs[key]['optimization_input']) == 0
    assert len(workflow.qcfractal_jobs[key]['torsiondrive_input']) == 3
    assert len(workflow.qcfractal_jobs[key]['torsiondrive_input']['(0, 2, 3, 1)']['initial_molecule']) == 8

def test_check_workflow_options():
    """Check that the workflow options are the same in the json and database"""
    pass

def test_versions():
    """Chack that the versions given in JSON is the version of the software used"""
    pass

def test_enumerate_states():
    pass

def test_enumerate_fragments_combinatorial():
    pass

def test_enumerate_fragments_wbo():
    pass

def test_generate_torsiondrive_input():
    pass






