""" Test launch fragmenter """

import pytest
import qcfractal.interface as portal
import warnings
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
    workflow_id = 'wbo'
    workflow_json = get_fn('workflows.json')

    workflow = workflow_api.WorkFlow(client=client, workflow_json=workflow_json, workflow_id=workflow_id)
    workflow.workflow(molecules_smiles=smiles, molecule_titles=['pentane'])

    assert len(workflow.qcfractal_jobs) == 1
    key = list(workflow.qcfractal_jobs.keys())[0]
    assert len(workflow.qcfractal_jobs[key]) == 2
    assert len(workflow.qcfractal_jobs[key]['torsiondrive_input']) == 3
    assert len(workflow.qcfractal_jobs[key]['torsiondrive_input']['[CH3:1][CH2:4][CH2:3][CH3:2]']['initial_molecule']) == 1

@testing.using_rdkit
@testing.using_geometric
@testing.using_torsiondrive
def test_check_workflow_options(fractal_compute_server):
    """Check that the workflow options are the same in the json and database"""

    client = portal.FractalClient(fractal_compute_server)
    workflow_id = 'combinatorial'
    workflow_json_f = get_fn('workflows.json')
    with open(workflow_json_f) as file:
        workflow_json = json.load(file)[workflow_id]['fragmenter']

    workflow = workflow_api.WorkFlow(client=client, workflow_json=workflow_json_f, workflow_id=workflow_id)
    #workflow = portal.collections.OpenFFWorkflow(workflow_id, client, **workflow_json)
    assert workflow_api._check_workflow(workflow_json, workflow.off_workflow)

    # change some items in workflow_json
    workflow_json['enumerate_states']['options']['protonation'] = True
    print(workflow_json['enumerate_states']['options']['protonation'])

    with pytest.raises(ValueError):
        assert workflow_api._check_workflow(workflow_json, workflow.off_workflow)
    with open(workflow_json_f) as file:
            workflow_json = json.load(file)[workflow_id]['fragmenter']
    workflow_json['enumerate_fragments']['scheme'] = 'wbo'
    with pytest.raises(ValueError):
        assert workflow_api._check_workflow(workflow_json, workflow.off_workflow)


@testing.using_rdkit
@testing.using_geometric
@testing.using_torsiondrive
def test_versions():
    """Chack that the versions given in JSON is the version of the software used"""
    pass

def test_enumerate_states(fractal_compute_server):
    """Test enumerating states"""
    client = portal.FractalClient(fractal_compute_server)
    worfklow_json = get_fn('workflows.json')
    wf = workflow_api.WorkFlow(workflow_id='combinatorial', client=client)
    states = wf.enumerate_states('CCCC')
    keys = list(states.keys())
    assert 'provenance' in keys
    assert 'states' in keys
    assert len(states['states']) == 1
    assert list(states['states'])[0] == '[H]C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]'


def test_enumerate_fragments_combinatorial(fractal_compute_server):
    client = portal.FractalClient(fractal_compute_server)
    wf = workflow_api.WorkFlow(workflow_id='combinatorial', client=client)

    frags =  wf.enumerate_fragments('CCCC')
    assert not frags

    frags = wf.enumerate_fragments('CCCCC')
    assert len(frags) == 1

def test_enumerate_fragments_wbo(fractal_compute_server):
    client = portal.FractalClient(fractal_compute_server)
    workflow_json = get_fn('workflows.json')
    wf = workflow_api.WorkFlow(workflow_id='wbo', workflow_json=workflow_json, client=client)

    frags = wf.enumerate_fragments('CCCCC')
    assert 'CCCC' in frags
    assert len(frags) == 1

def test_generate_torsiondrive_input():
    pass

def test_provenance():
    pass






