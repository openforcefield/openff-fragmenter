""" Tests generating files for qm torsion scan """

import pytest
import qcfractal.interface as portal
from qcfractal import testing
#from qcfractal.testing import fractal_compute_server

# @testing.using_rdkit
# @testing.using_geometric
# @testing.using_torsiondrive
# def test_torsiondrive_run(fractal_compute_server):
#
#     # Cannot use this fixture without these services. Also cannot use `mark` and `fixture` decorators
#     pytest.importorskip("torsiondrive")
#     pytest.importorskip("geometric")
#     pytest.importorskip("rdkit")
#
#     client = portal.FractalClient(fractal_compute_server)
#
#     # Add a HOOH
#     hooh = {
#         'symbols': ['H', 'O', 'O', 'H'],
#         'geometry': [
#              1.84719633,  1.47046223,  0.80987166,
#              1.3126021,  -0.13023157, -0.0513322,
#             -1.31320906,  0.13130216, -0.05020593,
#             -1.83756335, -1.48745318,  0.80161212
#         ],
#         'name': 'HOOH',
#         'connectivity': [[0, 1, 1], [1, 2, 1], [2, 3, 1]],
#     }
#
#     # Geometric options
#     tdinput = {
#     "initial_molecule": [hooh],
#     "keywords": {
#         "dihedrals": [[0, 1, 2, 3]],
#         "grid_spacing": [90]
#     },
#     "optimization_spec": {
#         "program": "geometric",
#         "keywords": {
#             "coordsys": "tric",
#         }
#     },
#     "qc_spec": {
#         "driver": "gradient",
#         "method": "UFF",
#         "basis": None,
#         "keywords": None,
#         "program": "rdkit",
#     },
# }
#     ret = client.add_service([tdinput])
#     fractal_compute_server.await_results()
#     assert len(fractal_compute_server.list_current_tasks()) == 0
#




