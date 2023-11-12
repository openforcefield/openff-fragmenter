"""
Regression tests for the WBO and Pfizer fragmentation schemes.
"""
import pytest
from openff.fragmenter.fragment import WBOFragmenter
from openff.toolkit.topology import Molecule


@pytest.mark.parametrize(
    "parent, fragments",
    [
        pytest.param(
            "[H:30][c:6]1[c:7]([c:8]([c:16]([c:4]([c:5]1[H:29])[C@@:3]2([C:22](=[O:23])[N:20]([C:18](=[N+:17]2[H:37])[N:19]([H:38])[H:42])[C:21]([H:39])([H:40])[H:41])[C:2]([H:27])([H:28])[C:1]([H:24])([H:25])[H:26])[H:36])[c:9]3[c:10]([c:11]([c:12]([c:13]([c:15]3[H:35])[Cl:14])[H:34])[H:33])[H:32])[H:31]",
            (
                "[H:14][c:1]1[c:2]([c:5]([c:10]([c:6]([c:3]1[H:16])[H:19])[c:11]2[c:7]([c:4]([c:8]([c:12]([c:9]2[H:22])[Cl:13])[H:21])[H:17])[H:20])[H:18])[H:15]",
                "[H:17][c:1]1[c:2]([c:4]([c:6]([c:5]([c:3]1[H:19])[H:21])[C@@:9]2([C:7](=[O:16])[N:13]([C:8](=[N+:14]2[H:30])[N:15]([H:31])[H:32])[C:11]([H:25])([H:26])[H:27])[C:12]([H:28])([H:29])[C:10]([H:22])([H:23])[H:24])[H:20])[H:18]",
            ),
        ),
        pytest.param(
            "[H:30][c:6]1[c:7]([c:8]([c:16]([c:4]([c:5]1[H:29])[C@@:3]2([C:22](=[O:23])[N:20]([C:18](=[N+:17]2[H:37])[N:19]([H:38])[H:42])[C:21]([H:39])([H:40])[H:41])[C:2]([H:27])([H:28])[C:1]([H:24])([H:25])[H:26])[H:36])[c:9]3[c:10]([c:11]([c:12]([c:13]([c:15]3[H:35])[Cl:14])[H:34])[H:33])[H:32])[H:31]",
            (
                "[H:14][c:1]1[c:2]([c:5]([c:10]([c:6]([c:3]1[H:16])[H:19])[c:11]2[c:7]([c:4]([c:8]([c:12]([c:9]2[H:22])[Cl:13])[H:21])[H:17])[H:20])[H:18])[H:15]",
                "[H:17][c:1]1[c:2]([c:4]([c:6]([c:5]([c:3]1[H:19])[H:21])[C@@:9]2([C:7](=[O:16])[N:13]([C:8](=[N+:14]2[H:30])[N:15]([H:31])[H:32])[C:11]([H:25])([H:26])[H:27])[C:12]([H:28])([H:29])[C:10]([H:22])([H:23])[H:24])[H:20])[H:18]",
            ),
        ),
        pytest.param(
            "[H:30][c:6]1[c:7]([c:8]([c:16]([c:4]([c:5]1[H:29])[C@@:3]2([C:22](=[O:23])[N:20]([C:18](=[N+:17]2[H:37])[N:19]([H:38])[H:42])[C:21]([H:39])([H:40])[H:41])[C:2]([H:27])([H:28])[C:1]([H:24])([H:25])[H:26])[H:36])[c:9]3[c:10]([c:11]([c:12]([c:13]([c:15]3[H:35])[Cl:14])[H:34])[H:33])[H:32])[H:31]",
            (
                "[H:14][c:1]1[c:2]([c:5]([c:10]([c:6]([c:3]1[H:16])[H:19])[c:11]2[c:7]([c:4]([c:8]([c:12]([c:9]2[H:22])[Cl:13])[H:21])[H:17])[H:20])[H:18])[H:15]",
                "[H:17][c:1]1[c:2]([c:4]([c:6]([c:5]([c:3]1[H:19])[H:21])[C@@:9]2([C:7](=[O:16])[N:13]([C:8](=[N+:14]2[H:30])[N:15]([H:31])[H:32])[C:11]([H:25])([H:26])[H:27])[C:12]([H:28])([H:29])[C:10]([H:22])([H:23])[H:24])[H:20])[H:18]",
            ),
        ),
    ],
)
@pytest.mark.slow
def test_wbo_regression(parent, fragments):
    """Make sure we can produce the expected fragments using both openeye and AT+RDKIT."""
    parent = Molecule.from_mapped_smiles(parent)
    engine = WBOFragmenter(keep_non_rotor_ring_substituents=True, threshold=0.03)
    result = engine.fragment(molecule=parent)
    result_fragments = [fragment.molecule for fragment in result.fragments]
    assert len(result_fragments) == 3
    for fragment in fragments:
        fragment_mol = Molecule.from_mapped_smiles(fragment)
        assert fragment_mol in result_fragments


def test_stereo_error():
    """
    Make sure the molecule can be fragmented, before PR123 this molecule would result in a runtime error as new
    stereocenters would be created.
    """
    parent = Molecule.from_smiles(
        "[H][O][C]([H])([H])[C]([H])([H])[C]([H])([N]([H])[C]1=[N][C]([H])=[C]2[C](=[N]1)[N]([C]([H])([H])[H])[C](=[O])[C]([O][C]1=[C]([F])[C]([H])=[C]([F])[C]([H])=[C]1[H])=[C]2[H])[C]([H])([H])[C]([H])([H])[O][H]"
    )
    engine = WBOFragmenter()
    _ = engine.fragment(parent)
