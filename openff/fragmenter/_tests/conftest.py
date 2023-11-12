import pytest
from openff.toolkit.topology import Molecule


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture()
def potanib() -> Molecule:
    return Molecule.from_smiles(
        "[H:40][c:3]1[c:8]([c:20]2[n:30][c:12]([c:14]([n:32]2[n:31][c:11]1[H:48])[C:2]#"
        "[C:1][c:13]3[c:9]([c:15]([c:4]([c:5]([c:16]3[C:26]([H:58])([H:59])[H:60])"
        "[H:42])[H:41])[C:21](=[O:36])[N:35]([H:66])[c:19]4[c:7]([c:6]([c:17]([c:18]"
        "([c:10]4[H:47])[C:29]([F:37])([F:38])[F:39])[C:28]([H:64])([H:65])[N:34]5"
        "[C:24]([C:22]([N:33]([C:23]([C:25]5([H:56])[H:57])([H:52])[H:53])[C:27]"
        "([H:61])([H:62])[H:63])([H:50])[H:51])([H:54])[H:55])[H:43])[H:44])[H:46])"
        "[H:49])[H:45]",
        allow_undefined_stereo=True,
    )


@pytest.fixture()
def abemaciclib() -> Molecule:
    return Molecule.from_smiles(
        "[H:38][c:1]1[c:2]([c:14]([n:28][c:5]([c:8]1[C:25]([H:64])([H:65])[N:33]2"
        "[C:17]([C:19]([N:34]([C:20]([C:18]2([H:46])[H:47])([H:50])[H:51])[C:26]"
        "([H:66])([H:67])[C:22]([H:55])([H:56])[H:57])([H:48])[H:49])([H:44])[H:45])"
        "[H:42])[N:35]([H:69])[c:16]3[n:29][c:6]([c:12]([c:13]([n:31]3)[c:7]4[c:3]"
        "([c:10]5[c:9]([c:11]([c:4]4[H:41])[F:36])[n:30][c:15]([n:32]5[C:27]([H:68])"
        "([C:23]([H:58])([H:59])[H:60])[C:24]([H:61])([H:62])[H:63])[C:21]([H:52])"
        "([H:53])[H:54])[H:40])[F:37])[H:43])[H:39]",
        allow_undefined_stereo=False,
    )


@pytest.fixture()
def dasatanib() -> Molecule:
    return Molecule.from_smiles(
        "[H:34][c:1]1[c:2]([c:6]([c:7]([c:8]([c:3]1[H:36])[Cl:33])[N:28]([H:57])[C:14]"
        "(=[O:30])[c:9]2[c:5]([n:23][c:13]([s:32]2)[N:29]([H:58])[c:11]3[c:4]([c:10]"
        "([n:24][c:12]([n:25]3)[C:20]([H:50])([H:51])[H:52])[N:26]4[C:15]([C:17]([N:27]"
        "([C:18]([C:16]4([H:41])[H:42])([H:45])[H:46])[C:21]([H:53])([H:54])[C:22]"
        "([H:55])([H:56])[O:31][H:59])([H:43])[H:44])([H:39])[H:40])[H:37])[H:38])"
        "[C:19]([H:47])([H:48])[H:49])[H:35]",
        allow_undefined_stereo=True,
    )


@pytest.fixture()
def butane() -> Molecule:
    return Molecule.from_smiles(
        "[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])"
        "[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]",
        allow_undefined_stereo=True,
    )
