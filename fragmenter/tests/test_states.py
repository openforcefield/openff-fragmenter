import pytest
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import OpenEyeToolkitWrapper, RDKitToolkitWrapper

from fragmenter.states import _enumerate_stereoisomers
from fragmenter.tests.utils import global_toolkit_wrapper


@pytest.mark.parametrize(
    "smiles, force_flip, expected",
    [
        (
            "FC(Br)(Cl)[C@@H](Br)(Cl)",
            True,
            [
                "F[C@](Br)(Cl)[C@@H](Br)(Cl)",
                "F[C@](Br)(Cl)[C@H](Br)(Cl)",
                "F[C@@](Br)(Cl)[C@@H](Br)(Cl)",
                "F[C@@](Br)(Cl)[C@H](Br)(Cl)",
            ],
        ),
        (
            "FC(Br)(Cl)[C@@H](Br)(Cl)",
            False,
            ["F[C@](Br)(Cl)[C@@H](Br)(Cl)", "F[C@@](Br)(Cl)[C@@H](Br)(Cl)"],
        ),
        ("F[C@H](Cl)Br", True, ["F[C@H](Cl)Br", "F[C@@H](Cl)Br"]),
        ("F[C@H](Cl)Br", False, ["F[C@H](Cl)Br"]),
    ],
)
@pytest.mark.parametrize(
    "toolkit_wrapper", [OpenEyeToolkitWrapper(), RDKitToolkitWrapper()]
)
def test_enumerate_stereoisomers(smiles, force_flip, expected, toolkit_wrapper):

    if (
        "@" in smiles
        and force_flip
        and isinstance(toolkit_wrapper, RDKitToolkitWrapper)
    ):

        pytest.skip(
            "the rdkit wrapper cannot flip existing stereo. see #892 of the OFF TK."
        )

    with global_toolkit_wrapper(toolkit_wrapper):

        molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

        stereoisomers = _enumerate_stereoisomers(
            molecule.to_openeye(),
            force_flip=force_flip,
            enum_nitrogen=True,
            verbose=False,
        )

    assert len(stereoisomers) == len(expected)

    expected = {
        Molecule.from_smiles(stereoisomer, allow_undefined_stereo=True).to_smiles(
            explicit_hydrogens=False, isomeric=True, mapped=False
        )
        for stereoisomer in expected
    }
    actual = {
        Molecule.from_openeye(stereoisomer, allow_undefined_stereo=True).to_smiles(
            explicit_hydrogens=False, isomeric=True, mapped=False
        )
        for stereoisomer in stereoisomers
    }

    assert expected == actual
