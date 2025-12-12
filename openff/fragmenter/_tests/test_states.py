import pytest
from openff.toolkit.topology import Molecule

from openff.fragmenter.states import _enumerate_stereoisomers


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
        ("FC(Cl)Br", False, ["F[C@H](Cl)Br", "F[C@@H](Cl)Br"]),
        ("F[C@H](Cl)Br", True, ["F[C@H](Cl)Br", "F[C@@H](Cl)Br"]),
        ("F[C@H](Cl)Br", False, ["F[C@H](Cl)Br"]),
        ("ClC=CBr", False, [r"Cl/C=C/Br", r"Cl\C=C/Br"]),
        ("Cl/C=C/Br", False, [r"Cl/C=C/Br"]),
        ("Cl/C=C/Br", True, [r"Cl/C=C/Br", r"Cl\C=C/Br"]),
    ],
)
def test_enumerate_stereoisomers(smiles, force_flip, expected):
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

    stereoisomers = _enumerate_stereoisomers(molecule, force_flip=force_flip)

    assert len(stereoisomers) == len(expected)

    expected = {
        Molecule.from_smiles(stereoisomer, allow_undefined_stereo=True).to_smiles(
            explicit_hydrogens=False, isomeric=True, mapped=False
        )
        for stereoisomer in expected
    }
    actual = {
        stereoisomer.to_smiles(explicit_hydrogens=False, isomeric=True, mapped=False) for stereoisomer in stereoisomers
    }

    assert expected == actual
