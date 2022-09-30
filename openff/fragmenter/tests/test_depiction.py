import os

import pytest
from openff.fragmenter.depiction import (
    _oe_render_fragment,
    _oe_render_parent,
    _rd_render_fragment,
    _rd_render_parent,
    depict_fragmentation_result,
    depict_fragments,
)
from openff.fragmenter.fragment import Fragment, FragmentationResult
from openff.toolkit.topology import Molecule


@pytest.mark.parametrize("draw_function", [_oe_render_parent, _rd_render_parent])
def test_xx_render_parent(draw_function):

    try:
        svg_contents = draw_function(
            Molecule.from_smiles(
                "[C:1]([H:5])([H:6])([H:7])[C:2]([H:8])([H:9])[C:3]([H:10])([H:11])[C:4]([H:12])([H:13])([H:14])"
            )
        )
    except ModuleNotFoundError as e:
        pytest.skip(str(e))

    assert isinstance(svg_contents, str)
    assert "svg" in svg_contents


@pytest.mark.parametrize("draw_function", [_oe_render_fragment, _rd_render_fragment])
def test_xx_render_fragment(draw_function):

    try:

        svg_contents = draw_function(
            Molecule.from_smiles(
                "[C:1]([H:5])([H:6])([H:7])[C:2]([H:8])([H:9])[C:3]([H:10])([H:11])[C:4]([H:12])([H:13])([H:14])"
            ),
            Molecule.from_smiles(
                "[C:1]([H:3])([H:4])([H:5])[C:2]([H:6])([H:7])([H:8])"
            ),
            (1, 2),
        )

    except ModuleNotFoundError as e:
        pytest.skip(str(e))

    assert isinstance(svg_contents, str)
    assert "svg" in svg_contents


def test_depict_fragments(tmpdir):

    output_file = os.path.join(tmpdir, "report.html")

    depict_fragments(
        Molecule.from_smiles(
            "[C:1]([H:5])([H:6])([H:7])[C:2]([H:8])([H:9])[C:3]([H:10])([H:11])[C:4]([H:12])([H:13])([H:14])"
        ),
        {
            (1, 2): Molecule.from_smiles(
                "[C:1]([H:3])([H:4])([H:5])[C:2]([H:6])([H:7])([H:8])"
            )
        },
        output_file,
    )

    with open(output_file) as file:
        contents = file.read()

    assert "<html>" in contents
    assert "<img src=" in contents


def test_depict_fragmentation_result(tmpdir):

    output_file = os.path.join(tmpdir, "report.html")

    depict_fragmentation_result(
        FragmentationResult(
            parent_smiles="[C:1]([H:5])([H:6])([H:7])[C:2]([H:8])([H:9])[C:3]([H:10])([H:11])[C:4]([H:12])([H:13])([H:14])",
            fragments=[
                Fragment(
                    smiles="[C:1]([H:3])([H:4])([H:5])[C:2]([H:6])([H:7])([H:8])",
                    bond_indices=(1, 2),
                )
            ],
            provenance={},
        ),
        output_file,
    )

    with open(output_file) as file:
        contents = file.read()

    assert "<html>" in contents
    assert "<img src=" in contents
