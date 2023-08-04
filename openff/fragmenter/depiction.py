import base64
import os
from typing import Collection, Dict, Optional

from jinja2 import Template
from openff.fragmenter.fragment import BondTuple, FragmentationResult
from openff.fragmenter.utils import get_map_index
from openff.toolkit.topology import Molecule
from openff.utilities import MissingOptionalDependencyError, requires_oe_module
from pkg_resources import resource_filename


@requires_oe_module("oechem")
def _oe_fragment_predicates(map_indices: Collection[int]):
    """Returns an atom and bond predicate which matches atoms whose map index
    appears in the specified ``map_indices`` collection"""

    from openeye import oechem

    class PredicateAtoms(oechem.OEUnaryAtomPred):
        def __init__(self, inner_map_indices):
            oechem.OEUnaryAtomPred.__init__(self)
            self.map_indices = inner_map_indices

        def __call__(self, atom):
            return atom.GetMapIdx() in self.map_indices

        def CreateCopy(self):
            return PredicateAtoms(self.map_indices).__disown__()

    class PredicateBonds(oechem.OEUnaryBondPred):
        def __init__(self, inner_map_indices):
            oechem.OEUnaryBondPred.__init__(self)
            self.map_indices = inner_map_indices

        def __call__(self, bond: oechem.OEBondBase):
            return (
                bond.GetBgn().GetMapIdx() in self.map_indices
                and bond.GetEnd().GetMapIdx() in self.map_indices
            )

        def CreateCopy(self):
            return PredicateBonds(self.map_indices).__disown__()

    return PredicateAtoms(map_indices), PredicateBonds(map_indices)


@requires_oe_module("oedepict")
def _oe_wbo_label_display(bond_tuples: Collection[BondTuple]):
    """Returns a ``OEDisplayBondPropBase`` subclass which will label bonds with the
    specified map indices with their WBO values if present.
    """

    from openeye import oedepict

    class LabelWibergBondOrder(oedepict.OEDisplayBondPropBase):
        def __init__(self, inner_bond_tuples):
            oedepict.OEDisplayBondPropBase.__init__(self)
            self.bond_tuples = inner_bond_tuples

        def __call__(self, bond):
            map_tuple = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())

            if (
                map_tuple not in self.bond_tuples
                and tuple(reversed(map_tuple)) not in self.bond_tuples
            ):
                return " "

            if "fractional_bond_order" not in bond.GetData():
                return " "

            return f"{bond.GetData('fractional_bond_order'):.2f}"

        def CreateCopy(self):
            return LabelWibergBondOrder(self.bond_tuples).__disown__()

    return LabelWibergBondOrder(bond_tuples)


@requires_oe_module("oedepict")
def _oe_render_parent(
    parent: Molecule,
    rotor_bonds: Optional[Collection[BondTuple]] = None,
    image_width: int = 572,
    image_height: int = 198,
) -> str:
    from openeye import oedepict

    rotor_bonds = [] if rotor_bonds is None else rotor_bonds

    # Map the OpenFF molecules into OE ones, making sure to explicitly set the atom
    # map on the OE object as this is not handled by the OpenFF toolkit.
    oe_parent = parent.to_openeye()

    for atom in oe_parent.GetAtoms():
        atom.SetMapIdx(get_map_index(parent, atom.GetIdx(), False))

    oedepict.OEPrepareDepiction(oe_parent)

    # Set-up common display options.
    image = oedepict.OEImage(image_width, image_height)

    display_options = oedepict.OE2DMolDisplayOptions(
        image_width, image_height, oedepict.OEScale_AutoScale
    )
    display_options.SetTitleLocation(oedepict.OETitleLocation_Hidden)
    display_options.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
    display_options.SetAtomLabelFontScale(1.2)
    display_options.SetBondPropertyFunctor(_oe_wbo_label_display(rotor_bonds))

    display = oedepict.OE2DMolDisplay(oe_parent, display_options)

    oedepict.OERenderMolecule(image, display)

    svg_contents = oedepict.OEWriteImageToString("svg", image)
    return svg_contents.decode()


@requires_oe_module("oedepict")
@requires_oe_module("oechem")
def _oe_render_fragment(
    parent: Molecule,
    fragment: Molecule,
    bond_indices: BondTuple,
    image_width: int = 283,
    image_height: int = 169,
) -> str:
    from openeye import oechem, oedepict

    # Map the OpenFF molecules into OE ones, making sure to explicitly set the atom
    # map on the OE object as this is not handled by the OpenFF toolkit.
    oe_parent = parent.to_openeye()

    for atom in oe_parent.GetAtoms():
        atom.SetMapIdx(get_map_index(parent, atom.GetIdx(), False))

    oedepict.OEPrepareDepiction(oe_parent)

    oe_fragment = fragment.to_openeye()

    for atom in oe_fragment.GetAtoms():
        atom.SetMapIdx(get_map_index(fragment, atom.GetIdx(), False))

    oe_parent_bond = oe_parent.GetBond(
        oe_parent.GetAtom(oechem.OEHasMapIdx(bond_indices[0])),
        oe_parent.GetAtom(oechem.OEHasMapIdx(bond_indices[1])),
    )

    # Set-up common display options.
    image = oedepict.OEImage(image_width, image_height)

    display_options = oedepict.OE2DMolDisplayOptions(
        image_width, image_height, oedepict.OEScale_AutoScale
    )

    display_options.SetTitleLocation(oedepict.OETitleLocation_Hidden)
    display_options.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
    display_options.SetAtomLabelFontScale(1.2)

    # display_options.SetBondPropertyFunctor(_oe_wbo_label_display({bond_indices}))

    display = oedepict.OE2DMolDisplay(oe_parent, display_options)

    fragment_atom_predicate, fragment_bond_predicate = _oe_fragment_predicates(
        {atom.GetMapIdx() for atom in oe_fragment.GetAtoms()}
    )

    not_fragment_atoms = oechem.OENotAtom(fragment_atom_predicate)
    not_fragment_bonds = oechem.OENotBond(fragment_bond_predicate)

    oedepict.OEAddHighlighting(
        display,
        oedepict.OEHighlightByColor(oechem.OEGrey, 0.75),
        not_fragment_atoms,
        not_fragment_bonds,
    )

    rotatable_bond = oechem.OEAtomBondSet()

    rotatable_bond.AddBond(oe_parent_bond)
    rotatable_bond.AddAtom(oe_parent_bond.GetBgn())
    rotatable_bond.AddAtom(oe_parent_bond.GetEnd())

    oedepict.OEAddHighlighting(
        display,
        oechem.OEColor(oechem.OELimeGreen),
        oedepict.OEHighlightStyle_BallAndStick,
        rotatable_bond,
    )

    oedepict.OERenderMolecule(image, display)

    svg_contents = oedepict.OEWriteImageToString("svg", image)
    return svg_contents.decode()


def _rd_render_parent(
    parent: Molecule,
    image_width: int = 572,
    image_height: int = 198,
) -> str:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.rdDepictor import Compute2DCoords

    rd_parent: Chem.Mol = parent.to_rdkit()
    rd_parent = Chem.RemoveHs(rd_parent)
    Compute2DCoords(rd_parent)

    drawer = Draw.MolDraw2DSVG(image_width, image_height)
    drawer.DrawMolecule(rd_parent)
    drawer.FinishDrawing()

    svg_contents = drawer.GetDrawingText()

    return svg_contents


def _rd_render_fragment(
    parent: Molecule,
    fragment: Molecule,
    bond_indices: BondTuple,
    image_width: int = 283,
    image_height: int = 169,
) -> str:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.rdDepictor import Compute2DCoords

    rd_parent: Chem.Mol = parent.to_rdkit()

    for atom in rd_parent.GetAtoms():
        atom.SetAtomMapNum(get_map_index(parent, atom.GetIdx(), False))

    rd_parent = Chem.RemoveHs(rd_parent)
    Compute2DCoords(rd_parent)

    map_indices = {*fragment.properties["atom_map"].values()} - {0}

    fragment_atom_indices = [
        atom.GetIdx()
        for atom in rd_parent.GetAtoms()
        if atom.GetAtomMapNum() in map_indices
    ]
    fragment_bond_indices = [
        bond.GetIdx()
        for bond in rd_parent.GetBonds()
        if bond.GetBeginAtom().GetAtomMapNum() in map_indices
        and bond.GetEndAtom().GetAtomMapNum() in map_indices
    ]

    rotatable_bond_index = [
        bond.GetIdx()
        for bond in rd_parent.GetBonds()
        if bond.GetBeginAtom().GetAtomMapNum() in bond_indices
        and bond.GetEndAtom().GetAtomMapNum() in bond_indices
    ]

    for atom in rd_parent.GetAtoms():
        atom.SetAtomMapNum(0)

    drawer = Draw.MolDraw2DSVG(image_width, image_height)

    draw_options = drawer.drawOptions()
    draw_options.useBWAtomPalette()

    drawer.DrawMolecule(
        rd_parent,
        highlightAtoms=fragment_atom_indices,
        highlightAtomColors={
            index: (52.0 / 255.0, 143.0 / 255.0, 235.0 / 255.0)
            for index in fragment_atom_indices
        },
        highlightBonds=fragment_bond_indices + rotatable_bond_index,
        highlightBondColors={
            index: (239.0 / 255.0, 134.0 / 255.0, 131.0 / 255.0)
            if index in rotatable_bond_index
            else (52.0 / 255.0, 143.0 / 255.0, 235.0 / 255.0)
            for index in fragment_bond_indices + rotatable_bond_index
        },
    )
    drawer.FinishDrawing()

    svg_contents = drawer.GetDrawingText()

    return svg_contents


def _compress_svg(svg_contents: str) -> str:
    encoded_image = base64.b64encode(svg_contents.encode()).decode()
    return f"data:image/svg+xml;base64,{encoded_image}"


def depict_fragmentation_result(result: FragmentationResult, output_file: str):
    """Generates a HTML report of fragments for a parent molecule with the rotatable
    bond highlighted.

    Parameters
    ----------
    result
        The result of fragmenting a molecule with a fragmentation engine.
    output_file : str
        The name of the file to write out to.
    """

    depict_fragments(
        result.parent_molecule,
        {fragment.bond_indices: fragment.molecule for fragment in result.fragments},
        output_file,
    )


def depict_fragments(
    parent: Molecule, fragments: Dict[BondTuple, Molecule], output_file: str
):
    """Generates a HTML report of fragments for a parent molecule with the rotatable
    bond highlighted.

    Parameters
    ----------
    parent
        The parent molecule that was fragmented.
    fragments
        The fragments generated around each rotatable bond.
    output_file : str
        The name of the file to write out to.
    """

    try:
        header_svg = _oe_render_parent(parent, [*fragments])
        fragment_svg = [
            _oe_render_fragment(parent, fragment, bond_tuple)
            for bond_tuple, fragment in fragments.items()
        ]

    except (ModuleNotFoundError, MissingOptionalDependencyError):
        header_svg = _rd_render_parent(parent)
        fragment_svg = [
            _rd_render_fragment(parent, fragment, bond_tuple)
            for bond_tuple, fragment in fragments.items()
        ]

    template_path = resource_filename(
        "openff.fragmenter", os.path.join("data", "report-template.html")
    )

    with open(template_path) as file:
        template = Template(file.read())

    rendered = template.render(
        parent=_compress_svg(header_svg),
        fragments=[_compress_svg(svg) for svg in fragment_svg],
    )

    with open(output_file, "w") as file:
        file.write(rendered)
