from fragmenter.utils import get_fgroup_smarts, get_fgroup_smarts_comb


def test_get_fgroup_smarts():
    """Tests that the `get_fgroup_smarts` utility returns correctly."""

    smarts = get_fgroup_smarts()

    assert "hydrazine" in smarts
    assert smarts["hydrazine"] == "[NX3][NX3]"

    assert "phosphon" not in smarts

    assert len(smarts) == 21


def test_get_fgroup_smarts_comb():
    """Tests that the `get_fgroup_smarts_comb` utility returns correctly."""

    smarts = get_fgroup_smarts_comb()

    assert "phosphon" in smarts
    assert smarts["phosphon"] == "[PX4](=[OX1])([OX2])([OX2])"

    assert "amide_2" not in smarts

    assert len(smarts) == 12
