import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture
def thioacid(carbon, oxygen, sulfur, hydrogen, atom):
    bonders = ()
    deleters = ()
    return _TestCase(
        functional_group=stk.Thioacid(
            carbon=carbon,
            oxygen=oxygen,
            sulfur=sulfur,
            hydrogen=hydrogen,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, oxygen, sulfur, hydrogen, atom),
        bonders=bonders,
        deleters=deleters,
    )
