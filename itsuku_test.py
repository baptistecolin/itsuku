from itsuku import *
import pytest

def test_phi():
    # shoud fail if the seed is not 4 bytes long
    seed = int_to_4bytes(123)
    seed_3bytes = seed[:3]
    seed_8bytes = seed+seed

    assert len(seed_3bytes) == 3
    assert len(seed_8bytes) == 8

    with pytest.raises(AssertionError):
        phi(seed_3bytes, 4)
    with pytest.raises(AssertionError):
        phi(seed_5bytes, 4)


def test_phis():
    # TODO : write test
    return None

def test_H():
    # TODO : write test
    return None

def test_int_to_4bytes():
    # TODO : write test
    return None

def test_memory_build():
    # TODO : write test
    return None

def test_merkle_tree():
    # TODO : write test
    return None

def test_compute_Y():
    # TODO : write test
    return None

def test_opening():
    # TODO : write test
    return None

def test_PoW():
    # TODO : write test
    return None
