from itsuku import *
import pytest
from hashlib import sha512

def test_phi():
    # it shoud fail if the seed is not 4 bytes long
    seed = int_to_4bytes(123)
    seed_3bytes = seed[:3]
    seed_8bytes = seed+seed

    assert len(seed_3bytes) == 3
    assert len(seed_8bytes) == 8

    with pytest.raises(AssertionError):
        phi(seed_3bytes, 4)
    with pytest.raises(AssertionError):
        phi(seed_8bytes, 4)

    # it should return the same result using the high-level or low level algorithm
    assert phi(seed, 4, method='high-level') == phi(seed, 4, method='low-level')
    assert phi(seed, 28, method='high-level') == phi(seed, 28, method='low-level')
    assert phi(seed, 1024, method='high-level') == phi(seed, 1024, method='low-level')

    # it should return a value phi(i) inferior to i
    assert phi(seed, 2) < 2
    assert phi(seed, 4) < 4
    assert phi(seed, 256) < 256
    assert phi(seed, 1024) < 1024


def test_phis():
    seed = int_to_4bytes(256)

    # it should output an array of length n
    for n in range(1,12):
        assert len(phis(seed, 10, n)) == n

    # TODO : more tests, probably

def test_H():
    # it should return a bytes array of length M
    x = int_to_4bytes(123456)
    for i in range(1,15):
        assert len(H(i,x)) == i
    
    # the H function should output the last M bytes of the sha512 hash of i
    for i in range(1,15):
        sha = sha512() # resetting the sha function
        sha_input = int_to_4bytes(123)
        sha.update(sha_input)
        assert sha.digest()[:10] == H(10,sha_input)
        


def test_int_to_4bytes():
    # it should always return a 4 bytes string
    assert len(int_to_4bytes(0)) == 4
    assert len(int_to_4bytes(1)) == 4
    assert len(int_to_4bytes(10)) == 4
    assert len(int_to_4bytes(1024)) == 4
    assert len(int_to_4bytes(4294967295)) == 4
    assert isinstance(int_to_4bytes(0), bytes)

    # it should return the expected hex bytes string, in a big endian order 
    assert int_to_4bytes(0) == b"\x00\x00\x00\x00"
    assert int_to_4bytes(1) == b"\x00\x00\x00\x01"
    assert int_to_4bytes(2) == b"\x00\x00\x00\x02"
    assert int_to_4bytes(16) == b"\x00\x00\x00\x10"
    assert int_to_4bytes(256) == b"\x00\x00\x01\x00"

@pytest.mark.skip(reason="to be filled")
def test_memory_build():
    # TODO : write test
    return None

@pytest.mark.skip(reason="to be filled")
def test_merkle_tree():
    # TODO : write test
    return None

@pytest.mark.skip(reason="to be filled")
def test_compute_Y():
    # TODO : write test
    return None

@pytest.mark.skip(reason="to be filled")
def test_opening():
    # TODO : write test
    return None

@pytest.mark.skip(reason="to be filled")
def test_PoW():
    # TODO : write test
    return None
