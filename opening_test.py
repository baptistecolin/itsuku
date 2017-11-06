from opening import *
from math import log2
import pytest

def test_opening():
    # testing on examples that feature 8 leaves
    T = 8
    H = int(log2(T))+1

    left_leaves = [i for i in range(T//2,T)]
    right_leaves = [i for i in range(0,T//2)]
    half_of_leaves = [2*i for i in range(0,T//2)]
    all_leaves = [i for i in range(0,T)]

    left_leaves_expected_opening = [(1,0)]
    right_leaves_expected_opening = [(1,1)]
    half_of_leaves_expected_opening = [(H-1, 2*i+1) for i in range(0,T//2)]
    all_leaves_expected_opening = []

    assert opening(T, left_leaves) == left_leaves_expected_opening
    assert opening(T, right_leaves) == right_leaves_expected_opening
    assert opening(T, half_of_leaves) == half_of_leaves_expected_opening
    assert opening(T, all_leaves) == all_leaves_expected_opening

@pytest.mark.skip(reason="to be filled")
def test_openingForOneArray():
    # TODO : write test
    return None

@pytest.mark.skip(reason="to be filled")
def test_opening_2():
    # TODO : write test
    return None

