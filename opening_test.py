from opening import *
import pytest

def test_opening():
    # testing on examples that feature 8 leaves
    T = 8

    left_leaves = [4,5,6,7]
    right_leaves = [0,1,2,3]
    half_of_leaves = [2*i for i in range(0,4)]

    left_leaves_expected_opening = [(1,0)]
    right_leaves_expected_opening = [(1,1)]
    half_of_leaves_expected_opening = [(3, 2*i+1) for i in range(0,4)]

    assert opening(T, left_leaves) == left_leaves_expected_opening
    assert opening(T, right_leaves) == right_leaves_expected_opening
    assert opening(T, half_of_leaves) == half_of_leaves_expected_opening

@pytest.mark.skip(reason="to be filled")
def test_openingForOneArray():
    # TODO : write test
    return None

@pytest.mark.skip(reason="to be filled")
def test_opening_2():
    # TODO : write test
    return None

