from opening import *
from math import log2
import pytest

def test_opening():
    # testing on examples that feature different amount of leaves
    for T in [4,8,32,1024]:
        H = int(log2(T))+1

        left_leaves = [i for i in range(T//2,T)]
        right_leaves = [i for i in range(0,T//2)]
        half_of_leaves = [2*i for i in range(0,T//2)]
        all_leaves = [i for i in range(0,T)]
        one_leaf = [0]

        left_leaves_expected_opening = [(1,0)] # only one leaf, the right son of the root
        right_leaves_expected_opening = [(1,1)] # only one leaf, the left son of the root
        half_of_leaves_expected_opening = [(H-1, 2*i+1) for i in range(0,T//2)] # remaining leaves
        all_leaves_expected_opening = [] # nothing ! everything can already be computed
        one_leaf_expected_opening = [(i,1) for i in range(H-1,0,-1)] # comb shaped opening

        assert opening(T, left_leaves) == left_leaves_expected_opening
        assert opening(T, right_leaves) == right_leaves_expected_opening
        assert opening(T, half_of_leaves) == half_of_leaves_expected_opening
        assert opening(T, all_leaves) == all_leaves_expected_opening
        assert opening(T, one_leaf) == one_leaf_expected_opening

def test_openingForOneArray():
    # testing for the same examples as above
    for T in [4,8,32,1024]:
        H = int(log2(T))+1

        left_leaves = [i for i in range(T//2,T)]
        right_leaves = [i for i in range(0,T//2)]
        half_of_leaves = [2*i for i in range(0,T//2)]
        all_leaves = [i for i in range(0,T)]
        one_leaf = [0]

        left_leaves_expected_opening = [1] # only one leaf, the right son of the root
        right_leaves_expected_opening = [2] # only one leaf, the left son of the root
        half_of_leaves_expected_opening = [2**(H-1)-1 + 2*i+1 for i in range(0,T//2)] # remaining leaves
        all_leaves_expected_opening = [] # nothing ! everything can already be computed
        one_leaf_expected_opening = [2**i for i in range(H-1,0,-1)] # comb shaped opening

        assert openingForOneArray(T, left_leaves) == left_leaves_expected_opening
        assert openingForOneArray(T, right_leaves) == right_leaves_expected_opening
        assert openingForOneArray(T, half_of_leaves) == half_of_leaves_expected_opening
        assert openingForOneArray(T, all_leaves) == all_leaves_expected_opening
        assert openingForOneArray(T, one_leaf) == one_leaf_expected_opening


@pytest.mark.skip(reason="to be filled")
def test_opening_2():
    # TODO : write test
    return None

