from opening import *
from itsuku import compute_merkle_tree_node
from math import log2
import pytest
import os
import random

def test_opening():
    # testing on examples that feature different amount of leaves
    for T in [4,8,32,1024]:
        H = int(log2(T))+1

        left_leaves = [i for i in range(T//2,T)]
        right_leaves = [i for i in range(0,T//2)]
        half_of_leaves = [2*i for i in range(0,T//2)]
        all_leaves = [i for i in range(0,T)]
        one_leaf = [0]
        empty = []

        left_leaves_expected_opening = [(1,0)] # only one leaf, the right son of the root
        right_leaves_expected_opening = [(1,1)] # only one leaf, the left son of the root
        half_of_leaves_expected_opening = [(H-1, 2*i+1) for i in range(0,T//2)] # remaining leaves
        all_leaves_expected_opening = [] # nothing ! everything can already be computed
        one_leaf_expected_opening = [(i,1) for i in range(H-1,0,-1)] # comb shaped opening
        empty_expected_opening = [(H-1,t) for t in range(T)] # all the leaves
        
        assert opening(T, left_leaves) == left_leaves_expected_opening
        assert opening(T, right_leaves) == right_leaves_expected_opening
        assert opening(T, half_of_leaves) == half_of_leaves_expected_opening
        assert opening(T, all_leaves) == all_leaves_expected_opening
        assert opening(T, one_leaf) == one_leaf_expected_opening
        assert opening(T, empty) == empty_expected_opening

def test_openingForOneArray():
    # testing for the same examples as above
    for T in [4,8,32,1024]:
        H = int(log2(T))+1

        left_leaves = [i for i in range(T//2,T)]
        right_leaves = [i for i in range(0,T//2)]
        half_of_leaves = [2*i for i in range(0,T//2)]
        all_leaves = [i for i in range(0,T)]
        one_leaf = [0]
        empty = []

        left_leaves_expected_opening = [1] # only one leaf, the right son of the root
        right_leaves_expected_opening = [2] # only one leaf, the left son of the root
        half_of_leaves_expected_opening = [2**(H-1)-1 + 2*i+1 for i in range(0,T//2)] # remaining leaves
        all_leaves_expected_opening = [] # nothing ! everything can already be computed
        one_leaf_expected_opening = [2**i for i in range(H-1,0,-1)] # comb shaped opening
        empty_expected_opening = [t+T-1 for t in range(T)] # all the leaves

        assert openingForOneArray(T, left_leaves) == left_leaves_expected_opening
        assert openingForOneArray(T, right_leaves) == right_leaves_expected_opening
        assert openingForOneArray(T, half_of_leaves) == half_of_leaves_expected_opening
        assert openingForOneArray(T, all_leaves) == all_leaves_expected_opening
        assert openingForOneArray(T, one_leaf) == one_leaf_expected_opening
        assert openingForOneArray(T, empty) == empty_expected_opening
        
        # The opening + the initial leaves should be enough to enable us to compute the merkle tree root
        # Therefore, the following instructions shouldn't fail
        for t in [ left_leaves, right_leaves, half_of_leaves, all_leaves, one_leaf ]:
            known_nodes = {k: b'\x00'*64 for k in [i + (T-1) for i in t] + openingForOneArray(T, t) }
            I = os.urandom(64)
            compute_merkle_tree_node(0, known_nodes, I, T, 64)

        # Now we're even going to do it with randomly generated lists of indexes
        def random_list_of_indexes(T):
            res = []
            for index in range(T):
                if bool(random.getrandbits(1)):
                    res.append(index)

            return res

        for t in [
                    random_list_of_indexes(T),
                    random_list_of_indexes(T),
                    random_list_of_indexes(T),
                    random_list_of_indexes(T),
                    random_list_of_indexes(T)
                 ]:
            known_nodes = {k: b'\x00' for k in [i + (T-1) for i in t] + openingForOneArray(T, t) }
            I = os.urandom(64)
            compute_merkle_tree_node(0, known_nodes, I, T, 64)

def test_opening_2():
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

