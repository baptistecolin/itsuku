#!/usr/bin/env python3

import sys
import os
import struct
import json
from hashlib import sha512
from math import floor, ceil, log
from opening import openingForOneArray as opening
from collections import OrderedDict
# TODO : consider adding typing (import typing)

HASH = 'sha512' # hash function

# compute the Argon2 phi function
def phi(seed, i, byte_order='big'):
    # Will only work as expected if the seed is 4 bytes long
    assert type(seed) == bytes and len(seed) == 4

    j = int.from_bytes(seed, byte_order)

    # operations page 7 in https://www.cryptolux.org/images/0/0d/Argon2.pdf
    x = (j**2) // (2**32)
    y = ((i-1)*x) // (2**32)
    res = i - 1 - y

    return res

# list of phi functions, used up to n (number of dependencies)
PHI_K = [
    lambda i, pi: i-1,
    lambda i, pi: pi,
    lambda i, pi: pi//2,
    lambda i, pi: (i-1)//2,
    lambda i, pi: (pi + i)//2,
    lambda i, pi: (3*pi)//4,
    lambda i, pi: 3*i//4,
    lambda i, pi: pi//4,
    lambda i, pi: i//4,
    lambda i, pi: pi*7//8,
    lambda i, pi: i*7//8
]

# return all X[i] dependencies as a set of indexes
def phis(seed, i, n):
    assert type(seed) == bytes
    assert 1 <= n and n <= len(PHI_K)
    phi_i = phi(seed, i)
    return [ phi(i, phi_i) for phi in PHI_K[:n] ]

# return a M bytes hash of x
def H(M, x, method=HASH):
    # manual type check:-)
    assert type(M) == int and type(x) == bytes
    # Encapsulate hashing operations such as digest, update ... for better readability
    if method == 'sha512':
        h = sha512()
        h.update(x)
        return h.digest()[:M]
    else:
        raise Exception("unexpected hash '%s'" % method)

# ??? TODO: implement function F

# return int n as a 4 byte string, for hashing purposes
def int_to_4bytes(n):
    assert type(n) == int
    return struct.pack('>I', n)

# help, some redundancy
def _direct_X_i(x, I, p, k, l, n):
    assert k < n and n <= l
    return H(x, int_to_4bytes(k) + int_to_4bytes(p) + I)

# ??? FIXME this is not the expected formula
def _indirect_X_i(x, I, p, k, l, n, X):
    assert n <= k and k < l
    i = p*l + k
    assert ((type(X) is list and i-1 < len(X)) or \
            (type(X) is dict and i-1 in X))
    seed = X[i-1][:4]
    data = b''
    for phi in phis(seed, k, n):
        assert phi <= i-1 and \
            ((type(X) is list and p*l+phi < len(X)) or \
             (type(X) is dict and p*l+phi in X))
        data += X[p*l+phi]
    data += I
    return H(x, data)

# computation of X[i]
def compute_X_i(x, I, i, l, n):
    p, k = i // l, i % l
    if k < n:
        return _direct_X_i(x, I, p, k, l, n)
    else:
        return _indirect_X_i(x, I, p, k, l, n, X)

# build and return array X
# ??? this probably does not work if T is not a 2**.
def build_X(I, T, l, n, x):
    X = [None] * T
    P = (T + (l - 1)) // l
    # ??? particular case
    #assert float(l) == T / P
    # parallel segments
    for p in range(P):

        # Step 1.a: build initial elements out of i, p and I
        for k in range(n):
            X[p*l+k] = _direct_X_i(x, I, p, k, l, n)

        # Step 1.b: build elements that depend on antecedents using phi functions
        for k in range(n, l):
            X[p*l+k] = _indirect_X_i(x, I, p, k, l, n, X)

    return X


# rebuild a partial X
def rebuild_X(rL, I, l, n, x):
    X = {}
    for i, xs in rL.items():
        p, k = i // l, i % l
        if k < n:
            X_i = _direct_X_i(x, I, p, k, l, n)
        else:
            assert type(xs) is list
            seed = xs[0][:4]
            for j, v in zip(phis(seed, k, n), xs):
                if p*l+j in X:
                    assert X[p*l+j] == v
                else:
                    X[p*l+j] = v
            X_i = _indirect_X_i(x, I, p, k, l, n, X)
        if i in X:
            assert X[i] == X_i
        else:
            X[i] = X_i
    #print("nX=%s" % { k:v.hex() for k,v in X.items() })
    return X

def _cmp_MT_leaf(I, Xi, M):
    return H(M, Xi + I)

def _cmp_MT_node(I, X1, X2, M):
    return H(M, X1 + X2 + I)

# build merkle tree
# ??? TODO should work for non 2**
def build_MT(I, X, M):
    T = len(X)

    # Step 2.a. : build Merkle-tree as an array
    B = [None] * (2*T-1)

    # Step 2.b. : Compute leaf elements out of hashes of X
    for i in range(T):
        B[i + T - 1] = _cmp_MT_leaf(I, X[i], M)

    # Step 2.c. : Compute intermediate elements as hashes of their sons
    for i in range(T-2, -1, -1): # Downward iteration from T-2 to 0, both included
        B[i] = _cmp_MT_node(I, B[2*i+1], B[2*i+2], M)

    return B

# sigh... this should be in Python
from sortedcontainers import SortedSet

# rebuild partial Merkle Tree from available informations
# this is a bottom-up version of recursive "compute_MT_node",
# which checks that all values are used as expected
def rebuild_MT(rZ, I, X, M, T):
    B = {}
    for i, v in X.items():
        B[i + T - 1] = _cmp_MT_leaf(I, v, M)
    for i, v in rZ.items():
        assert i not in B
        B[i] = v
    indexes = SortedSet(B.keys())
    while len(indexes) >= 2:
        i2, i1 = indexes.pop(), indexes.pop()
        assert i1 + 1 == i2 and i2 % 2 == 0
        i0 = i1 // 2
        assert not i0 in indexes
        indexes.add(i0)
        B[i0] = _cmp_MT_node(I, B[i1], B[i2], M)
    assert len(indexes) == 1 and indexes.pop() == 0
    return B

# Let's use a recursive function !
# T = number of leaves of the MT
# known_nodes = dictionnary { i: b } of known nodes of the MT
#               where 0 <= i < 2T-1 is the index of the node in the array representation of the tree
#               and b is the hash stored in the corresponding node
# index = index in the array representation of the MT of the hash we want to compute
def compute_MT_node(index, known_nodes, I, T, M):
    assert index < 2*T-1 , "Out of bound index : %i" % index
    if index in known_nodes:
        return known_nodes[index]
    else:
        return H(M,
                compute_MT_node(2*index+1, known_nodes, I, T, M) +
                compute_MT_node(2*index+2, known_nodes, I, T, M) + I )

# Surprisingly, there is no XOR operation for bytearrays, so this has to been done this way.
# See : https://bugs.python.org/issue19251
def xor(a,b):
    # ensure that len(a) >= len(b)
    if len(a) < len(b):
        a, b = b, a
    return bytes(x ^ y for x, y in zip(a, b'\x00' * (len(a) - len(b)) + b))

# compute the Y sequence from nonce and other stuff
# X maybe a full or partial array
def compute_Y(I, X, T, L, S, N, Psi, byte_order='big'):
    # build array Y of length L+1
    Y = [None] * (L+1)

    # initialization
    Y[0] = H(S, N + Psi + I)

    # build array Y and keep used X indexes
    i = [None] * L
    for j in range(1, L+1):
        # Step 5.a
        # ??? a full modulo is expensive for non 2**
        # should it rather be on a few bytes?
        i[j-1] = int.from_bytes(Y[j-1], byte_order) % T
        # Step 5.b
        Y[j] = H(S, Y[j-1] + xor(X[i[j-1]], I))

    # Compute final Omega
    Omega = H(S, xor(b''.join(Y[:0:-1] if len(Y) % 2 == 1 else Y[::-1]), I))

    return Y, Omega, i

# check PoW solution wrt expected difficulty
def is_PoW_solved(d, x, S):
    assert len(x) == S and len(d) == S
    # lexicographic order
    return x > d

# debug helper
def print_rL(name, rL, l, n):
    print("%s = {" % name)
    for i, v in rL.items():
        if type(v) is list:
            print("  %d:%s," % (i, [x.hex() for x in v]))
        else:
            print("  %d:'%s'," % (i, v.hex()))
    print("}")

# return the roundL structure
# which maps selected indexes with their antecedents value in X so that they can be recomputed
def build_rL(rI, X, l, n):
    rL = {}
    for ij in rI:
        p, k = ij // l, ij % l
        if k < n:
            # i[j] is such that X[i[j]] was built at step 1.a
            # ??? useless, can be recomputed
            #round_L[ij] = X[ij]
            rL[ij] = []
        else :
            # i[j] is such that X[i[j]] was built at step 1.b
            seed = X[ij-1][:4]
            # ??? we could skip those which can be recomputed?
            rL[ij] = [ X[p*l + phi] for phi in phis(seed, k, n) ]

    return rL


# Given a round_L object, such as the one that is returned at the end of the Itsuku PoW,
# this function computes the list of the indexes of all the element
# of the corresponding array X that are known and stored in round_L
#
# Computing this information turns out to be necessary before computing the opening of a merkle tree.

# returns the set of indexes provided directly or indirectly with in roundL
def get_provided_indexes(rL, T, l, n):
    res = set()
    for i in rL:
        p, k = i // l, i % l
        res.add(i)
        if k >= n:
            # Case when round_L[i_j] items have been built at step (1.b)
            seed = rL[i][0][:4]
            # argh, union is a function...
            res = res.union((p*l + phi) for phi in phis(seed, k, n))
    return res

def build_rZ(rL, MT, T, l, n):
    rZ = {}
    for k in opening(T, get_provided_indexes(rL, T, l, n)):
        rZ[k] = MT[k]
    return rZ

# ???
def clean_Z(Z):
    return { k: v.hex() for k,v in Z.items() }

# ???
def trim_round_L(round_L, P, T, n):
    l = T//P
    assert l == T/P

    # Only keeping elements built a step 1.b, since elements
    # built at step 1.a can be recomputed knowing only I
    # Also, converting the bytearray elements of X to a format
    # that can be JSON-serialized
    trimmed_round_L = OrderedDict.fromkeys(round_L.keys())

    for k in trimmed_round_L:
        if k % l >= n:
            trimmed_round_L[k]= [ item.hex() for item in round_L[k] ]
        else:
            trimmed_round_L[k] = []

    return trimmed_round_L

# full version
def build_JSON_output(N, round_L, Z, P, T, n, I, M, L, S, x, d):
    data = {'answer':{}, 'params':{}}
    # needed PoW
    data['answer']['N'] = N.hex()
    data['answer']['rL'] = trim_round_L(round_L, P, T, n)
    data['answer']['Z'] = clean_Z(Z)
    # no need to add i to the data because it can be obtain by extracting the keys of round_L
    # informational
    data['params']['I'] = I.hex()
    data['params']['T'] = T # P * l <= T
    data['params']['P'] = P
    data['params']['l'] = l
    data['params']['n'] = n
    data['params']['M'] = M
    data['params']['L'] = L
    data['params']['S'] = S
    data['params']['x'] = x
    data['params']['d'] = d.hex()

    return json.dumps(data, separators=(',',':'))

# minimal json export
def exportPoW(N, rL, rZ):
    data = {
        'N': N.hex(),
        'L': { i: [ j.hex() for j in v ] for i, v in rL.items() },
        'Z': { i: v.hex() for i, v in rZ.items() }
    }
    return json.dumps(data)

# reverse of exportPoW
def importPoW(s):
    data = json.loads(s)
    N = bytes.fromhex(data['N'])
    rL = { int(i): [ bytes.fromhex(j) for j in v ] for i, v in data['L'].items() }
    rZ = { int(i): bytes.fromhex(v) for i, v in data['Z'].items() }
    return N, rL, rZ

# nonce size?
def solvePoW(I, T, l, n, M, L, S, x, d):
    X = build_X(I, T, l, n, x)
    B = build_MT(I, X, M)
    Psi = B[0]
    counter = 0
    while True:
        counter += 1
        # Choose nonce, could be a counter.
        N = os.urandom(8)
        Y, Omega, rI = compute_Y(I, X, T, L, S, N, Psi)
        # sigh, Python is still missing a do/while loop
        if Omega < d:
            break
    rL = build_rL(rI, X, l, n)
    rZ = build_rZ(rL, B, T, l, n)
    return exportPoW(N, rL, rZ), Omega, counter

def checkPoW(I, T, l, n, M, L, S, x, d, json_in):
    nN, nrL, nrZ = importPoW(json_in)
    nX = rebuild_X(nrL, I, l, n, x)
    nB = rebuild_MT(nrZ, I, nX, M, T)
    nPsi = nB[0]
    nY, nOmega, nrI = compute_Y(I, nX, T, L, S, nN, nPsi)
    return nOmega < d, nOmega

# UNUSED
# test values?
n = 4 # number of dependencies
T = 2**5 # length of the main array, non power of 2 should be allowed
x = 64 # size of elements in the main array
M = 64 # size of elements in the Merkel Tree
S = 64 # size of elements of Y
L = 9 # length of one search
d = b'\x00' + b'\xff' * (S-1) # PoW difficulty (or strength)
#l = 2**15 # length of segment
l = 2**5
P = T/l # number of independent sequences
I = os.urandom(64) # initial challenge (randomly generated M bytes array)

# needed seed size is 4 bytes
assert 4 <= x
assert M <= x
