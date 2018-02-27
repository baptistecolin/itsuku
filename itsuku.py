#!/usr/bin/env python3

import os
import struct
import json
from hashlib import sha512
from math import floor, ceil, log
from opening import openingForOneArray as opening
from collections import OrderedDict
# TODO : consider adding typing (import typing)

n = 4 # number of dependencies
T = 2**5 # length of the main array, non power of 2 should be allowed
x = 64 # size of elements in the main array
M = 64 # size of elements in the Merkel Tree
S = 64 # size of elements of Y
L = 9 # length of one search
d = b'\x00'*(S-1) + b'\xff' # PoW difficulty (or strength)
#l = 2**15 # length of segment
l = 2**5
P = T/l # number of independent sequences
I = os.urandom(M) # initial challenge (randomly generated M bytes array)

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

# TODO: implement function F

# return int n as a 4 byte string, for hashing purposes
def int_to_4bytes(n):
    assert type(n) == int
    return struct.pack('>I', n)

# build and return array X
# ??? this probably does not work if T is not a 2**.
def memory_build(I, T, n, P, x):
    X = [None] * T
    # ??? we should provide l, and expect l * P >= T
    l = T // P
    # ??? particular case
    assert float(l) == T / P

    # parallel segments
    for p in range(P):

        # Step 1.a: build initial elements out of i, p and I
        for i in range(n):
            X[p*l+i] = H(x, int_to_4bytes(i) + int_to_4bytes(p) + I)

        # Step 1.b: build elements that depend on antecedents using phi functions
        for i in range(n, l):
            seed = X[p*l + i-1][:4]
            hinput = b''
            # ??? this is a simplified version
            for phi_k_i in phis(seed, i, n):
                hinput += X[p*l + phi_k_i]

            X[p*l+i] = H(x, hinput)

    return X

# build merkle tree
def merkle_tree(I, X, M):
    T = len(X)

    # Step 2.a. : build Merkle-tree as an array
    B = [None] * (2*T-1)

    # Step 2.b. : Compute leaf elements out of hashes of X
    for i in range(T):
        B[i + T - 1] =  H(M, X[i] + I)

    # Step 2.c. : Compute intermediate elements as hashes of their sons
    for i in range(T-2, -1, -1): # Downward iteration from T-2 to 0, both included
        B[i] = H(M, B[2*i+1] + B[2*i+2] + I )

    return B

# Let's use a recursive function !
# T = number of leaves of the MT
# known_nodes = dictionnary { i: b } of known nodes of the MT
#               where 0 <= i < 2T-1 is the index of the node in the array representation of the tree
#               and b is the hash stored in the corresponding node
# index = index in the array representation of the MT of the hash we want to compute
def compute_merkle_tree_node(index, known_nodes, I, T, M):
    assert index < 2*T-1 , "Out of bound index : %i" % index
    if index in known_nodes:
        return known_nodes[index]
    else:
        return H(M,
                compute_merkle_tree_node(2*index+1, known_nodes, I, T, M) +
                compute_merkle_tree_node(2*index+2, known_nodes, I, T, M) + I )


# Surprisingly, there is no XOR operation for bytearrays, so this has to been done this way.
# See : https://bugs.python.org/issue19251
def xor(a,b):
    return bytes(x ^ y for x, y in zip(a,b))


def compute_Y(I, X, T, L, S, N, PSI, byte_order='big'):
    # Build array Y of length L+1
    Y = [None]*(L+1)

    # Initialization
    Y[0] = H(S, N + PSI + I)

    # Building the array
    i = [None]*L
    for j in range(1, L+1):
        # Step 5.a
        i[j-1] = int.from_bytes(Y[j-1], byte_order) % T
        # Step 5.b
        Y[j] = H(S, Y[j-1] + xor(X[i[j-1]], I))

    # computing OMEGA
    if len(Y)%2==1:
        OMEGA_input = b''.join(Y[:0:-1])
    else:
        OMEGA_input = b''.join(Y[::-1])
    OMEGA = H(S, OMEGA_input)

    return Y, OMEGA, i

def is_PoW_solved(d, x, S=S):
    assert len(x) == S
    assert len(d) == S
    # using the lexicographic order
    return x > d

def build_L(i, X, P, n):
    round_L = OrderedDict.fromkeys(i) # will associate each index with the corresponding leaf and antecedent leaves

    # computing l
    l = len(X)//P
    assert l == floor(len(X)//P)

    for i_j in round_L:
        p = i_j // l

        if i_j % l < n:
            # i[j] is such that X[i[j]] was built at step 1.a
            round_L[i_j] = X[p*l:p*l+n]
        else :
            # i[j] is such that X[i[j]] was built at step 1.b
            seed = X[i_j-1][:4]
            round_L[i_j] = [ X[p*l + phi_k_i] for phi_k_i in phis(seed, i_j%l , n) ]

    return round_L


# Given a round_L object, such as the one that is returned at the end of the Itsuku PoW,
# this function computes the list of the indexes of all the element
# of the corresponding array X that are known and stored in round_L
#
# Computing this information turns out to be necessary before computing the opening of a merkle tree.

def provided_indexes(round_L, P, T, n):
    l = T//P
    assert l == floor(T/P)

    res = list(round_L.keys())

    for i_j in round_L :

        p = i_j // l

        if i_j % l < n :
            # Case when round_L[i_j] items have been built at step (1.a)
            res += range(p*l, p*l+n)
        else :
            # Case when round_L[i_j] items have been built at step (1.b)
            seed = round_L[i_j][0][:4]
            res += [p*l + phi_k_i for phi_k_i in phis(seed, i_j % l, n) ] + [i_j]
            # One may say that adding i_j to the list is wrong because round_L does
            # not encapsulate X[i_j]. Whereas it is true that X[i_j] is not witholded
            # by round_L, X[i_j] can be recomputed from the elements of round_L,
            # thus making it an element which can be considered as known if round_L is known

    # This is intended to remove the likely duplicates
    # It turns out it is the fastest way to achieve deduplication
    res = list(dict.fromkeys(res))

    return res

def build_Z(round_L, MT, P, T, n):
    indexes = provided_indexes(round_L, P, T, n)
    opening_indexes = opening(T, indexes)
    Z = dict.fromkeys(opening_indexes)
    for k in Z :
        Z[k] = MT[k]

    return Z

def clean_Z(Z):
    return  { k: v.hex() for k,v in Z.items() }

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

def build_JSON_output(N, round_L, Z, P, T, n, I, M, L, S, x, d):
    data = {'answer':{}, 'params':{}}
    data['answer']['N'] = N.hex()
    data['answer']['round_L'] = trim_round_L(round_L, P, T, n)
    data['answer']['Z'] = clean_Z(Z)
    # no need to add i to the data because it can be obtain by extracting the keys of round_L
    data['params']['n'] = n
    data['params']['P'] = P
    data['params']['T'] = T
    data['params']['I'] = I.hex()
    data['params']['M'] = M
    data['params']['L'] = L
    data['params']['S'] = S
    data['params']['x'] = x
    data['params']['d'] = d.hex()

    return json.dumps(data, separators=(',',':'))


def PoW(I=I, T=T, n=n, P=P, M=M, L=L, S=S, x=x, d=d, debug=False):
    X = memory_build(I, T, n, P, x, M)
    MT = merkle_tree(I, X, M)
    PSI = MT[0]

    # Choosing a nonce. Could be a counter.
    N = os.urandom(32)

    Y, OMEGA, i = compute_Y(I, X, T, L, S, N, PSI)
    counter = 0
    while not(is_PoW_solved(d, OMEGA, S)):
        N = os.urandom(32) # Choosing a new nonce
        Y, OMEGA, i = compute_Y(I,X,T,L,S,N,PSI)

        counter += 1
        if counter % 25 == 0:
            print("attempt n°"+str(counter))

    print("success on attempt #" + str(counter))

    round_L = build_L(i, X, P, n)
    Z = build_Z(round_L, MT, P, T, n)

    json_output = build_JSON_output(
            N=N,
            round_L=round_L,
            Z=Z,
            P=P,
            T=T,
            n=n,
            I=I,
            M=M,
            L=L,
            S=S,
            x=x,
            d=d
        )
    if debug:
        return json_output, X, MT, PSI, N, I, Y, OMEGA, i, round_L, Z
    else:
        return json_output
