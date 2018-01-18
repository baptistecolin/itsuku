#!/usr/bin/env python3

import os
import struct
from hashlib import sha512
from math import floor, ceil, log
from opening import openingForOneArray as opening

n = 4 # number of dependencies
T = 2**5 # length of the main array
x = 64 # size of elements in the main array
M = 4 # size of elements in the Merkel Tree
S = 64 # size of elements of Y
L = 9 # length of one search
d = b'\x00'*63 + b'\xff' # PoW difficulty (or strength)
#l = 2**15 # length of segment
l = 2**5
P = T/l # number of independent sequences
I = os.urandom(M) # initial challenge (randomly generated M bytes array)

HASH = "sha512" # hash function


def phi(seed, i, byte_order='big', method='high-level'):
    # Will only work as expected if the seed is 4 bytes long
    assert len(seed) == 4

    J = int.from_bytes(seed, byte_order)
    R = i-1
    
    if method=='high-level':
        res = R*(1-((J*J)//(2**64)))
    else:
        # We are using the operations suggested at page 7 in https://www.cryptolux.org/images/0/0d/Argon2.pdf
        x = (J**2)//(2**32)
        y = ((i-1)*x)//(2**32)
        res = i - 1 - y
    
    return res

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

def phis(seed,i,n):
    assert n>=1 and n<=11
    phi_i = phi(seed, i)
    return [ phi(i, phi_i) for phi in PHI_K[:n] ]

def H(M,x,method=HASH):
    # Encapsulate hashing operations such as digest, update ... for better readability
    
    if method == "sha512":
        hashfunc = sha512() # it is important that a new hash function is instanciated every time
                            # otherwise, the output would depend on the previous inputs ...
        hashfunc.update(x)
        output = hashfunc.digest()
        return output[:M]

# Turns the int 1024 into the byte string b'\x00\x00\x04\x00', that is fit for hashing
def int_to_4bytes(n):
    return struct.pack('>I', n)

def memory_build(I, T, n, P, M):
    # Step (1)
    # Building a challenge dependent memory
    X = [None]*T
    
    assert T//P == floor(T/P)
    l = T//P
    # Step (1.a)
    for p in range(P):
        for i in range(n):
            hash_input = int_to_4bytes(i) + int_to_4bytes(p) + I
    
            X[p*l+i] = H(x, hash_input)
    
    # Step (1.b)
    for p in range(P):
        for i in range(n,l):
            # The seed that is used by phi is the 4 first bytes of X[p*l+i-1]
            seed = X[p*l+i-1][:4]

            # computing phi_{k}(i) for all k up until n
            phis_list = phis(seed,i,n)
    
            # building the input of the hash function
            hash_input = b""
            for phi in phis_list:
                hash_input += X[p*l+ phi]
            
            # inserting the computed hash in the array
            X[p*l+i] = H(x, hash_input)

    return X


def merkle_tree(I, X, M):
    # Building the Merkle Tree
    # It will be implemented as an array, each element being a node
    # The node at index i has its left son at index 2*i+1, and its right son at index 2*i+2
    # The array is of length 2T-1, with T being the length of X (full binary tree)
    
    # The leaves of the tree are the elements of X. Thus, MT[-T:] == hash(X). 
    MT = [None]*(2*len(X)-1)
    MT[-T:] = [ H(M,x) for x in X ] 

    # Building the non-leaf nodes
    for i in range(len(X)-2,-1,-1): # Decreasing iteration from len(X)-1 to 0, both included
        MT[i] = H(M, MT[2*i+1] + MT[2*i+2] + I ) #Hash of left son + right son + challenge
    
    return MT

    # TODO : add exceptions so it manages to build a merkle tree that has a number of leaves that is not 2^n

# Surprisingly, there is no XOR operation for bytearrays, so this has to been done this way.
# See : https://bugs.python.org/issue19251
def xor(a,b):
    return bytes(x ^ y for x, y in zip(a,b))


def compute_Y(I, X, L, S, N, PSI, byte_order='big'):
    # Build array Y of length L+1
    Y = [None]*(L+1)

    # Initialization
    Y[0] = H(S, N + PSI + I)
    
    # Building the array
    i = [None]*L
    for j in range(1, L+1):
        # Step 5.a
        i[j-1] = int.from_bytes(Y[j-1], byte_order) % len(X)
        # Step 5.b
        Y[j] = H(S, Y[j-1] + xor(X[i[j-1]], I))

    # computing OMEGA
    if len(Y)%2==1:
        OMEGA_input = b''.join(Y[:0:-1])
    else:
        OMEGA_input = b''.join(Y[::-1])
    OMEGA = H(S, OMEGA_input)
    
    return Y, OMEGA, i

def is_PoW_solved(d, x):
    assert len(x) == 64
    assert len(d) == 64
    # using the lexicographic order
    return x > d

def build_L(i, X, P, n):
    round_L = {} # will associate each index with the corresponding leaf and antecedent leaves
    
    # computing l
    l = len(X)//P
    assert l == floor(len(X)//P)

    for j in range(len(i)):
        p = i[j] // l
        
        if i[j] % l < n:
            # i[j] is such that X[i[j]] was built at step 1.a
            round_L[i[j]] = X[p*l:p*l+n]
        else :
            # i[j] is such that X[i[j]] was built at step 1.b
            seed = X[i[j]-1][:4]
            round_L[i[j]] = [ X[p*l + phi_k_i] for phi_k_i in phis(seed, i[j]%l , n) ]
        
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
            
            res += [p*l + phi_k_i for phi_k_i in phis(seed, i_j, n) ] + [i_j]
            # One may say that adding i_j to the list is wrong because round_L does
            # not encapsulate X[i_j]. Whereas it is true that X[i_j] is not witholded
            # by round_L, X[i_j] can be recomputed from the elements of round_L,
            # thus making it an element which can be considered as known if round_L is known
    
    # This is intended to remove the likely dupicates
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

def build_JSON_output():
    # TODO: write it !
    return None

def PoW(I, T, n, P, M, L, S, d):
    X = memory_build(I, T, n, P, M)
    MT = merkle_tree(I, X, M)
    
    PSI = MT[0]
    
    # Choosing a nonce
    N = os.urandom(32)
   

    Y, OMEGA, i = compute_Y(I, X, L, S, N, PSI)
    counter = 0
    while not(is_PoW_solved(d, OMEGA)):
        Y, OMEGA = compute_Y(I,X,L,S,N,PSI)

        counter += 1
        if counter % 25 == 0:
            print("attempt n°"+str(counter))

    print("success on attempt #" + str(counter))
    
    round_L = build_L(i, X, P, T, n)
    Z = build_Z(round_L, MT, P, T, n)
    
    # TODO : rest of the protocol

    return N, Y
