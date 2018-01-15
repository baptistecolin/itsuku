#!/usr/bin/env python3

import os
import struct
from hashlib import sha512
from math import ceil, log

n = 4 # number of dependencies
T = 2**5 # length of the main array
x = 64 # size of elements in the main array
M = 4 # size of elements in the Merkel Tree
S = 64 # size of elements of Y
L = 9 # length of one search
d = 1 # PoW difficulty (or strength)
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

def phis(seed,i,n):
    assert n>=1 and n<=11
    
    # The ouput is a list of length n, where the k-th element is phi_{k}(i)
    res = []
    res.append(i-1) # phi_{0}(i)
    if n>=2:
        phi_i = phi(seed,i)
        res.append(phi_i) # phi_{1}
        if n>=3:
            res.append(phi_i // 2)
            if n>=4:
                res.append((i - 1) // 2)
                if n>=5:
                    res.append((phi_i + i) // 2)
                    if n>=6:
                        res.append(3 * phi_i // 4)
                        if n>=7:
                            res.append(3 * i // 4)
                            if n>=8:
                                res.append(phi_i // 4)
                                if n>=9:
                                    res.append(i // 4)
                                    if n>=10:
                                        res.append(phi_i * 7 // 8)
                                        if n>=11:
                                            res.append(i * 7 // 8)
    return res

    # TODO : add better phi

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

def leading_zeros(d, x):
    # the input is a byte string x
    # it is converted to an int (big endian convertion)
    # which is converted to a string, corresponding to the binary representation of the int
    # the initial '0b' is stripped from the beginning of the string
    # then we add as much zeros as needed at the beginning of the string for it to be of length 64
    # if the string was already of length 64 then no zero is added
    # 
    # in the end, we get a string which is the binary representation of the int that x represents,
    # in big endian mode, of length 64, so that it makes sense to extract the first d digits with d up to 64
    #
    # The 64 limit is due to the fact that the output of a SHA512 hash is of length 64
    assert d <= 64
    return bin(int.from_bytes(x,'big')).lstrip('-0b').zfill(64)[:d]=='0'*d

def build_L(i, X, l, n=n):
    round_L = {} # will associate each index with the corresponding leaf and antecedent leaves

    for j in range(len(i)):
        
        if i[j] % l < n:
            # i[j] is such that X[i[j]] was built at step 1.a
            p = i[j] // l
            round_L[i[j]] = X[p*l:p*l+n]
        else :
            # i[j] is such that X[i[j]] was built at step 1.b
            seed = X[i[j]-1][:4]
            p = i[j] // l
            round_L[i[j]] = [ X[p*l + phi_k_i] for phi_k_i in phis(seed, i[j], n) ]
        
    return round_L

def PoW(I, T, n, P, M, L, S, d):
    X = memory_build(I, T, n, P, M)
    MT = merkle_tree(I, X, M)
    
    PSI = MT[0]
    
    # Choosing a nonce
    N = os.urandom(32)
   

    Y, OMEGA, i = compute_Y(I, X, L, S, N, PSI)
    counter = 0
    while not(leading_zeros(d, OMEGA)):
        Y, OMEGA = compute_Y(I,X,L,S,N,PSI)

        counter += 1
        if counter % 25 == 0:
            print("attempt n°"+str(counter))

    print("success on attempt #" + str(counter))
    
    round_L = build_L(i, X, l)
    
    # TODO : rest of the protocol

    return N, Y
