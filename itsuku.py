#!/usr/bin/env python

from phis import phis
from pyblake2 import blake2b
from hashlib import sha512
from math import ceil, log

n = 2 # number of dependencies
T = 2**10 # length of the main array
x = 64 # size of elements in the main array
M = 64 # size of elemets in the Merkel Tree
L = ceil(3.3*log(T,2)) # length of one search
d = 5 # difficulty of the PoW
P = 1 # number of independent sequences
l=T/P # length of one independent sequence

X = [None]*T # memory
H = sha512() # hash function

if M == 64:
    I = 'N\x8a\xc3\x9c\x83w\x1e\xd4t\xb6\x90\xb0\x10f\xda\xd5F@f"$\x12\x89\x7fN\xf74\x86\xcf^\xf3/\xbc\x14\xea\xc4\x88w\x04\x0bP\xe4\xa8bL\x95Z)\xf8\x9f\x87\t\x14iR,\x0e\x8e\xdc\xd1\xce^\xc3U'  # initial challenge (randomly generated m bytes array)
else:
    I = os.urandom(M)

def memory_build(I,n,P,H):
    # Step (1)
    # Building a challenge dependent memory
        
    # Step (1.a)
    for p in range(P):
        for i in range(n):
            # TODO: properly turn the int p into a 4 bytes hex string
            hash_input = str(p) + I
    
            H.update(hash_input)
            X[p*l+i] = H.digest()
    
    # Step (1.b)
    for p in range(P):
        for i in range(n,l):
            # computing phi_{k}(i) for all k up until n
            phis_i = phis(i,n)
    
            # building the input of the hash function
            hash_input = ""
            for phi in phis_i:
                hash_input += X[p*l+ phi]
            H.update(hash_input)
            
            # inserting the computed hash in the array
            X[p*l+i] = H.digest()

    return X

print(memory_build(I,n,P,H))

