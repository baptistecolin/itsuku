#!/usr/bin/env python

import phis
from pyblake2 import blake2b
from hashlib import sha512
from math import ceil, log

bias = 2 # quadratic bias
n = 2 # number of dependencies
T = 2**10 # length of the main array
x = 64 # size of elements in the main array
M = 64 # size of elemets in the Merkel Tree
L = ceil(3.3*log(T,2)) # length of one search
d = 5 # difficulty of the PoW
P = 1 # number of independent sequences
l=T/P # length of one independent sequence

X = [] # memory
H = sha512() # hash function

if M == 64:
    I = 'N\x8a\xc3\x9c\x83w\x1e\xd4t\xb6\x90\xb0\x10f\xda\xd5F@f"$\x12\x89\x7fN\xf74\x86\xcf^\xf3/\xbc\x14\xea\xc4\x88w\x04\x0bP\xe4\xa8bL\x95Z)\xf8\x9f\x87\t\x14iR,\x0e\x8e\xdc\xd1\xce^\xc3U'  # initial challenge (randomly generated m bytes array)
else:
    I = os.urandom(M)

