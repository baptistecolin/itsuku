#! /usr/bin/env python3

import os
from itsuku import *

I = os.urandom(8)
T = 64
n = 4
l = 32
x = 8
L = 4
S = x
M = 3

# for information
# P = (T + (l-1)) // l

# expect 2**12 = 4096 iterations
d = b'\x00' + b'\x0f' + b'\xff' * (S-2)

print("I=%s T=%d l=%d n=%d M=%d L=%d S=%d x=%d d=%s" % (I.hex(), T, l, n, M, L, S, x, d.hex()))

print("solving...")
pow, Omega, cnt = solvePoW(I, T, l, n, M, L, S, x, d)

print("found Omega=%s after %d interations" % (Omega.hex(), cnt))
print("POW is %s" % pow)

print("checking...")
ok, nOmega = checkPoW(I, T, l, n, M, L, S, x, d, pow)

print("ok=%s recomputed Omega=%s" % (ok, nOmega.hex()))

assert Omega == nOmega
assert ok is True
