#!/usr/bin/env python3

import phis
from math import log2, floor

bias = 2 # quadratic bias
n = 2 # number of dependencies
T = 2**10 # length of the main array
x = 64 # size of elements in the main array
L = floor(3.3*log2(T)) # length of one search
d = 5 # difficulty of the PoW
P = 1 # number of independent sequences