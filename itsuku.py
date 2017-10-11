#!/usr/bin/env python3

import phis
from pyblake2 import blake2b
from math import floor, log

bias = 2 # quadratic bias
n = 2 # number of dependencies
T = 2**10 # length of the main array
x = 64 # size of elements in the main array
L = floor(3.3*log(T,2)) # length of one search
d = 5 # difficulty of the PoW
P = 1 # number of independent sequences

