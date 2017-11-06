#! /usr/bin/env python3

from math import ceil

f = 0.5
#f = 0.75
for L in range(1, 16):
    print("m(%f,%d) = %s" %
          (f, L, (ceil(0.5 * L) + (f**-L - f) / (1 - f)) / (ceil(1.5 * L) + 1)))
