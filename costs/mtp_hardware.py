#! /usr/bin/env python3
#
# $Id: mtp_hardware.py 2046 2017-11-23 15:30:39Z fabien $
#

import sys
import math
from math import ceil

# units
k, M, G  = 1000, 1000**2, 1000**3
Ki, Mi, Gi = 1024, 1024**2, 1024**3

# convergence
epsilon = 0.000001

# gate-equivalent areas
Ha = 100 * k # BLAKE2 core
F = 30 * M # BLAKE2 core frequency in Hz
Ma = 8 # 4T SRAM
tr = 4 # transistors per GE

for a in sys.argv[1:]:
    if a in ['-0', '-0a', '-0b' ]: # MTP-Argon2
        algo = 'MTP-Argon2'
        T, x, n, S, L, theta, bs = 2**21, 1024, 2, 16, 70, 9, False
        cX, cF, cR, cDN = 9, 11, 1.0, 68.4
        if a == '-0b':
            cDN = 116.2
    elif a == '-1': # prop1
        algo = 'prop1'
        T, x, n, S, L, theta, bs = 2**25, 64, 0, 64, 16, 9, True
        cX, cF, cR, cDN = 1, math.inf, math.inf, math.inf
    elif a == '-2': # prop2
        algo = 'prop2'
        T, x, n, S, L, theta, bs = 2**25, 64, 2, 64, 84, 9, True
        cX, cF, cR, cDN = 1, 1, 1.0, 98.0
    elif a == '-3': # prop3
        algo = 'Itsuku'
        T, x, n, S, L, theta, bs= 2**25, 64, 4, 64, 9, 9, True
        cX, cF, cR, cDN = 1, 1, 278.4, None
    elif a == '-21.1' or a == '-V100':
        hw, Ga, Bw = 'V100', 21.1 * G / tr, 1000.0 * G
    elif a == '-30':
        hw, Ga, Bw = '30bn', 30.0 * G / tr, 1300.0 * G
    elif a == '-50':
        hw, Ga, Bw = '50bn', 50.0 * G / tr, 2000.0 * G
    elif a == '-75':
        hw, Ga, Bw = '75bn', 75.0 * G / tr, 2500.0 * G
    elif a == '-100':
        hw, Ga, Bw = '100bn', 100.0 * G / tr, 3000.0 * G
    else:
        raise Exception("not implemented yet: %s" % a)

array_size = T * x

# (average) search state size
# 11 for N, 1 for L, hash size S
SS = (12 + S) * (((L + 1) // 2) if bs else 1)

print("S=%d SS=%d" % (S, SS))

def unit(v):
    if v is None:
        return '-'
    if v > 1000 * G:
        return "%.2fT" % (v / (1000 * G))
    elif v > G:
        return "%.2fG" % (v / G)
    elif v > M:
        return "%.2fM" % (v / M)
    elif v > G:
        return "%.2fk" % (v / k)
    else:
        return "%.2f" % v

##
## Direct method
##

# number of core for a fully pipelined PoW solver
C0 = 1 + cX * L + (((L + 1) // 2) if bs else 0)
Ca = C0 * Ha + theta * L * SS * Ma

# maximum number of solvers on die
Ndie = Ga / Ca
# maximum number of pipelined PoW solvers for the bandwidth
# print("Bw=%s x=%s L=%s F=%s" % (Bw, x, L, F))
Nbw = Bw / (x * L * F)

assert Nbw < Ndie, "bandwidth limited"

if array_size * Ma < Ga:
    # array in SRAM, threads not nedded
    imp = 'Array in SRAM'
    N = (Ga - array_size * Ma) / (C0 * Ha)
    Nbw, hit = None, None
else: # array in DRAM + on die cache
    # find fix point
    imp = 'Array in DRAM + cache'
    hit = 1.0
    new_hit = 0.0
    while abs(hit - new_hit) > epsilon:
        hit = new_hit
        N = Bw / ((1.0 - new_hit) * x * L * F)
        new_hit = (Ga - N * Ca) / (T * x * Ma)

print("%s on %s (%s): C=%d Ndie=%s Nbw=%s N=%.2f throuput=%s (%.1fM)" %
      (algo, hw, imp, C0, unit(Ndie), unit(Nbw), N, unit(N * F), N * F / M))

##
## Half Array
##

CR = C0 + cR * cF * L
nloads = (n-1) * cR + 1
CRa = CR * Ha + nloads * theta * L * SS * Ma

#print("cR=%d nloads=%d" % (cR, nloads))

NRdie = Ga / CRa
# maximum number of pipelined PoW solvers for the bandwidth
# print("Bw=%s x=%s L=%s F=%s" % (Bw, x, L, F))
NRbw = Bw / (nloads * x * L * F)

if 0.5 * array_size * Ma < Ga:
    # half array in SRAM, threads not really nedded
    imp = 'Half array in SRAM'
    NR = (Ga - 0.5 * array_size * Ma) / (CR * Ha)
    NRbw = None
else: # array in DRAM + on die cache
    # find fix point
    imp = 'Half array in DRAM + cache'
    hit = 1.0
    new_hit = 0.0
    while abs(hit - new_hit) > epsilon:
        hit = new_hit
        NR = Bw / ((1.0 - new_hit) * nloads * x * L * F)
        new_hit = (Ga - NR * CRa) / (0.5 * array_size * Ma)

print("%s on %s (%s): C=%.0f Ndie=%s Nbw=%s N=%.2f throuput=%s (%.1fM)" %
      (algo, hw, imp,
       CR, unit(NRdie), unit(NRbw), NR, unit(NR * F), NR * F / M))

##
## Transposed Search
##

imp = 'Transposed search in SRAM'

#CT = C0
#CTa = C0 * Ha # no threads
#NTdie = Ga / CTa

# number of parallel searches
nps = Ga // ((SS + 4) * T * Ma)
#print("nps=%d" % nps)
if nps > 0:
    # number of array elements transfered per second
    nTbw = Bw / x
    # number of needed cores to process the searches
    nTx = nps * cX * nTbw / F
    CT = ceil(nTx / (cX * L) * C0)
    # check that there is room enough for the hash cores...
    while CT * Ha + nps * (SS+4) * T * Ma > Ga:
        print("decrementing number of parallel searches")
        nps -= 1
        nTx = nps * cX * nTbw / F
        CT = ceil(nTx / (cX * L) * C0)
    assert CT * Ha + nps * (SS+4) * T * Ma <= Ga
    NTc = CT / C0
else:
    CT, NTc = 0, 0

imp += " [nps=%d]" % nps
print("%s on %s (%s): C=%.0f Ndie=%s Nbw=%s N=%.2f throuput=%s (%.1fM)" %
      (algo, hw, imp,
       CT, unit(None), unit(None), NTc, unit(NTc * F), NTc * F / M))


# alternatively, send many search states from DRAM...
# which is assume large enough so that we can neglect the array elements
# tranfers, which are amortized over a large number of searches
imp = 'Transposed search in DRAM'

# NOTE in & out (& in)
# NOTE the external DRAM storage is much larger with bs!
mu = 3 if bs else 2
# number of searches for which states elements transfered per second
nTbw = Bw / (mu * (SS + 4))
# number of cores needed to process the corresponding searches
nTx = cX * nTbw / F
# total number of needed cores including other stuff (init)
CT = ceil(nTx / (cX * L) * C0)
# production per tick
NTc = CT / C0
# check that there is room enough for the hash cores...
assert CT * Ha <= Ga, "large enough for needed cores"
print("die usage is %.1f%%" % (100.0 * CT * Ha / Ga))

print("%s on %s (%s): C=%.0f Ndie=%s Nbw=%s N=%.2f throuput=%s (%.1fM)" %
      (algo, hw, imp,
       CT, unit(None), unit(None), NTc, unit(NTc * F), NTc * F / M))

##
## Dinur-Nadler Attack
##

imp = "Dinur-Nadler attack [cDN=%s]" % cDN

if cDN is not None:
    CDN = cDN * C0
    NDN = Ga / (CDN * Ha)
else:
    CDN, NDN = None, 0.0

print("%s on %s (%s): C=%.0f Ndie=%s Nbw=%s N=%.2f throuput=%s (%.1fM)" %
      (algo, hw, imp,
       CDN if CDN is not None else 0,
       unit(NDN), unit(None), NDN, unit(NDN * F), NDN * F / M))

##
## pseudo random array attack, if applies
##
imp = "PRA attack"
CRA = C0 + cF * L
NRA = Ga / ( CRA * Ha)
print("%s on %s (%s): C=%0.f Ndie=%s throuput=%s (%.1fM)"
      % (algo, hw, imp, CRA, unit(NRA), unit(NRA * F), NRA * F / M))
