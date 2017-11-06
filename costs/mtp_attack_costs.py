#! /usr/bin/env python3
#
# $Id: mtp_attack_costs.py 1893 2017-10-14 20:20:01Z fabien $
#
# expected number of H calls to process when looking for a special X sequence
#

# defaults values
p = 21
L = 70
c_F = 11
c_X = 9
back_sweep = False
with_precomp = False
tmax = 20
n = 2 # how many dependencies, including previous element
dmax = 70 # max d

# get arguments
import sys

# -s & cold-blooded settings
for a in sys.argv[1:]:
    if a == '-h' or a == '--help':
        print("%s [-hs01] [var=val]*" % sys.argv[0], file=sys.stderr)
        sys.exit(0)
    elif a == '-s':
        back_sweep = True
    elif a == '-0':
        # defaults
        p, L, c_F, c_X, back_sweep, with_precomp, n, dmax = \
            21, 70, 11, 9, False, False, 2, 100
    elif a == '-1':
        # new defaults
        p, L, c_F, c_X, back_sweep, with_precomp, n, dmax = \
            25, 31, 1, 1, True, True, 6, 100
    else:
        exec(a)

assert n >= 2, "attack requires at least one other dependency"
assert dmax > 15

from math import log2

# define T unless explicitely
try:
    p = log2(T)
except NameError:
    T = 2 ** p

# helper functions
def sequence_weighted_cost(t, n):
    return sum(i * (1.0 - t ** (1-n)) * (t ** -((n-1)*(i-1))) \
               for i in range(t-2, 0, -1)) + (t-2) * t ** -((n-1)*(t-2))

def search_weighted_cost(t, L, bs):
    return sum((1 + c_X * i + 0.5 * c_F * t * i) * ((t-1)/t) ** i * (1/t)
               for i in range(0, L)) + \
                   ( 1 + c_X * L + 0.5 * c_F * t * L + \
                     ((0.5 * L) if back_sweep else 0)) * ((t-1)/t) ** L

print("p=%d (T=2^%d=%d) L=%d c_F=%d c_X=%d back_sweep=%s with_precomp=%s tmax=%d n=%d dmax=%d" %
      (p, p, T, L, c_F, c_X, back_sweep, with_precomp, tmax, n, dmax))

# non cheating (nc) cost of computing X
sq_nc = c_F * T

# non cheating (nc) cost of computing an omega
se_nc = 1 + c_X * L + ((0.5 * L) if back_sweep else 0)

print("non cheating cost: 2^%.2f + 2^{d+%.2f} | 2^{%.2f}" %
      (log2(sq_nc), log2(se_nc), log2(sq_nc + se_nc * 2**dmax)))

print(" t | DN-FC & DN-DN | TOTAL for dmax (multiplier)")
for t in range(3, tmax + 1):
    # number of F computation per attempt
    sq_wc = sequence_weighted_cost(t, n)
    # sequence expected number of attempts
    sq_na = t ** ((n-1)*(t-2))
    # log2 cost of computing one sequence
    sq_lc = log2(c_F * sq_wc * sq_na)
    # number of needed sequences
    sq_nb = T / t
    # total cost
    sq_ltc = c_F * sq_wc * sq_na * sq_nb
    # cost in H calls of one search with t compression
    se_wc = search_weighted_cost(t, L, back_sweep)
    # number of expected attempts for ONE omega,
    # because of missing array X elements...
    se_na = (t/(t-1)) ** L
    # Dinur & Nadler approximated evaluations of their attack, with n...
    sq_dn = c_F * T * t ** ((n-1)*(t-2)-1)
    se_dn = c_F * 0.5 * t * t * (1-1/t) ** -L
    # total cost
    total = (sq_ltc if with_precomp else 0) + se_wc * se_na * 2 ** dmax
    # show result
    print("%2d | 2^%.2f + 2^{d+%.2f} | 2^%.2f + 2^{d+%.2f} | 2^{%.2f} (x %.1f)" %
          ( t,
            # my Dinur & Nadler cost evaluation
            log2(sq_ltc), log2(se_wc*se_na),
            # Dinur & Nadler evaluations
            log2(sq_dn), log2(se_dn),
            # total cost for dmax
            log2(total),
            total / ((sq_nc if with_precomp else 0) + se_nc * 2 ** dmax)
          ))
