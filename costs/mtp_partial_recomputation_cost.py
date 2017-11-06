#! /usr/bin/env python3
#
# recomputation cost for all cells depending on #dependencies
#  - t>0: if 1/t are available
#  - t<0: 1/-t are removed
#  - distribution:
#    alt (regularly spaced), fist & last, oq/eq (odd/even quarters)
#
# if tmp:
#   space is available for keeping recomputed values, so that if the
#   same value is needed twice it is reused directly. However it seems
#   that the reuse rate is pretty low, so this does not change cost much.
# else:
#   This does assume that a temporary memory is NOT allocated to store
#   intermediate results, so all deps must be recomputed if required
#   for distinct array elements on the same recomputation of an element.
#
# NOTE
#   for very small p and many segments, the memory for one segment could
#   be allocated and the intermediate values kept, so that the recomputation
#   cost is only 2^p / 2.
#

# > 64 multiplier for small p and n
# p, t, n, N, c_F = 13, 2, 4, 100, 1
# p=13 t=2 n=4 N=1000 c_F=1

# defaults, see "-0"
p, t, n, N, c_F, tmp, bias = 21, 2, 2, 1, 11, False, 2

store = 'alt'
debug = False
seed = None

import re
is_direct_eval = re.compile('\w+=([-+]?\d|true|false|none)', re.I).match

# handle command line options
import sys
for a in sys.argv[1:]:
    if a == '-0': # MTP-Argon2 proposal
        p, t, n, N, c_F, tmp, bias = 21, 2, 2, 1, 11, False, 2
    elif a == '-1': # report proposal
        p, t, n, N, c_F, tmp, bias = 25, 2, 2, 1, 1, False, 2
    elif a == '-x': # testing...
        p, t, n, N, c_F, tmp, bias = 13, 2, 6, 5, 1, True, 2
    elif a == '-d':
        debug = True
    elif a == '-a':
        store = 'alt'
    elif a == '-f':
        store = 'first'
    elif a == '-l':
        store = 'last'
    elif is_direct_eval(a):
        exec(a)
    else: # put quotes on a string
        exec(re.sub("=(.*)",r"='\1'", a))

# cleanup phi & cost caches
def reset_cache():
    global phi_x, cost_F
    phi_x = [-1 for i in range(0, T)]
    phi_x[0:n] = [None] * n
    cost_F = [-1 for i in range(0, T)]
    for i in range(n):
        cost_F[i] = 0
    # precompute phis for determinism
    for i in range(n, T):
        phis(i, n)

# possibly biased phi random function implementation
import random
if seed is not None:
    assert N == 1, "only one run with a seed"
    random.seed(seed)

# biased random generator
def phi(i, bias=bias):
    assert 0 <= i < T
    if phi_x[i] == -1:
        # get a possibly biased number in [0,1)
        u = random.random()
        if bias == 2:
            r = 1.0 - u * u
        elif bias == 1:
            r = 1.0 - u
        elif bias == 3:
            r = 1.0 - u * u * u
        else: # handle any power...
            r = 1.0 - pow(u, bias)
        # compute corresponding index
        phi_x[i] = int((i - 1) * r)
        assert 0 <= phi_x[i] and phi_x[i] <= i - 1
        # (2 if bias == 2 else 1)
    return phi_x[i]

# \phi_k functions for k \in [0, 11]
PHI_K = [
    lambda i, pi, n: i-1,
    lambda i, pi, n: pi,
    lambda i, pi, n: (pi + i) // 2,
    lambda i, pi, n: 7 * i // 8,
    lambda i, pi, n: (pi + 3*i) // 4,
    lambda i, pi, n: (3*pi + i) // 4,
    # n > 6
    lambda i, pi, n: 7 * pi // 8,
    lambda i, pi, n: 3 * pi // 4,
    lambda i, pi, n: 3 * (i - 1) // 4,
#    lambda i, pi, n: pi // 4,
#    lambda i, pi, n: (i - 1) // 4,
]

# return a set of needed elements, up to n
def phis(i, n, bias=bias):
    assert 1 <= n and n <= len(PHI_K)
    phi_i = phi(i, bias) if n >= 2 else None
    return [ phi(i, phi_i, n) for phi in PHI_K[:n] ]

def isStored(i, store):
    assert 0 <= i < T
    if store == 'alt':
        return \
            (t > 0 and i % t == 0 or i < 2) or \
            (t < 0 and i % -t != 0 or i < 2)
    elif store == 'first':
        # first elements
        return \
            (t > 0 and i * t < T) or \
            (t < 0 and i * -t < (-t - 1) * T)
    elif store == 'oq':
        # odd quarters
        return \
            (t > 0 and 2 * i * t < T) or \
            (t > 0 and 2 * i >= T and 2 * (i - T/2) * t < T) or \
            (t < 0 and 2 * i * -t < (-t - 1) * T) or \
            (t < 0 and 2 * i >= T and 2 * (i - T/2) * -t < (-t - 1) * T)
    elif store == 'last':
        # last elements are available
        return \
            (t > 0 and i * t <= (t-1) * T) or \
            (t < 0 and i * -t >= T)
    elif store == 'oq':
        raise Exception("store='oq' not implemented yet")
    elif store == 'eq':
        raise Exception("store='eq' not implemented yet")
    else:
        raise Exception("invalid store=%s" % store)


# cost in F calls of accessing X[i] (no temporary memory version)
def costFnt(i, store):
    #if debug:
    #    print("costFnt(%d) %s" % (i, phis(i, n)))
    assert 0 <= i < T
    # we need to evaluate
    if cost_F[i] == -1:
        if isStored(i, store):
            cost_F[i] = 0
        else:
            # first call, record computed items to avoid them
            # recursive evaluation
            cost_F[i] = 1 + sum(costFnt(j, store) for j in phis(i, n))
            # not available, recompute
            assert cost_F[i] != -1, "some cost"
    return cost_F[i]

# recursive evaluation in F calls for computing X[i] with temporary memory
def costFt(i, store, done=None):
    start = False
    if done is None:
        start = True
        # already computed, use memoized value
        if cost_F[i] >= 0:
            return cost_F[i]
        # else compute recursively
        import numpy as np
        done = np.zeros(i+1, dtype=np.bool)
        # VERY INEFFICIENT: [ False for x in range(i+1) ]
    if done[i]: # déjà vu
        if debug:
            print("# already computed: %d" % i)
        return 0
    if cost_F[i] == -1 and isStored(i, store):
        cost_F[i] = 0
    if cost_F[i] == 0:
        return 0
    else:
        #print("cost_F[%d] = %d" % (i, cost_F[i]))
        assert not done[i]
        done[i] = True
        cost = 1 + sum(costFt(j, store, done) for j in phis(i, n))
        if debug:
            print("cost(%d) = %d %s (%s)" %
                  (i, cost, phis(i, n),
                   ''.join('*' if done[j] else '.' for j in range(len(done)))))
        if start:
            cost_F[i] = cost
        return cost

# cost in F calls of accessing X[i]
def costF(i, store, tmp):
    if tmp: # with tmp space
        return costFt(i, store)
    else: # without tmp space
        return costFnt(i, store)

# cost in X accesses of accessing X[i]
def costX(i):
    if t > 0 and i % t == 0 or i < 2: # available
        return 1
    elif t < 0 and i % -t != 0 or i < 2: # available
        return 1
    else: # recompute... hmmm, significant collisions?
        return n * costF(i, store, tmp)

# total cost accumulated over all cells
def totals(T):
    reset_cache()
    totF, totX = 0, 0
    for i in range(n, T):
        if debug:
            print("* computing i=%d" % i)
        totF += costF(i, store, tmp)
        totX += costX(i)
    return (totF, totX)

#
# COMPUTE AVERAGE COSTS
#

T = 2 ** p

print("p=%d t=%d n=%d N=%d store=%s tmp=%s bias=%d seed=%s" %
      (p, t, n, N, store, tmp, bias, seed))

if N == 1:
    (totF, totX) = totals(T)
else:
    print("computing: ", file=sys.stderr, end='', flush=True)
    assert N > 1
    totF, totX = 0.0, 0.0
    for i in range(N):
        (tF, tX) = totals(T)
        totF += tF
        totX += tX
        print("*", file=sys.stderr, end='', flush=True)
    totF /= N
    totX /= N
    print("", file=sys.stderr)

# compute saving ratio alpha
alpha = ((t-1.0) / t) if t > 0 else (1.0 / -t)
# but at least n to recompute a value!
if alpha == 1.0:
    alpha = 1.0 - n / T

# show parameters and result
print("p=%d (T=2^p=%d) t=%s (α=%.3f) n=%d N=%d c_F=%d store=%s tmp=%s bias=%d seed=%s" %
      (p, T, t, alpha, n, N, c_F, store, tmp, bias, seed))

if N > 1:
    print("average over %d" % N)

try:
    print("F calls per cell: %f%s" % ((totF / T), " ?" if totF >= T*T else ""))
    print("cost multiplier: %f" % (c_F * totF / T))
    print("X accesses per cell: %f" % (totX / T))
except OverflowError:
    print("F calls per cell: ~ 10^%d" % (len(str(totF)) - len(str(T)) - 1))
    print("X accesses per cell: ~ 10^%d" % (len(str(totX)) - len(str(T)) - 1))

if (debug):
    print("F = %s" % cost_F)
