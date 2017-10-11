import random

bias = 2

# biased random generator
def phi(i, bias=bias):
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
    result = int((i - 1) * r)
    return result

# return a set of needed elements, up to n
def phis(i, n, bias=bias):
    assert n >= 1 and n <= 11
    phis_i = { (i-1) }
    if n >= 2:
        phi_i = phi(i)
        # set has -=, but not += for adding an element
        # |= is for sets
        phis_i.add(phi_i)
        # YUK! I want a nice switch...
        if n >= 3:
            phis_i.add(phi_i // 2)
            if n >= 4:
                phis_i.add((i - 1) // 2)
                if n >= 5:
                    phis_i.add((phi_i + i) // 2)
                    if n >= 6:
                        phis_i.add(3 * phi_i // 4)
                        if n >= 7:
                            phis_i.add(3 * i // 4)
                            if n >= 8:
                                phis_i.add(phi_i // 4)
                                if n >= 9:
                                    phis_i.add(i // 4)
                                    if n >= 10:
                                        phis_i.add(phi_i * 7 // 8)
                                        if n >= 11:
                                            phis_i.add(i * 7 // 8)
    return phis_i
