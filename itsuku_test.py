import pytest
from itsuku import *

def test_phi():
    # it shoud fail if the seed is not 4 bytes long
    seed = int_to_4bytes(123)
    seed_3bytes = seed[:3]
    seed_8bytes = seed+seed

    assert len(seed_3bytes) == 3
    assert len(seed_8bytes) == 8

    with pytest.raises(AssertionError):
        phi(seed_3bytes, 4)
    with pytest.raises(AssertionError):
        phi(seed_8bytes, 4)

    # it should return the same result using the high-level or low level algorithm
    assert phi(seed, 4, method='high-level') == phi(seed, 4, method='low-level')
    assert phi(seed, 28, method='high-level') == phi(seed, 28, method='low-level')
    assert phi(seed, 1024, method='high-level') == phi(seed, 1024, method='low-level')

    # it should return a value phi(i) inferior to i
    assert phi(seed, 2) < 2
    assert phi(seed, 4) < 4
    assert phi(seed, 256) < 256
    assert phi(seed, 1024) < 1024


def test_phis():
    seed = int_to_4bytes(256)

    # it should output an array of length n
    for n in range(1,12):
        assert len(phis(seed, 10, n)) == n

    # TODO : more tests, probably

def test_H():
    # it should return a bytes array of length M
    x = int_to_4bytes(123456)
    for i in range(1,15):
        assert len(H(i,x)) == i
    
    # the H function should output the last M bytes of the sha512 hash of i
    for i in range(1,15):
        sha = sha512() # resetting the sha function
        sha_input = int_to_4bytes(i)
        sha.update(sha_input)
        assert sha.digest()[:10] == H(10,sha_input)
        


def test_int_to_4bytes():
    # it should always return a 4 bytes string
    assert len(int_to_4bytes(0)) == 4
    assert len(int_to_4bytes(1)) == 4
    assert len(int_to_4bytes(10)) == 4
    assert len(int_to_4bytes(1024)) == 4
    assert len(int_to_4bytes(4294967295)) == 4
    assert isinstance(int_to_4bytes(0), bytes)

    # it should return the expected hex bytes string, in a big endian order 
    assert int_to_4bytes(0) == b"\x00\x00\x00\x00"
    assert int_to_4bytes(1) == b"\x00\x00\x00\x01"
    assert int_to_4bytes(2) == b"\x00\x00\x00\x02"
    assert int_to_4bytes(16) == b"\x00\x00\x00\x10"
    assert int_to_4bytes(256) == b"\x00\x00\x01\x00"

def test_memory_build():
    M = 64
    x = 64
    T = 2**5
    I = os.urandom(M)

    for P in [1,2,4]:
        l = T//P
        for n in range(2,min(12,l)): # it should work for different values of n. n can't get bigger than l, otherwise the n "seeds" cannot fit in a slice
            X = memory_build(I, T, n, P, M)

            # Initialization steps
            for p in range(P):
                for i in range(n):
                    hash_input = int_to_4bytes(i) + int_to_4bytes(p) + I
                    assert X[p*l+i] == H(x, hash_input)

            # Construction steps
            for p in range(P):
                for i in range(n,l):
                    seed = X[p*l+i-1][:4] # seed that is used at each step is the 
                                          # 4 first bytes of the previous array item
                    
                    # asserting that the 0<=phi_k(i)<i condition is actually verified
                    phi_k = phis(seed,i,n)
                    for phi_k_i in phi_k:
                        assert phi_k_i < i
                        assert 0 < phi_k_i

                    hash_input = b""
                    for phi_k_i in phi_k:
                        hash_input += X[p*l+phi_k_i]
                    
                    # asserting the validity of the constructed item
                    assert X[p*l+i] == H(x, hash_input)

def test_merkle_tree():
    M = 64
    T = 2**5
    n = 2
    P = 1
    I = os.urandom(M)
    l = T//P
    
    for n in range(2,12): # should work for different values of n
        X = memory_build(I, T, n, P, M)
        
        MT = merkle_tree(I, X, M)

        # asserting the length is 2*T-1
        assert len(MT) == 2*T-1
        # asserting the end of the MT is actually the hashed original array
        assert MT[-T:] == [H(M,x) for x in X]
        # asserting the constructed items are the hash of their sons
        for i in range(T-1):
            assert MT[i] == H(M, MT[2*i+1]+MT[2*i+2]+I)
    
    # test on a particular case : if the initial array is constant,
    # then each "floor" of the merkle tree should be constant
    X0 = [int_to_4bytes(0)]*T
    MT0 = merkle_tree(I, X0, M)

    # iterating over the floors
    for i in range(1,int(log(T,2))):
        #remembering the value that we expect to find all over the floor
        value = MT0[(2**i)-1]
        # iterating inside the floor
        for j in range((2**i)-1, (2**(i-1))-2):
            assert MT0[i] == value

def test_xor():
    assert xor(b"\x00", b"\x00") == b"\x00"
    assert xor(b"\x01", b"\x00") == b"\x01"
    assert xor(b"\x00", b"\x01") == b"\x01"
    assert xor(b"\x01", b"\x01") == b"\x00"

def test_compute_Y():
    M = 64
    T = 2**5
    P = 1
    S = 64
    L = ceil(3.3*log(T,2))
    I = os.urandom(M)
    l = T//P
    
    for n in range(2,12): # should work for different values of n
        X = memory_build(I, T, n, P, M)
        MT = merkle_tree(I, X, M)
        PSI = MT[0]
        N = os.urandom(32) # nounce

        Y, OMEGA, i = compute_Y(I, X, L, S, N, PSI)

        # asserting length
        assert len(Y) == L+1
        assert len(i) == L
        # verifying Y[0] is built as expected
        assert Y[0] == H(S, N + PSI + I)
        # verifying Y is correctly constructed
        for j in range(1,L+1):
            assert i[j-1] == int.from_bytes(Y[j-1], 'big') % T
            assert Y[j] == H(S, Y[j-1] + xor(X[i[j-1]], I))

def test_is_PoW_solved():
    assert is_PoW_solved(b'\x00'*64, b'\x00'*63 + b'\x01') == True
    assert is_PoW_solved(b'\x00'*64, b'\x00'*64) == False
    assert is_PoW_solved(b'\xff'*63 + b'\xfe', b'\xff'*64) == True

    with pytest.raises(AssertionError):
        is_PoW_solved(b'\x00', b'\x00'*64)
    with pytest.raises(AssertionError):
        is_PoW_solved(b'\x00', b'\x00')
    with pytest.raises(AssertionError):
        is_PoW_solved(b'\x00'*64, b'\x00')

def test_build_L():
    M = 64
    T = 2**5
    P = 1
    S = 64
    L = ceil(3.3*log(T,2))
    I = os.urandom(M)
    l = T//P
    
    for n in range(2,12): # should work for different values of n
        X = memory_build(I, T, n, P, M)
        MT = merkle_tree(I, X, M)
        PSI = MT[0]
        N = os.urandom(32) # nounce
        Y, OMEGA, i = compute_Y(I, X, L, S, N, PSI)

        round_L = build_L(i, X, P, n)
        
        for i_j in i:
            assert len(round_L[i_j]) == n
            p = i_j // l
            if i_j % l < n:
                # assert correct construction
                assert round_L[i_j] == X[p*l:p*l+n]
                
                # by construct, X[i_j] should be part of round_L[i_j]
                assert X[i_j] in round_L[i_j]

                # assert that the elements of round_L are actually computable
                for k in range(0,n):
                    stuff_to_hash = int_to_4bytes(k) + int_to_4bytes(p) + I
                    assert round_L[i_j][k] == H(M, stuff_to_hash)

            else:
                seed = X[i_j-1][:4]
                # assert correct construction
                assert round_L[i_j] == [ X[p*l + phi_k_i] for phi_k_i in phis(seed, i_j, n) ]

                # make sure X[i_j] can be recomputed using the elements of round_L[i_j]
                hash_input = b""
                for item in round_L[i_j]:
                    hash_input += item
                assert H(M,hash_input) == X[i_j]

def test_provided_indexes():
    M = 64
    T = 2**5
    P = 1
    S = 64
    L = ceil(3.3*log(T,2))
    I = os.urandom(M)
    l = T//P
    
    for n in range(2,12): # should work for different values of n
        X = memory_build(I, T, n, P, M)
        MT = merkle_tree(I, X, M)
        PSI = MT[0]
        N = os.urandom(32) # nounce
        Y, OMEGA, i = compute_Y(I, X, L, S, N, PSI)
        round_L = build_L(i, X, P, n)

        indexes = provided_indexes(round_L, P, T, n)

        for i_j in i:
            # all the i_j should be in indexes as X[i_j] can alwas be recomputed
            # using only the elements of round_L[i_j]. Thus, X[i_j] can and must
            # be considered as a given if round_L is known
            assert i_j in indexes
            
            p = i_j // l
            
            if i_j % l < n :
                # case where the elements of round[i_j] were computed at step (1.a)
                for k in range(p*l, p*l+n):
                    assert k in indexes
            else :
                # case where the elements of round[i_j] were computed at step (1.b)

                seed = round_L[i_j][0][:4]

                # The seed computed using round_L should be the same seed as the
                # one used during the construction of X
                assert seed == X[i_j - 1][:4]
                for index in [p*l+phi_k_i for phi_k_i in phis(seed, i_j, n)]:
                    assert index in indexes
            
        # Asserting there are no duplicates
        assert len(indexes) == len(set(indexes))

def test_build_Z():
    M = 64
    T = 2**5
    P = 1
    S = 64
    L = ceil(3.3*log(T,2))
    I = os.urandom(M)
    l = T//P
    
    for n in range(2,12): # should work for different values of n
        X = memory_build(I, T, n, P, M)
        MT = merkle_tree(I, X, M)
        PSI = MT[0]
        N = os.urandom(32) # nounce
        Y, OMEGA, i = compute_Y(I, X, L, S, N, PSI)
        round_L = build_L(i, X, P, n)
        
        # A shift has to be applied so the indexes match those of the
        # Merkle Tree and not those of X.
        indexes = [ index + T - 1 for index in  provided_indexes(round_L, P, T, n)]

        Z = build_Z(round_L, MT, P, T, n)

        for k in Z:
            assert k not in indexes
            assert Z[k] == MT[k]
            if k >= T-1:
                assert Z[k] == H( M, X[k-(T-1)] )
        
        assert set(Z.keys()) == set(opening(T, provided_indexes(round_L, P, T, n)))


@pytest.mark.skip(reason="to be filled")
def test_build_JSON_output():
    # TODO : write test
    return None

@pytest.mark.skip(reason="to be filled")
def test_PoW():
    # TODO : write test
    return None
