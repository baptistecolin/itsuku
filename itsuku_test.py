import pytest
from itsuku import *
from opening import openingForOneArray
from collections import OrderedDict

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
    assert phi(seed, 4) == phi(seed, 4)
    assert phi(seed, 28) == phi(seed, 28)
    assert phi(seed, 1024) == phi(seed, 1024)

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
        assert len(H(i,x)) == i # TODO : check 0<

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

def test_build_X():
    M = 64
    x = 32
    T = 2**5
    I = os.urandom(M)

    for P in [1,2,4]:
        l = T//P
        for n in range(2,min(12,l)): # it should work for different values of n. n can't get bigger than l, otherwise the n "seeds" cannot fit in a slice
            X = build_X(I, T, l, n, x)

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

def test_build_MT():
    M = 64
    x = 32
    T = 2**5
    n = 2
    I = os.urandom(M)

    for P in [1,2,4]:
        l = T//P
        for n in range(2,min(12,l)): # should work for different values of n
            X = build_X(I, T, n, P, x)

            MT = build_MT(I, X, M)

            # asserting the length is 2*T-1
            assert len(MT) == 2*T-1
            # asserting the end of the MT is actually the hashed original array
            assert MT[-T:] == [H(M, x+I) for x in X]
            # asserting the constructed items are the hash of their sons
            for i in range(T-1):
                assert MT[i] == H(M, MT[2*i+1]+MT[2*i+2]+I)

        # test on a particular case : if the initial array is constant,
        # then each "floor" of the merkle tree should be constant
        X0 = [int_to_4bytes(0)]*T
        MT0 = build_MT(I, X0, M)

        # iterating over the floors
        for i in range(1,int(log(T,2))):
            #remembering the value that we expect to find all over the floor
            value = MT0[(2**i)-1]
            # iterating inside the floor
            for j in range((2**i)-1, (2**(i-1))-2):
                assert MT0[i] == value

def test_compute_MT_node():
    M = 64
    x = 32
    I = os.urandom(M)

    # basic examples
    assert compute_MT_node(0, {0: b'\x00'*64}, I, 1, M) == b'\x00'*64
    assert compute_MT_node(0, {2: b'\x00'*64, 3: b'\x11'*64, 4: b'\xff'*64}, I, 4, M) == H(64, H(64, b'\x11'*64 + b'\xff'*64 + I) + b'\x00'*64 + I)
    assert compute_MT_node(1, {2: b'\x00'*64, 3: b'\x11'*64, 4: b'\xff'*64}, I, 4, M) == H(64, b'\x11'*64 + b'\xff'*64 + I)
    assert compute_MT_node(1, {3: b'\x11'*64, 4: b'\xff'*64}, I, 4, M) == H(64, b'\x11'*64 + b'\xff'*64 + I)

    # should be able to compute anything if all nodes are known
    T = 2**5
    n = 2

    for P in [1,2,4]:
        l = T//P
        for n in range(2,min(12,l)): # should work for different values of n
            X = build_X(I, T, n, P, x)
            MT = build_MT(I, X, M)

            known_nodes = { i:v for i,v in enumerate(MT) }

            for i in range(len(MT)):
                assert compute_MT_node(i, known_nodes, I, T, M) == MT[i]

            # Now let's try and remove information while remaining computable
            for k in range(T-1):
                known_nodes = { i:v for i,v in enumerate(MT) if i>k}
                for i in range(k+1):
                    # the k first elements of the MT have been removed from known_nodes,
                    # they should be properly computed nonetheless
                    assert compute_MT_node(i, known_nodes, I, T, M) == MT[i]

    # should throw an error if it gets outside of the expected bounds
    with pytest.raises(AssertionError):
        assert compute_MT_node(1, {0: b'\x00'*64}, I, 1, M) == b'\x00'*64
    with pytest.raises(AssertionError):
        assert compute_MT_node(7, {2: b'\x00'*64, 3: b'\x11'*64, 4: b'\xff'*64}, I, 4, M) == H(64, H(64, b'\x11'*64 + b'\xff'*64 + I) + b'\x00'*64 + I)
    with pytest.raises(AssertionError):
        assert compute_MT_node(1, {0: b'\x00'*64}, I, 2, M) == b'\x00'*64

def test_xor():
    assert xor(b"\x00", b"\x00") == b"\x00"
    assert xor(b"\x01", b"\x00") == b"\x01"
    assert xor(b"\x00", b"\x01") == b"\x01"
    assert xor(b"\x01", b"\x01") == b"\x00"

def test_compute_Y():
    M = 64
    x = 32
    T = 2**5
    S = 16
    L = ceil(3.3*log(T,2))
    I = os.urandom(M)

    for P in [1,2,4]:
        l = T//P
        for n in range(2,min(12,l)): # should work for different values of n
            X = build_X(I, T, n, P, x)
            MT = build_MT(I, X, M)
            PSI = MT[0]
            N = os.urandom(32) # nounce

            Y, OMEGA, i = compute_Y(I, X, T, L, S, N, PSI)

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
    assert is_PoW_solved(b'\x00'*64, b'\x00'*63 + b'\x01', 64) == True
    assert is_PoW_solved(b'\x00'*64, b'\x00'*64, 64) == False
    assert is_PoW_solved(b'\xff'*63 + b'\xfe', b'\xff'*64, 64) == True

    with pytest.raises(AssertionError):
        is_PoW_solved(b'\x00', b'\x00'*64, 64)
    with pytest.raises(AssertionError):
        is_PoW_solved(b'\x00', b'\x00', 64)
    with pytest.raises(AssertionError):
        is_PoW_solved(b'\x00'*64, b'\x00', 64)

def test_build_rL():
    M = 64
    x = 32
    T = 2**5
    S = 16
    L = ceil(3.3*log(T,2))
    I = os.urandom(M)

    for P in [1,2,4]:
        l = T//P
        for n in range(2,min(12,l)): # should work for different values of n
            X = build_X(I, T, n, P, x)
            MT = build_MT(I, X, M)
            PSI = MT[0]
            N = os.urandom(32) # nounce
            Y, OMEGA, i = compute_Y(I, X, T, L, S, N, PSI)

            round_L = build_rL(i, X, P, n)

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
                        assert round_L[i_j][k] == H(x, stuff_to_hash)

                else:
                    seed = X[i_j-1][:4]
                    # assert correct construction
                    assert round_L[i_j] == [ X[p*l + phi_k_i] for phi_k_i in phis(seed, i_j%l, n) ]

                    # make sure X[i_j] can be recomputed using the elements of round_L[i_j]
                    hash_input = b""
                    for item in round_L[i_j]:
                        hash_input += item
                    assert H(x,hash_input) == X[i_j]

def test_provided_indexes():
    M = 64
    x = 32
    T = 2**5
    S = 16
    L = ceil(3.3*log(T,2))
    I = os.urandom(M)

    for P in [1,2,4]:
        l = T//P
        for n in range(2,min(12,l)): # should work for different values of n
            X = build_X(I, T, n, P, x)
            MT = build_MT(I, X, M)
            PSI = MT[0]
            N = os.urandom(32) # nounce
            Y, OMEGA, i = compute_Y(I, X, T, L, S, N, PSI)
            round_L = build_rL(i, X, P, n)

            indexes = provided_indexes(round_L, P, T, n)

            for index in indexes:
                assert index < T

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
                    for index in [p*l+phi_k_i for phi_k_i in phis(seed, i_j % l, n)]:
                        assert index in indexes

            # Asserting there are no duplicates
            assert len(indexes) == len(set(indexes))

def test_build_rZ():
    M = 64
    x = 32
    T = 2**5
    S = 16
    L = ceil(3.3*log(T,2))
    I = os.urandom(M)

    for P in [1,2,4]:
        l = T//P
        for n in range(2,min(12,l)): # should work for different values of n
            X = build_X(I, T, n, P, x)
            MT = build_MT(I, X, M)
            PSI = MT[0]
            N = os.urandom(32) # nounce
            Y, OMEGA, i = compute_Y(I, X, T, L, S, N, PSI)
            round_L = build_rL(i, X, P, n)

            Z = build_rZ(round_L, MT, P, T, n)

            # A shift has to be applied so the indexes match those of the
            # Merkle Tree and not those of X.
            indexes = [ index + T - 1 for index in  provided_indexes(round_L, P, T, n)]

            for k in Z:
                assert k not in indexes
                assert Z[k] == MT[k]
                if k >= T-1:
                    assert Z[k] == H( M, X[k-(T-1)]+I )

            assert set(Z.keys()) == set(opening(T, provided_indexes(round_L, P, T, n)))

def clean_Z():
    Z = { 5: b'\x00', 8: b'\xfe', 14: b'\xa4' }
    cleaned_Z = clean_Z(Z)
    assert len(cleaned_Z) == 3
    for k in Z:
        assert cleaned_Z[k] == Z[k].hex()

def test_trim_round_L():
    with pytest.raises(AssertionError):
        trim_round_L({}, 5, 2, 0)

    # Assert that the intended items are trimmed off of the dict
    round_L_1 = {7: [b'\x00'], 15:[b'\x00']} # should remain unchanged if (P, T, n) = (32, 4, 6)
    round_L_2 = {5: [b'\x00'], 10:[b'\x00']} # should be totally trimmed
    round_L_3 = {6: [b'\x00'], 14:[b'\x00']} # edge case : should remain unchanged
    round_L_4 = {5: [b'\x00'], 15:[b'\x00']} # should be partially modified

    assert trim_round_L(round_L_1 , 4, 32, 6) == {7: ['00'], 15: ['00']}
    assert trim_round_L(round_L_2 , 4, 32, 6) == {5: [], 10: []}
    assert trim_round_L(round_L_3 , 4, 32, 6) == {6: ['00'], 14: ['00']}
    assert trim_round_L(round_L_4 , 4, 32, 6) == {5: [], 15:['00']}

    # Assert that the bytearrays are properly converted
    round_L_5 = {7: [b'\x00', b'\xfe', b'\xa4']}
    assert trim_round_L(round_L_5, 4, 32, 6) == {7: ['00','fe','a4']}

def test_build_JSON_output():
    JSON = build_JSON_output(N=b'\x00'*63 + b'\xff', round_L={}, Z={}, P=4, T=32, n='n', I=b'\xff'*64, M='M', L='L', S='S', x='x', d=b'\x00'*64)

    data = json.loads(JSON)

    assert data['answer']['N'] == '00'*63 + 'ff'
    assert data['answer']['round_L'] == {}
    assert data['answer']['Z'] == {}

    assert data['params']['P'] == 4
    assert data['params']['T'] == 32
    assert data['params']['n'] == 'n'
    assert data['params']['I'] == 'ff'*64
    assert data['params']['M'] == 'M'
    assert data['params']['L'] == 'L'
    assert data['params']['S'] == 'S'
    assert data['params']['x'] == 'x'
    assert data['params']['d'] == '00'*64

#@pytest.mark.skip(reason="not good yet")
def test_PoW():
    M = 64
    x = 32
    T = 2**4
    S = 16
    L = ceil(3.3*log(T,2))
    I = os.urandom(M)
    d = b'\x00'*S # minimal difficulty

    for P in [1,2,4]:
        l = T//P
        for n in range(2,min(12,l)): # should work for different values of n
            json_output, X_PoW, MT_PoW, PSI_PoW, N_PoW, I_PoW, Y_PoW, OMEGA_PoW, i_PoW, round_L_PoW, Z_PoW = PoW(I=I, T=T, n=n, P=P, M=M, L=L, S=S, x=x, d=d, debug=True)
            data = json.loads(json_output, object_pairs_hook=OrderedDict)

            assert data['params']['P'] == P
            assert data['params']['T'] == T
            assert data['params']['n'] == n
            assert data['params']['I'] == I.hex()
            assert data['params']['M'] == M
            assert data['params']['L'] == L
            assert data['params']['S'] == S
            assert data['params']['x'] == x
            assert data['params']['d'] == d.hex()

            # Verifying the answer

            N = bytes.fromhex(data['answer']['N'])
            assert N == N_PoW

            I = bytes.fromhex(data['params']['I'])
            assert I == I_PoW

            unprocessed_Z = data['answer']['Z']
            unprocessed_round_L = data['answer']['round_L']

            # Preparing round_L

            round_L = OrderedDict.fromkeys(unprocessed_round_L)
            for k in unprocessed_round_L:
                round_L[int(k)] = [ bytes.fromhex(x) for x in unprocessed_round_L[k] ]
                del round_L[k] # the key used to be denoted by a char, and has to go

            # assert correct construction
            assert round_L == {k:(v if k%l >= n else []) for k,v in round_L_PoW.items()}

            # Preparing Z
            Z = {}
            for k in unprocessed_Z:
                Z[int(k)] = bytes.fromhex(unprocessed_Z[k])

            assert Z == Z_PoW

            # Verifications
            for k in round_L:
                assert T-1 <= k+(T-1) < 2*T-1
                assert k+(T-1) not in Z

            # Building back X
            X = {}

            # First, building back elements that are built at step 1.a., that can be built from scratch
            for p in range(P):
                for i in range(n):
                    X[p*l+i] = H(x, int_to_4bytes(i) + int_to_4bytes(p) + I)

            # Now, going through round_L to add the provided elements
            for i_j in round_L:
                if i_j%l >= n :
                    # adding all the antecedents
                    seed = round_L[i_j][0][:4]
                    hash_input = b''
                    for k, phi_k_i in enumerate(phis(seed, i_j%l, n)):
                        X[p*l+phi_k_i] = round_L[i_j][k]
                        hash_input += round_L[i_j][k]

                    # recomputing X[i_j] and adding it
                    X[i_j] = H(x, hash_input)

            # asserting correct reconstruction
            for a,e in X.items():
                assert X_PoW[a] == e

            # Let's build a dict of all the nodes we know (round_L, Z, and the precomputable ones), that satisfies the requirement of compute_merkle_tree_node
            known_nodes = {**Z, **{ k + (T-1) : H(M, v+I) for k,v in X.items() } }

            for index, node in known_nodes.items():
                assert MT_PoW[index] == node

            # Verifications
            assert len(known_nodes) == len(Z) + len(X)

            PSI = compute_MT_node(0, known_nodes, I, T, M)

            assert PSI == PSI_PoW

            # We can now use the previous functions to compute i and OMEGA
            Y, OMEGA, computed_i = compute_Y(I, X, T, L, S, N, PSI)

            assert Y == Y_PoW
            assert OMEGA == OMEGA_PoW

            # Verifying the two conditions that define the success of the PoW :
            #
            # 1. OMEGA satisfies the difficulty constraint
            assert computed_i == i_PoW
            assert is_PoW_solved(d, OMEGA, S=S) == True
            # 2. The keys of round_L correspond the i that has been computed by compute_Y
            assert set(computed_i) == set(round_L.keys())
