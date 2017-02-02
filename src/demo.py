import random
import datetime
import sys

import crs
import encdec
import prover
import verifier


system_random = random.SystemRandom()


def secure_shuffle(lst):
    random.shuffle(lst, random=system_random.random)


def random_permutation(n):
    s = list(range(n))
    secure_shuffle(s)
    return s


def initialize(n):
    gk = crs.mk_gk()
    Chi = crs.mk_Chi(gk.q)
    CRS, td = crs.mk_crs(n, gk, Chi)
    return gk, Chi, CRS


def mk_s_randoms(n, q):
    return [[q.random() for j in range(2)] for i in range(n)]


def demo(n, messages):
    gk, Chi, CRS = initialize(n)
    secret = Chi.gamma
    pk1 = CRS.pk1
    pk2 = CRS.pk2
    ciphertexts = encrypt_messages(gk.q, pk1, pk2, messages)

    start = datetime.datetime.now()
    sigma = random_permutation(n)
    print("SIGMA = %s" % sigma)
    s_randoms = mk_s_randoms(n, gk.q)
    proof = prover.prove(n, CRS, ciphertexts, sigma, s_randoms)
    shuffled_ciphertexts, \
        (A1, A2), pi_1sp, pi_c1_1, pi_c1_2, pi_c2_1, pi_c2_2 = proof

    perm_ok, valid, consistent = verifier.verify(
        n, CRS, ciphertexts, shuffled_ciphertexts,
        A1, A2, pi_1sp, pi_c1_1, pi_c1_2, pi_c2_1, pi_c2_2)
    print("VERIFY: %s %s %s" % (perm_ok, valid, consistent))
    end = datetime.datetime.now()

    TABLES = encdec.make_tables(pk1, pk2, n)
    shuffled_ms = decrypt_messages(gk.q, secret, TABLES, shuffled_ciphertexts)
    print(shuffled_ms)
    print("elapsed: %s" % (end - start))


def encrypt_messages(order, pk1, pk2, messages):
    return [encdec.encrypt(order, pk1, pk2, message) for message in messages]


def decrypt_messages(order, secret, tables, ciphertexts):
    return [encdec.decrypt(order, cs, secret, tables) for cs in ciphertexts]


if __name__ == '__main__':
    n = 10
    if len(sys.argv) >= 2:
        n = int(sys.argv[1])
    messages = list(range(n))
    demo(len(messages), messages)
