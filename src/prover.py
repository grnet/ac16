from bplib.bp import G1Elem, G2Elem
from ilin2 import enc


def inverse_perm(s):
    r = [None] * len(s)
    for index, value in enumerate(s):
        r[value] = index
    return r


def step1a(n, gk, sigma, g1_polys, g1rho, g2_polys, g2rho):
    randoms = [gk.q.random() for i in range(n - 1)]  # n-1 values
    inverted_sigma = inverse_perm(sigma)
    A1 = []
    A2 = []
    for inv_i, ri in zip(inverted_sigma, randoms):
        p1_value = g1_polys[inv_i]
        p2_value = g2_polys[inv_i]
        A1.append(p1_value + ri * g1rho)
        A2.append(p2_value + ri * g2rho)
    return randoms, A1, A2


def step1b(gk, randoms):
    rand_n = - sum(randoms) % gk.q
    randoms.append(rand_n)
    return randoms


def get_infs(gk):
    inf1 = G1Elem.inf(gk.G)
    inf2 = G2Elem.inf(gk.G)
    return inf1, inf2


def step1c(gk, A1, A2, g1_sums, g2_sums):
    inf1, inf2 = get_infs(gk)
    prod1 = sum(A1, inf1)
    prod2 = sum(A2, inf2)
    A1.append(g1_sums - prod1)
    A2.append(g2_sums - prod2)
    return A1, A2


def step2a(sigma, A1, randoms, g1_poly_zero, g1rho, g1_poly_squares):
    pi_1sp = []
    inverted_sigma = inverse_perm(sigma)
    for inv_i, ri, Ai1 in zip(inverted_sigma, randoms, A1):
        g1i_poly_sq = g1_poly_squares[inv_i]
        v = (2 * ri) * (Ai1 + g1_poly_zero) - (ri * ri) * g1rho + g1i_poly_sq
        pi_1sp.append(v)
    return pi_1sp


def step3a(sigma, ciphertexts, s_randoms, pk1, pk2):
    v1s_prime = []
    v2s_prime = []
    for perm_i, s_random in zip(sigma, s_randoms):
        (v1, v2) = ciphertexts[perm_i]
        v1s_prime.append(tuple_add(v1, enc(pk1, s_random[0], s_random[1], 0)))
        v2s_prime.append(tuple_add(v2, enc(pk2, s_random[0], s_random[1], 0)))
    return list(zip(v1s_prime, v2s_prime))


def step4a(gk, s_randoms, g2_polys, g2rho):
    rs = tuple([gk.q.random() for i in range(2)])
    (rs1, rs2) = rs
    pi_c1_1 = rs1 * g2rho
    pi_c1_2 = rs2 * g2rho
    for si, g2_polyi in zip(s_randoms, g2_polys):
        si1, si2 = si
        pi_c1_1 += si1 * g2_polyi
        pi_c1_2 += si2 * g2_polyi
    return rs, (pi_c1_1, pi_c1_2)


def tuple_map(func, tpl):
    return tuple(map(func, tpl))


def tuple_add(tpl1, tpl2):
    zipped = zip(tpl1, tpl2)
    return tuple(z[0] + z[1] for z in zipped)


def step4b(gk, ciphertexts, rs, randoms, pk1, pk2):
    pi_c2_1 = enc(pk1, rs[0], rs[1], 0)
    pi_c2_2 = enc(pk2, rs[0], rs[1], 0)
    for ciphertext, ri in zip(ciphertexts, randoms):
        v1, v2 = ciphertext
        pi_c2_1 = tuple_add(pi_c2_1, tuple_map(lambda x: ri * x, v1))
        pi_c2_2 = tuple_add(pi_c2_2, tuple_map(lambda x: ri * x, v2))
    return pi_c2_1, pi_c2_2


def prove(n, crs, ciphertexts, sigma, s_randoms):
    randoms, A1, A2 = step1a(
        n, crs.gk, sigma, crs.g1_polys, crs.g1rho, crs.g2_polys, crs.g2rho)
    randoms = step1b(crs.gk, randoms)
    A1, A2 = step1c(crs.gk, A1, A2, crs.g1_sum, crs.g2_sum)
    pi_1sp = step2a(
        sigma, A1, randoms, crs.g1_poly_zero, crs.g1rho, crs.g1_poly_squares)
    v_prime = step3a(sigma, ciphertexts, s_randoms, crs.pk1, crs.pk2)
    rs, (pi_c1_1, pi_c1_2) = step4a(crs.gk, s_randoms, crs.g2_polys, crs.g2rho)
    pi_c2_1, pi_c2_2 = step4b(
        crs.gk, ciphertexts, rs, randoms, crs.pk1, crs.pk2)
    return (v_prime, (A1[:-1], A2[:-1]),
            pi_1sp, pi_c1_1, pi_c1_2, pi_c2_1, pi_c2_2)
