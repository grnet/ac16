import prover
from bplib.bp import GTElem


def step1(gk, A1, A2, g1_sum, g2_sum):
    return prover.step1c(gk, A1, A2, g1_sum, g2_sum)


def step2(gk, n):
    p1 = [gk.q.random() for i in range(n)]
    p2 = [gk.q.random() for j in range(3)]
    p3 = [[gk.q.random() for j in range(3)]
          for i in range(n)]
    p4 = [gk.q.random() for j in range(3)]
    return p1, p2, p3, p4


def get_infT(gk):
    return GTElem.one(gk.G)


def step3(gk, e, A1, A2, p1, pi_1sp, g1alpha, g2alpha, g2rho, pair_alpha):
    inf1, inf2 = prover.get_infs(gk)
    infT = get_infT(gk)
    prodT = infT
    prod1 = inf1
    sum_p = 0
    assert len(A1) == len(A2) == len(p1) == len(pi_1sp)
    for (Ai1, Ai2, p1i, pi_1spi) in zip(A1, A2, p1, pi_1sp):
        prodT *= e(p1i * (Ai1 + g1alpha), Ai2 + g2alpha)
        prod1 += p1i * pi_1spi
        sum_p = (sum_p + p1i) % gk.q
    right1 = e(prod1, g2rho)
    right2 = pair_alpha ** sum_p
    right = right1 * right2
    print 1, right1 == infT
    print 2, right2 == infT
    print 2, type(right2)
    print "INFS"
    print sum_p
    print prod1 == inf1
    return prodT == right


def step4(gk, e, p2, p3, pi_c2_1, pi_c2_2, v_primes, g1rho, g2beta):
    inf1, inf2 = prover.get_infs(gk)

    def pi_c_prod(inf, pi_c2_):
        prod_c2_ = inf
        for (p2j, pi_c2_j) in zip(p2, pi_c2_):
            prod_c2_ += p2j * pi_c2_j
        return prod_c2_

    def nested_prods(inf, flag):
        outer_prod = inf
        for vi_prime, p3i in zip(v_primes, p3):
            inner_prod = inf
            vi_f_prime = vi_prime[flag]
            for (vi_f_prime_j, p3ij) in zip(vi_f_prime, p3i):
                inner_prod += p3ij * vi_f_prime_j
            outer_prod += inner_prod
        return outer_prod

    left = e(g1rho, pi_c_prod(inf2, pi_c2_2) + nested_prods(inf2, 1))
    right = e(pi_c_prod(inf1, pi_c2_1) + nested_prods(inf1, 0), g2beta)
    return left == right


def step5(gk, pk1, g2rho, pi_c1_1, pi_c1_2, pi_c2_1, e, p4):
    inf1, _ = prover.get_infs(gk)
    g1hat, h1 = pk1
    pair1 = e(g1hat, p4[1] * pi_c1_2 + p4[2] * (pi_c1_1 + pi_c1_2))
    pair2 = e(h1, p4[0] * pi_c1_1 + p4[1] * pi_c1_2)
    prod = inf1
    for (p4j, pi_c2_1j) in zip(p4, pi_c2_1):
        prod += p4j * pi_c2_1j
    pair3 = e(prod, g2rho)
    return pair1 + pair2 - pair3


def step6(gk, e, vs, v_primes, A2, p4, R, g2_polys):
    def do_inner(vi):
        vi1 = vi[0]
        inf1, _ = prover.get_infs(gk)
        inner_prod = inf1
        for (p4j, vi1j) in zip(p4, vi1):
            inner_prod += p4j * vi1j
        return inner_prod

    infT = get_infT(gk)
    outer_numer = infT
    for (vi_prime, g2_poly_i) in zip(v_primes, g2_polys):
        outer_numer += e(do_inner(vi_prime), g2_poly_i)

    outer_denom = infT
    for (vi, Ai2) in zip(vs, A2):
        outer_denom += e(do_inner(vi), Ai2)

    return outer_numer - outer_denom == R


def verify(n, crs, vs, v_primes, A1, A2,
           pi_1sp, pi_c1_1, pi_c1_2, pi_c2_1, pi_c2_2):
    gk = crs.gk
    A1, A2 = step1(gk, A1, A2, crs.g1_sum, crs.g2_sum)
    p1, p2, p3, p4 = step2(gk, n)
    perm_ok = step3(gk, gk.e, A1, A2, p1, pi_1sp,
                    crs.g1alpha, crs.g2alpha, crs.g2rho, crs.pair_alpha)
    valid = step4(
        gk, gk.e, p2, p3, pi_c2_1, pi_c2_2, v_primes, crs.g1rho, crs.g2beta)
    R = step5(gk, crs.pk1, crs.g2rho, pi_c1_1, pi_c1_2, pi_c2_1, gk.e, p4)
    consistent = step6(gk, gk.e, vs, v_primes, A2, p4, R, crs.g2_polys)
    return perm_ok, valid, consistent
