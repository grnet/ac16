import prover
from bplib.bp import GTElem


def step1(gk, A1, A2, g1_sums, g2_sums):
    return prover.step1c(gk, A1, A2, g1_sums, g2_sums)


def step2(gk, n):
    p1 = [gk.q.random() for i in range(n)]
    p2 = [gk.q.random() for j in range(3)]
    p3 = [[gk.q.random() for j in range(3)]
          for i in range(n)]
    p4 = [gk.q.random() for j in range(3)]
    return p1, p2, p3, p4


def get_infT(gk):
    return GTElem.zero(gk.G)


def step3(gk, e, A1, A2, p1, pi_1sp, g1alpha, g2alpha, g2rho, pair_alpha):
    inf1, inf2 = prover.get_infs(gk)
    infT = get_infT(gk)
    prodT = infT
    prod1 = inf1
    sum_p = 0
    for (Ai1, Ai2, p1i, pi_1spi) in zip(A1, A2, p1, pi_1sp):
        prodT += e(p1i * (Ai1 + g1alpha), Ai2 + g2alpha)
        prod1 += p1i * pi_1spi
        sum_p += p1i
    right = e(prod1, g2rho) + sum_p * pair_alpha
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
    right = e(pi_c_prod(pi_c2_1) + nested_prods(inf1, 0), g2beta)
    return left == right


def verify(
        crs, vs, v_primes, A1, A2, pi_1sp, pi_c1_1, pi_c1_2, pi_c2_1, pi_c2_2):
    pass
