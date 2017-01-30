

def prover(crs, v, sigma, s):
    pass


def inverse_perm(s):
    pass


def step1(n, gk, sigma):
    randoms = [gk.q.random() for i in range(1, n)]  # n-1 values
    inverted_sigma = inverse_perm(sigma)
    expL = lambda i: poly(inverted_sigma(i), Chi.chi)
    exps = [expL(i) + random[i-1] * Chi.rho
            for i in range(1, n)]
    A1 = [gk.g1 * exp for exp in exps]
    A2 = [gk.g2 * exp for exp in exps]
    return randoms, A1, A2


def step2(n, gk, randoms):
    rand_n = - sum(randoms)
    randomds.append(rand_n)
    return randoms


def step3(n, gk, A1, A2, crs):
    inf1 = G1Elem.inf(gk.G)
    inf2 = G2Elem.inf(gk.G)
    prod1 = reduce(lambda x, y: x * y, A1, inf1)
    prod2 = reduce(lambda x, y: x * y, A2, inf2)
    inv_prod1 = prod1.neg()
    inv_prod2 = prod2.neg()
    A1last = crs.last_value * inv_prod1
    A2last = crs.last_value * inv_prod2
    return A1 + [A1last], A2 + [A2last]
