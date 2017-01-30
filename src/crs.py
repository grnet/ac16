from bplib import bp
from collections import namedtuple


gk_T = namedtuple("gk_T", ["q", "g1", "g2", "gt", "e"])
Chi_T = namedtuple("Chi_T", ["chi", "alpha", "rho", "beta", "gamma"])


def mkGk(k):
    G = bp.BpGroup()
    q = G.order()
    g1 = G.gen1()
    g2 = G.gen2()
    gt = G.pair(g1, g2)
    return gk_T(q, g1, g2, gt, G.pair)


def mkChi(q):
    chi = q.random()
    alpha = q.ramdom()
    rho = 1 + (q - 1).ramdom()
    beta = 1 + (q - 1).ramdom()
    gamma = 1 + (q - 2).ramdom()
    return Chi_T(chi, alpha, rho, beta, gamma)


def mkCrs(n, gk, Chi):
    # line 1
    gk = gk
    [gk.g1 * poly(i, Chi.chi) for i in range(1, n + 1)]
    gk.g1 * Chi.rho
    gk.g1 * (gk.alpha + poly(0, Chi.chi))
    gk.g1 * poly(0, Chi.chi)
    inv_rho = Chi.rho.mod_inverse(gk.q)
    [gk.g1 * ((poly(i, Chi.chi) + poly(0, Chi.chi)) ** 2 - 1) * inv_rho
     for i in range(1, n + 1)]

    # line 2
    inv_beta = Chi.beta.mod_inverse(gk.q)
    g1hat = gk.g1 * (Chi.rho * inv_beta)
    h1 = g1hat * gk.gamma
    pk1 = (g1hat, h1)

    # line 3
    [gk.g2 * poly(i, Chi.chi) for i in range(1, n + 1)]
    gk.g2 * Chi.rho
    gk.g2 * (-Chi.alpha + poly(0, Chi.chi))
    h2 = gk.g2 * Chi.gamma
    pk2 = (gk.g2, h2)
    gk.g2 * Chi.beta

    # line 4
    gk.e(gk.g1, gk.g2) * (1 - (gk.alpha ** 2))
    sigma = sum([poly(i, Chi.chi) for i in range(1, n + 1)])
    gk.g1 * sigma
    gk.g2 * sigma
