from bplib import bp
from collections import namedtuple
from pisgen import generate_pi


gk_T = namedtuple("gk_T", ["q", "g1", "g2", "gt", "e"])
Chi_T = namedtuple("Chi_T", ["chi", "alpha", "rho", "beta", "gamma"])
CRS_T = namedtuple(
    "CRS_T", ["gk", "g1_polys", "g1rho", "g1alpha", "g1_poly_zero",
              "g1_poly_squares", "pk1", "g2_polys", "g2rho", "g2alpha",
              "pk2", "g2beta", "pairing", "g1_sum", "g2_sum"])


def mk_gk(k):
    G = bp.BpGroup()
    q = G.order()
    g1 = G.gen1()
    g2 = G.gen2()
    gt = G.pair(g1, g2)
    return gk_T(q, g1, g2, gt, G.pair)


def mk_Chi(q):
    chi = q.random()
    alpha = q.random()
    rho = 1 + (q - 1).random()
    beta = 1 + (q - 1).random()
    gamma = 1 + (q - 2).random()
    return Chi_T(chi, alpha, rho, beta, gamma)


def mk_crs(n, gk, Chi):
    polys_all = generate_pi(Chi.chi, n + 1, gk.q)
    poly_zero = polys_all[0]
    polys = polys_all[1:]

    # line 1
    gk = gk
    g1_polys = [poly * gk.g1 for poly in polys]
    g1rho = Chi.rho * gk.g1
    g1alpha = (gk.alpha + poly_zero) * gk.g1
    g1_poly_zero = poly_zero * gk.g1
    inv_rho = Chi.rho.mod_inverse(gk.q)
    g1_poly_squares = []
    for poly in polys:
        nom = (poly + poly_zero) ** 2 - 1
        g1_poly_squares.append((nom * inv_rho) * gk.g1)

    # line 2
    inv_beta = Chi.beta.mod_inverse(gk.q)
    g1hat = (Chi.rho * inv_beta) * gk.g1
    h1 = gk.gamma * g1hat
    pk1 = (g1hat, h1)

    # line 3
    g2_polys = [poly * gk.g2 for poly in polys]
    g2rho = Chi.rho * gk.g2
    g2alpha = (-Chi.alpha + poly_zero) * gk.g2
    h2 = Chi.gamma * gk.g2
    pk2 = (gk.g2, h2)
    g2beta = Chi.beta * gk.g2

    # line 4
    pairing = (1 - gk.alpha ** 2) * gk.e(gk.g1, gk.g2)
    poly_sum = sum([poly for poly in polys])
    g1_sum = poly_sum * gk.g1
    g2_sum = poly_sum * gk.g2

    CRS = CRS_T(gk, g1_polys, g1rho, g1alpha, g1_poly_zero,
                g1_poly_squares, pk1, g2_polys, g2rho, g2alpha,
                pk2, g2beta, pairing, g1_sum, g2_sum)
    trapdoor = (Chi.chi, Chi.rho)
    return CRS, trapdoor
