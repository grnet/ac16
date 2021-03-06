from petlib.ec import EcGroup
from bplib.bp import BpGroup

def gen_params_ec_group(nid=713):
    G = EcGroup(nid)
    g = G.generator()
    o = G.order()
    return (G, g, o)

def gen_params_bp_g1(nid=713):
    G = BpGroup()
    g = G.gen1()
    o = G.order()
    return (G, g, o)

def gen_params_bp_g2(nid=713):
    G = BpGroup()
    g = G.gen2()
    o = G.order()
    return (G, g, o)

def key_gen(params):
    _, g, o = params
    x = 1 + (o - 2).random()
    h = x * g
    return ((g, h), x)
    
def enc(pk, s1, s2, m):
    g, h = pk
    return (s1*h, s2*(g+h), (m + s1 + s2) * g)

def dec(c, o, sk, table):
    c1, c2, c3 = c
    e1 = -(sk).mod_inverse(o)    
    e2 = -(sk + 1).mod_inverse(o)
    v = (c3 + e2*c2 + e1*c1)
    return table[v]

def make_table(g, n):
    table = {}
    for i in range(n):
        elem = (i * g)
        table[elem] = i
    return table

def test_encdec(params):
    table = make_table(params[1], 1000)
    G, g, o = params
    pk, sk = key_gen(params)
    s1 = o.random()
    s2 = o.random()
    c = enc(pk, s1, s2, 666)
    assert(dec(c, o, sk, table) == 666)
    c = enc(pk, s1, s2, 7)
    assert(dec(c, o, sk, table) == 7)
    import random
    ps = random.sample(range(1000), 100)
    for i in range(100):
        c = enc(pk, s1, s2, ps[i])
        assert(dec(c, o, sk, table) == ps[i])

if __name__ == '__main__':
    print("Testing EC group...")
    test_encdec(gen_params_ec_group())
    print("Testing BP group G1...")
    test_encdec(gen_params_bp_g1())
    print("Testing BP group G2...")
    test_encdec(gen_params_bp_g2())
