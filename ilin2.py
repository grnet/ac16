from petlib.ec import EcGroup

def params_gen(nid=713):
    G = EcGroup(nid)
    g = G.generator()
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

def dec(c, g, sk, table):
    _, g, o = params
    c1, c2, c3 = c
    e1 = -(sk).mod_inverse(o)    
    e2 = -(sk + 1).mod_inverse(o)
    v = (c3 + e2*c2 + e1*c1)
    return table[v]

def make_table(params):
    _, g, o = params
    table = {}
    for i in range(-1000, 1000):
        table[i * g] = i
    return table

params = params_gen()
table = make_table(params)
G, g, o = params
pk, sk = key_gen(params)

s1 = o.random()
s2 = o.random()
c = enc(pk, s1, s2, 666)
print(dec(c, g, sk, table))
