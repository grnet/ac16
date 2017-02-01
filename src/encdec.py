import ilin2


def encrypt(q, pk1, pk2, m):
    s1 = q.random()
    s2 = q.random()
    c1 = ilin2.enc(pk1, s1, s2, m)
    c2 = ilin2.enc(pk2, s1, s2, m)
    return c1, c2


def make_tables(pk1, pk2, n):
    table1 = ilin2.make_table(pk1[0], n)
    table2 = ilin2.make_table(pk2[0], n)
    return table1, table2


def decrypt(gk, ctext, secret, tables):
    m1 = ilin2.dec(ctext[0], gk.q, secret, tables[0])
    m2 = ilin2.dec(ctext[1], gk.q, secret, tables[1])
    return m1, m2


def test(gk, Chi, CRS):
    pk1 = CRS.pk1
    pk2 = CRS.pk2
    tables = make_tables(pk1, pk2, 100)
    secret = Chi.gamma
    m = 17
    ctexts = encrypt(gk.q, pk1, pk2, m)
    m1, m2 = decrypt(gk, ctexts, secret, tables)
    return m1, m2
