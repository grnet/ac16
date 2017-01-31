from petlib.bn import Bn

def lagrangian(i, n, chi, q):
    """Evaluates Lagrange basis polynomial l_i at point x.
    
    Distinct points are 1, 2, ..., n + 1.
    Returns nonzero field element.
    
    Arguments:
    i -- polynomial index in range 1, 2, ..., n + 1
    chi -- input value of the polynomial
    q -- the order of the group
    """
    numerator = Bn(1)
    denominator = Bn(1)
    for j in range(1, n+2):
        if i == j:
            continue

        numerator = numerator.mod_mul(chi - j, q)
        elem = i - j
        denominator = denominator.mod_mul(elem, q)

    return numerator.mod_mul(denominator.mod_inverse(q), q)

def compute_denominators(k, q):
    """Computes denominators for Lagrange basis polynomials.
    
    Uses distinct points 1, ...,k
    Arguments:
    k -- number of basis polynomials
    q -- the order of the group
    """
    denominators = []
    temp = Bn(1)
    for i in range(1, k+1):
        if i == 1:
            for j in range(2, k+1):
                elem = i - j;
                temp = temp.mod_mul(elem, q)
        elif i == k:
            elem = 1 - k;
            temp = temp.mod_mul(elem, q)
        else:            
            inverse = Bn(i - 1 - k)
            inverse = inverse.mod_inverse(q)
            elem = i - 1
            temp = temp.mod_mul(elem, q)
            temp = temp.mod_mul(inverse, q)
        denominators.append(temp)
    return denominators


def compute_denominators_slow(k, q):
    """Computes denominators for Lagrange basis polynomials.
    This function can be used to test compute_denominators(k, q).

    Uses distinct points 1, ...,k
    Arguments:
    k -- number of basis polynomials
    q -- the order of the group
    """
    denominators = []

    for i in range(1, k+1):
        denominator = Bn(1)
        for j in range(1, k+1):
            if i == j:
                continue
            elem = i - j
            denominator = denominator.mod_mul(elem, q)
        denominators.append(denominator)
    return denominators


def generate_pi(chi, n, q):
    """Computes vector of elements P_i(chi) for i = 1, ..., n.

    Uses Lagrange basis polynomials with distinct points 1, ..., n + 1
    Arguments:
    chi -- point of evaluation
    n -- the number of elements
    q -- the order of the group
    """
    
    if chi <= n + 1:
        raise ValueError("chi should be greater than n + 1, chi=%s n+1=%s"
                         % (chi, n+1))
    pis = []
            
    prod = Bn(1)
    # prod = (x - w_1) (x - w_2) ... (x - w_{n+1})
    for j in range(1, n + 2):
        prod = prod.mod_mul(chi - j, q)

    # denoms[0] = 1 / (w_1 - w_2) (w_1 - w_3) ... (w_1 - w_{n + 1})
    # denoms[1] = 1 / (w_2 - w_1) (w_2 - w_3) ... (w_2 - w_{n + 1})
    # denoms[n+1] = 1 / (w_{n+1}- w_1) (w_{n+1} - w_2) ... (w_{n+1} - w_n)
    denoms = compute_denominators(n + 1, q)
            
    missing_factor = Bn(chi - (n + 1))
            
    ln_1 = prod.mod_mul(missing_factor.mod_inverse(q), q)
    ln_1 = ln_1.mod_mul(denoms[n].mod_inverse(q), q)

    two = Bn(2)
    for i in range(1, n + 1):
        missing_factor = Bn(chi - i)
        l_i = prod.mod_mul(missing_factor.mod_inverse(q), q)
        l_i = l_i.mod_mul(denoms[i - 1].mod_inverse(q), q)
        pis.append(two.mod_mul(l_i, q).mod_add(ln_1, q))
        
    return pis

