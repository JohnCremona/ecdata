######################################################################
#
# Functions for Minkowski-reduction of generators, and naive reduction
# of torsion generators and of generators mod torsion.
#
######################################################################

pt_wt = lambda P: len(str(P))

def reduce_tgens(tgens, verbose=False):
    """
    tgens: list of torsion generators (if two, sorted by order)

    Return a new list of generators which is minimal with respect to string length.
    """
    r = len(tgens)
    if r == 0:
        return tgens
    if r == 1: # cyclic
        P1 = tgens[0]
        n1 = P1.order()
        if n1 == 2: # no choice
            return tgens
        Plist = [i * P1 for i in range(1, n1) if n1.gcd(i) == 1]
        Plist.sort(key=pt_wt)
        Q = Plist[0]
        if Q != P1 and verbose:
            print("Replacing torsion [{}] generator {} with {}".format(n1, P1, Q))
        return [Q]
    # now r=2 and P1 has order n1=2  while n2 = 2, 4, 6, 8.
    # -- we use brute force
    assert r == 2
    if tgens[0].order() > tgens[1].order():
        tgens.reverse()
    P1, P2 = tgens
    n1 = P1.order() # = 2
    n2 = P2.order() # = 2, 4, 6 or 8
    assert n1 == 2 and n2 in [2, 4, 6, 8]
    m = n2 // 2
    P1a = m * P2    # other 2-torsion
    P1b = P1 + P1a  # points
    gen_pairs = []
    for j in range(n2):
        jP2 = j * P2
        for i in range(n1):
            Q = i * P1 + jP2
            if Q.order() != n2:
                continue
            Q2 = m * Q
            for P in [P1, P1a, P1b]:
                if P != Q2:
                    gen_pairs.append((P, Q))
    pt_wt2 = lambda PQ: sum(pt_wt(P) for P in PQ)
    gen_pairs.sort(key=pt_wt2)
    rtgens = list(gen_pairs[0])
    # for structure [2,2] we swap the two gens over so that the first has smallest x-coordinate
    if n2 == 2:
        rtgens.sort(key=lambda P: list(P)[0])
    if rtgens != tgens and verbose:
        print("Replacing torsion [{},{}] generators {} with {}".format(n1, n2, tgens, rtgens))
    return rtgens

def check_minkowski(gens):
    """
    Check the points are Minkowski-reduced (rank up to 3 only)
    """
    r = len(gens)
    if r < 2 or r > 3:
        return True
    if r == 2:
        P1, P2 = gens
        h1 = P1.height()
        h2 = P2.height()
        h3 = (P1 + P2).height()
        return h1 < h2 and h2 < h3 and h3 < 2 * h1 + h2
    # r=3
    if not check_minkowski(gens[:2]):
        return False
    if not check_minkowski(gens[1:]):
        return False
    if not check_minkowski(gens[::2]):
        return False
    P1, P2, P3 = gens
    h3 = P3.height()
    P4 = P1 + P2
    if h3 > (P3 + P4).height() or h3 > (P3 - P4).height():
        return False
    P4 = P1 - P2
    return h3 < (P3 + P4).height() and h3 < (P3 - P4).height()

def reduce_mod_2d(P3, P1, P2, debug=False):
    """
    Assuming [P1,P2] reduced, return P3+n1*P1+n2*P2 of minimal height
    """
    if debug:
        print("Reducing {} mod [{},{}]".format(P3, P1, P2))
    h1 = P1.height()
    h2 = P2.height()
    assert h1 <= h2
    P12 = P1 + P2
    h12 = (P12.height() - h1 - h2) / 2
    assert 2 * h12.abs() <= h1

    # now the height of x*P1+y*P2 is ax^2+2bxy+cy^2

    h3 = P3.height()
    h13 = ((P1 + P3).height() - h1 - h3) / 2
    h23 = ((P2 + P3).height() - h2 - h3) / 2

    d = h1 * h2 - h12 * h12
    y1 = (h2 * h13 - h12 * h23) / d
    y2 = (h1 * h23 - h12 * h13) / d

    # now y1*P1+y2*p2 is the orthogonal projection of P3 onto the P1,P2-plane

    n1 = y1.round()
    n2 = y2.round()

    if debug:
        print("orthog proj has coords ({},{}), rounded to ({},{})".format(y1, y2, n1, n2))
    Q3 = P3 - (n1 * P1 + n2 * P2) # approximate answer
    if debug:
        print("base reduction is {}, height {}".format(Q3, Q3.height()))
    P21 = P1 - P2
    Q3list = [Q3, Q3 - P1, Q3 + P1, Q3 - P2, Q3 + P2, Q3 - P12, Q3 + P12, Q3 - P21, Q3 + P21]
    Q3list.sort(key=lambda P: P.height())
    if debug:
        print("candidates for reduction: {}".format(Q3list))
        print(" with heights: {}".format([Q.height() for Q in Q3list]))
        R = P3 - P1 - P2
        print(" P3-P1-P2 ={} has height {}".format(R, R.height()))
    return Q3list[0]

def mreduce_gens(gens, debug=False):
    r = len(gens)
    if r < 2 or r > 3:
        return gens
    if r == 2:
        P1, P2 = gens
        h1 = P1.height()
        h2 = P2.height()
        h12 = ((P1 + P2).height() - h1 - h2) / 2

        while True:
            x = (h12/h1).round()
            y = h12 - x * h1
            P1, P2 = P2 - x * P1, P1
            h1, h2 = h2, h1
            h1 = h1 - x * (y + h12)
            h12 = y
            if h1 > h2:
                return [P2, P1]

    # now r=3
    if check_minkowski(gens):
        return gens
    P1, P2, P3 = gens
    if debug:
        print("--------------------------------------------------------")
    while True:
        if debug:
            print("At top of loop: {}".format([P1, P2, P3]))
        P1, P2 = mreduce_gens([P1, P2]) # recursive
        if debug:
            print("After one 2D step: {}".format([P1, P2, P3]))
        P3 = reduce_mod_2d(P3, P1, P2)
        if debug:
            print("After reducing P3 mod [P1,P2]: {}".format([P1, P2, P3]))
        h3 = P3.height()
        if h3 >= P2.height():
            newgens = [P1, P2, P3]
            if not check_minkowski(newgens):
                label = gens[0].curve().label()
                print("{}: gens = {}, newgens = {} are not Minkowski-reduced!".format(label, gens, newgens))
                print("heights are {}".format([P.height() for P in newgens]))
                raise RuntimeError
            return newgens
        P1, P2, P3 = (P1, P3, P2) if P1.height() < h3 else (P3, P1, P2)
        if debug:
            print("After resorting: {}".format([P1, P2, P3]))

def reduce_gens(gens, tgens, verbose=False, label=None):
    """
    gens: list of generators mod torsion
    tgens: list of torsion generators

    (1) Reduce the torsion generators w.r.t. string length
    (2) Minkowki reduce the mod-torsion gens (or just LLL-reduce if rank>3)
    (3) Reduce the mod-torsion gens mod torsion w.r.t. string length

    """
    rtgens = reduce_tgens(tgens)
    if not gens:
        return [], rtgens
    E = gens[0].curve()
    newgens = mreduce_gens(E.lll_reduce(gens)[0]) # discard transformation matrix

    Tlist = E.torsion_points() if tgens else [E(0)]

    def reduce_one(P):
        mP = -P
        Plist = [P + T for T in Tlist] + [mP + T for T in Tlist]
        Plist.sort(key=pt_wt)
        return Plist[0]

    newgens = [reduce_one(P) for P in newgens]

    if verbose and len(gens) > 1 and newgens != gens:
        print("replacing {} generators ({}) {} with {}".format(len(gens), label, gens, newgens))

    return newgens, rtgens
