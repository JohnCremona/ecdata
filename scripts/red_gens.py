def reduce_tgens(tgens, verbose=False):
    """
    tgens: list of torsion generators
    """
    if not tgens:
        return tgens
    P1 = tgens[0]
    n1 = P1.order()
    r = len(tgens)
    pt_wt = lambda P: len(str(P))
    if r==1: # cyclic
        if n1==2: # no choice
            return tgens
        Plist = [i*P1 for i in sxrange(1,n1) if n1.gcd(i)==1]
        Plist.sort(key=pt_wt)
        Q = Plist[0]
        if Q!=P1 and verbose:
            print("Replacing torsion [{}] generator {} with {}".format(n1,P1,Q))
        return [Q]
    # now r=2 and P1 has order n1=2
    assert r==2
    assert n1==2
    P2 = tgens[1]
    n2 = P2.order() # = 2, 4, 6 or 8
    m = n2//2
    P1a = m*P2    # other 2-torsion
    P1b = P1+P1a  # points
    gen_pairs = []
    for j in range(n2):
        jP2 = j*P2
        for i in range(n1):
            Q = i*P1+jP2
            if Q.order() != n2:
                continue
            Q2 = m*Q
            for P in [P1,P1a,P1b]:
                if not P==Q2:
                    gen_pairs.append((P,Q))
    pt_wt2 = lambda PQ: sum(pt_wt(P) for P in PQ)
    gen_pairs.sort(key = pt_wt2)
    rtgens = list(gen_pairs[0])
    if rtgens != tgens and verbose:
        print("Replacing torsion [{},{}] generators {} with {}".format(n1,n2,tgens,rtgens))
    return rtgens

def check_minkowski(gens):
    r = len(gens)
    if r<2 or r>3:
        return True
    if r==2:
        P1, P2 = gens
        h1 = P1.height()
        h2 = P2.height()
        h3 = (P1+P2).height()
        return h1<h2 and h2<h3 and h3<2*h1+h2
    # r=3
    if not check_minkowski(gens[:2]):
        return False
    if not check_minkowski(gens[1:]):
        return False
    if not check_minkowski(gens[::2]):
        return False
    P1, P2, P3 = gens
    h3 = P3.height()
    P4 = P1+P2
    if h3 > (P3+P4).height() or h3 > (P3-P4).height():
        return False
    P4 = P1-P2
    return h3 < (P3+P4).height() and h3 < (P3-P4).height()

def mreduce_gens(gens):
    r = len(gens)
    if r!=2:
        return gens
    if r==2:
        P1, P2 = gens
        h1 = P1.height()
        h2 = P2.height()
        h12 = ((P1+P2).height()-h1-h2)/2

        while True:
            x = (h12/h1).round()
            y = h12 - x*h1
            P1, P2 = P2-x*P1, P1
            h1, h2 = h2, h1
            h1  = h1 - x*(y+h12)
            h12 = y
            if h1>h2:
                return [P2,P1]

def reduce_gens(gens, tgens, verbose=False, label=None):
    """
    gens: list of generators mod torsion
    tgens: list of torsion generators
    """
    rtgens = reduce_tgens(tgens)
    if not gens:
        return [], rtgens
    E = gens[0].curve()
    newgens = E.lll_reduce(gens)[0] # discard transformation matrix

    Tlist = E.torsion_points() if tgens else [E(0)]

    for i,P in enumerate(newgens):
        mP = -P
        Plist = [P+T for T in Tlist] + [mP+T for T in Tlist]
        Plist.sort(key = lambda Q: len(str(Q)))
        newgens[i] = Plist[0]
    if verbose and len(gens)>1 and newgens!=gens:
            print("replacing {} generators ({}) {} with {}".format(len(gens), label, gens, newgens))
    return newgens, rtgens

