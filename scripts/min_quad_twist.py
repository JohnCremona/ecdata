# May 2023 new function for cmputing minimal quadratic twist for any curve /Q
#
# - for j not 0 or 1728 this only depends on the j-invariant
# - for curves with CM (and j not 0, 1728) we use a lookup table for speed
# - for non-CM curves it's enough to
#   (1) find minimal conductor;
#   (2) sort into isogeny classes (if more than one curve);
#   (3) return the unique curve in the "first" class if not CM, or
#   (4) for j=0 or 1728, at step (3) there will be 2 curves in the first
#   class with distinct discriminants of the same sign, and we
#   return the one with smaller absolute discriminant.
#
# NB In step (4) for j=287496, CM discriminant -16 only, both curves
# in the first isogeny class have the same discriminant, and we return
# the one with positive c6, which is '32a3'=[0, 0, 0, -11, -14]

from sage.all import Set, infinity, QQ, EllipticCurve, cm_j_invariants

# dictionary to hold the list of a-invariants of the minimal twist for
# all CM j-invariants other than 0 and 1728.  Keys are j-invariants.

min_quad_disc_CM_dict = {}
min_quad_disc_CM_dict[-12288000] = [0, 0, 1, -30, 63] # '27a4' CM disc = -27
min_quad_disc_CM_dict[54000] = [0, 0, 0, -15, 22] # '36a2' CM disc = -12
min_quad_disc_CM_dict[287496] = [0, 0, 0, -11, -14] # '32a3' CM disc = -16
min_quad_disc_CM_dict[16581375] = [1, -1, 0, -37, -78] # '49a2' CM disc = -28
min_quad_disc_CM_dict[-3375] = [1, -1, 0, -2, -1] # '49a1' CM disc = -7
min_quad_disc_CM_dict[8000] = [0, 1, 0, -3, 1] # '256a1' CM disc = -8
min_quad_disc_CM_dict[-32768] = [0, -1, 1, -7, 10] # '121b1' CM disc = -11
min_quad_disc_CM_dict[-884736] = [0, 0, 1, -38, 90] # '361a1' CM disc = -19
min_quad_disc_CM_dict[-884736000] = [0, 0, 1, -860, 9707] # '1849a1' CM disc = -43
min_quad_disc_CM_dict[-147197952000] = [0, 0, 1, -7370, 243528] # '4489a1' CM disc = -67
min_quad_disc_CM_dict[-262537412640768000] = [0, 0, 1, -2174420, 1234136692] # '26569a1' CM disc = -163
assert Set(min_quad_disc_CM_dict.keys()) + Set([0,1728]) == Set(cm_j_invariants(QQ))

def isog_key(E):
    return E.aplist(100)

def min_quad_twist(E):
    """
    Input: E, an elliptic curve over Q.
    Output: Emqt,D where Emqt is the minimal quadratic twist of E, its twist by D
    """
    j = E.j_invariant()
    try:
        Emqt = EllipticCurve(min_quad_disc_CM_dict[j])
        return Emqt, E.is_quadratic_twist(Emqt)
    except KeyError:
        pass
    assert j in [0,1728] or not E.has_cm()

    # find minimal quadratic twist of E at primes >3

    E = E.global_minimal_model()
    A = -27*E.c4()
    B = -54*E.c6()
    for p in A.gcd(B).support():
        eA = A.valuation(p)//2 if A else infinity
        eB = B.valuation(p)//3 if B else infinity
        e = min(eA, eB)
        if e>0:
            pe = p**e
            A /= pe**2
            B /= pe**3

    # now twist by -1,2,3

    tw = [1,-1,2,-2,3,-3,6,-6]
    EAB = EllipticCurve([A,B])
    # EAB may not be minimal but the quadratic twist method returns minimal models
    Elist = [EAB.quadratic_twist(t) for t in tw]

    # find minimal conductor:

    min_cond = min(E.conductor() for E in Elist)
    Elist = [E for E in Elist if E.conductor() == min_cond]
    if len(Elist)==1:
        return Elist[0], E.is_quadratic_twist(Elist[0])

    # sort into isogeny classes and pick out first

    from itertools import groupby
    Elist.sort(key = isog_key)
    Elist = list(list(next(groupby(Elist, key=isog_key)))[1])

    # This should leave 1 curve or 2 when E has CM
    n = len(Elist)
    if n==1:
        return Elist[0], E.is_quadratic_twist(Elist[0])

    assert n==2 and j in [0, 1728]
    Elist.sort(key=lambda e: e.discriminant().abs())
    return Elist[0], E.is_quadratic_twist(Elist[0])

# Function to add (or overwrite) the fields 'min_quad_twist_ainvs' and
# 'min_quad_twist_disc' from a curve record obtained from
# e.g. read_data(file_types=['curvedata'], ranges=['00000-09999'])

def add_mqt(record):
    from codec import decoders, encode
    E = EllipticCurve(decoders['ainvs'](record['ainvs']))
    Emqt, D = min_quad_twist(E)
    Emqt = encode(list(Emqt.ainvs()))
    D = encode(D)
    if Emqt != record['min_quad_twist_ainvs']:
        record['min_quad_twist_ainvs'] = Emqt
        record['min_quad_twist_disc'] = D
    return record

def add_mqt_range(r, base_dir):
    from files import read_data, write_curvedata
    dat = read_data(base_dir=base_dir, file_types=['curvedata'], ranges=[r])
    for lab, rec in dat.items():
        dat[lab] = add_mqt(rec)
    write_curvedata(dat, r, base_dir=base_dir)


