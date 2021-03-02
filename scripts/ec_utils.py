# elliptic curve utility functions

from sage.all import (Magma, pari, QQ,
                      Set, prime_range,
                      mwrank_get_precision, mwrank_set_precision)

from twoadic import init_2adic
from galrep import init_galrep

mwrank_saturation_precision = 1000 # 500 not enough for 594594bf2
mwrank_saturation_maxprime = 1000

GP = '/usr/local/bin/gp'

magma = Magma()
init_2adic(magma)
init_galrep(magma)

magma_count = 0
magma_count_max = 100 # restart magma after this number of uses
MagmaEffort = 100000   # 1000 not enough for  282203479a1 (rank 2)

def get_magma(verbose=False):
    global magma, magma_count, magma_count_max
    if magma_count == magma_count_max:
        if verbose:
            print("Restarting Magma")
        magma.quit()
        magma = Magma()
        init_2adic(magma)
        init_galrep(magma)
        magma_count = 1
    else:
        if verbose:
            print("Reusing Magma (count={})".format(magma_count))
        magma_count += 1
    return magma

# Assuming that E is known to have rank 1, returns a point on E
# computed by Magma's HeegnerPoint command
def magma_rank1_gen(E, mE):
    mP = mE.HeegnerPoint(nvals=2)[1]
    P = E([mP[i].sage() for i in [1,2,3]])
    return P

# Assuming that E is known to have rank 1, returns a point on E
# computed by GP's ellheegner() command
def pari_rank1_gen_old(E, stacksize=1024000000):
    from os import system, getpid, unlink
    f = 'tempfile-'+str(getpid())
    comm = "LD_LIBRARY_PATH=/usr/local/lib; echo `echo 'ellheegner(ellinit("+str(list(E.ainvs()))+"))' | %s -q -f -s %s` > %s;" % (GP,stacksize,f)
    system(comm)
    P = open(f).read()
    #print(P)
    P = open(f).read().partition("[")[2].partition("]")[0]
    P = P.replace("\xb1","") # needed for 497805u1
    #print(P)
    unlink(f)
    P = E([QQ(c) for c in P.split(',')])
    #print(P)
    return P

def pari_rank1_gen(E):
    return E(pari(E).ellheegner().sage())

def get_magma_gens(E, mE):
    MS = mE.MordellWeilShaInformation(RankOnly=True, Effort=MagmaEffort, nvals=3)
    rank_bounds = [r.sage() for r in MS[0]]
    gens = [E(P.Eltseq().sage()) for P in MS[1]]
    return rank_bounds, gens

def get_gens_mwrank(E):
    return E.gens(algorithm='mwrank_lib', descent_second_limit=15, sat_bound=2)

def get_rank1_gens(E, mE, verbose=0):
    if verbose:
        print(" - trying a point search...")
    gens = E.point_search(15)
    if gens:
        if verbose:
            print("--success: P = {}".format(gens[0]))
        return gens
    if verbose:
        print("--failed.  Trying pari's ellheegner...")
    gens = [pari_rank1_gen(E)]
    if gens:
        if verbose:
            print("--success: P = {}".format(gens[0]))
        return gens
    if verbose:
        print("--failed.  Trying Magma's HeegnerPoint...")
    try:
        gens = [magma_rank1_gen(E, mE)]
        if gens:
            if verbose:
                print("--success: P = {}".format(gens[0]))
            return gens
    except:
        pass
    if verbose:
        print("-- failed. Trying Magma...")
    rb, gens = get_magma_gens(E, mE)
    if gens:
        if verbose:
            print("--success: P = {}".format(gens[0]))
        return gens
    if verbose:
        print("--failed.  Trying mwrank...")
    return get_gens_mwrank(E)

def get_gens_simon(E):
    E.simon_two_descent(lim3=5000)
    return E.gens()

def get_gens(E, ar, verbose=0):
    if ar==0:
        return []
    mag = get_magma()
    mE = mag(E)
    if ar==1:
        if verbose>1:
            print("{}: a.r.=1, finding a generator".format(E.ainvs()))
        gens = get_rank1_gens(E,mE, verbose)
    else: # ar >=2
        if verbose>1:
            print("{}: a.r.={}, finding generators using Magma".format(E.ainvs(),ar))
        rb, gens = get_magma_gens(E, mE)
        if verbose>1:
            print("gens = {}".format(gens))

    # Now we have independent gens, and saturate them

    prec0=mwrank_get_precision()
    mwrank_set_precision(mwrank_saturation_precision)
    if verbose>1:
        print("Starting saturation (p<{})...".format(mwrank_saturation_maxprime))
    gens, index, reg = E.saturation(gens, max_prime=mwrank_saturation_maxprime)
    mwrank_set_precision(prec0)
    if verbose>1:
        print("... finished saturation (index {}, new reg={})".format(index, reg))

    return gens

# Given a matrix of isogenies and a list of points on the initial
# curve returns a# list of their images on each other curve.  The
# complication is that the isogenies will only exist when they have
# prime degree.

# Here we assume that the points in Plist are saturated, and only
# resaturate their images at primes up to the masimum prime dividing
# an isogeny degree.

def map_points(maps, Plist, verbose=0):
    ncurves = len(maps)
    if len(Plist)==0:
        return [[] for _ in range(ncurves)]
    if ncurves==1:
        return [Plist]
    if verbose>1:
        print("in map_points with degrees {}".format([[phi.degree() if phi else 0 for phi in r] for r in maps]))
    maxp = max([max([max(phi.degree().support(), default=0) if phi else 0 for phi in r], default=0) for r in maps], default=0)
    if verbose>1:
        print("  maxp = {}".format(maxp))

    Qlists = [Plist] + [[]]*(ncurves-1)
    nfill = 1
    for i in range(ncurves):
        if nfill==ncurves:
          break
        for j in range(1,ncurves):
            if not (maps[i][j] == 0) and Qlists[j]==[]:
                Qlists[j] = [maps[i][j](P) for P in Qlists[i]]
                nfill += 1
    # now we saturate the points just computed at all primes up to maxp
    prec0=mwrank_get_precision()
    mwrank_set_precision(mwrank_saturation_precision)
    for i in range(1,ncurves):
        E = Qlists[i][0].curve()
        if verbose>1:
          print("Saturating curve {} (maxp={})...".format(i, maxp))
        Qlists[i], n, r = E.saturation(Qlists[i], max_prime=maxp)
        if verbose>1:
          print("--saturation index was {}".format(n))
    mwrank_set_precision(prec0)
    return Qlists

# Find integral points in a fail-safe way uing both Sage and Magma,
# comparing, returning the union in all cases and outputting a warning
# message if they disagree.
def get_integral_points_with_sage(E, gens):
    return [P[0] for P in E.integral_points(mw_base=gens)]

def get_integral_points_with_magma(E, gens):
    mag = get_magma()
    mE = mag(E)
    xs = [E(P.Eltseq().sage())[0] for P in mE.IntegralPoints(FBasis=[mE(list(P)) for P in gens])]
    return xs

def get_integral_points(E, gens, verbose=True):
    x_list_magma = get_integral_points_with_magma(E, gens)
    x_list_sage = get_integral_points_with_sage(E, gens)
    if x_list_magma != x_list_sage:
        if verbose:
            print("Curve {}: \n".format(E.ainvs))
            print("Integral points via Magma: {}".format(x_list_magma))
            print("Integral points via Sage: {}".format(x_list_sage))
    x_list = list(Set(x_list_sage)+Set(x_list_magma))
    x_list.sort()
    return x_list

def get_modular_degree(E, label):
    degphi_magma = 0
    degphi_sympow = 0
    #return E.modular_degree(algorithm='sympow')
    try:
        degphi_magma = E.modular_degree(algorithm='magma')
    except RuntimeError:
        print("{}: degphi via magma failed".format(label))
        try:
            degphi_sympow = E.modular_degree(algorithm='sympow')
        except RuntimeError:
            print("{}: degphi via sympow failed".format(label))
    if degphi_magma:
        if degphi_sympow:
            if degphi_magma==degphi_sympow:
                return degphi_magma
            else:
                print("{}: degphi = {} from magma but {} from sympow!".format(label, degphi_magma, degphi_sympow))
                return degphi_magma
        else:
            return degphi_magma
    else:
        if degphi_sympow:
            return degphi_sympow
        else:
            print("{}: no success in computing degphi via magma or sympow".format(label))
            return 0

# Sage's E.aplist(100) returns a list of the Fourier coefficients for
# p<100.  For the aplist files, we want to replace the coefficient for
# p|N with the W-eigenvalue (the root number) and append the
# W-eigenvalues for p|N, p>100.  Not relevant for making LMFDBupload
# files.

def wstr(n,w):  # str(n) with enough spaces prepended to give width w
    a = str(n)
    if len(a)<w:
        a = ' '*(w-len(a)) + a
    return a

def my_ap(E,D,p):
    if p.divides(D):
        return E.root_number(p)
    return E.ap(p)

def my_ap_str(E,D,p):
    if p.divides(D):
        a = E.root_number(p)
        if a==1: 
            if p>23:
                return '  +'
            return ' +'
        if p>23:
            return '  -'
        return ' -'
    if p>23: 
        return wstr(E.ap(p),3)
    return wstr(E.ap(p),2)

def my_aplist(E):
    D = E.discriminant()
    ap = [my_ap_str(E,D,p) for p in prime_range(100)]
    qlist = D.support()
    for q in qlist:
        if q>100:
            if E.root_number(q)==1:
                ap.append('+('+str(q)+')')
            else:
                ap.append('-('+str(q)+')')
    return ' '.join(ap)

