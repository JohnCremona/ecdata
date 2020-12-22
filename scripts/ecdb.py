import os
import sys
HOME = os.getenv("HOME")
sys.path.append(os.path.join(HOME, 'ecdata', 'scripts'))
from sage.all import EllipticCurve, Integer, ZZ, QQ, Set, Magma, prime_range, factorial, mwrank_get_precision, mwrank_set_precision, srange, pari, EllipticCurve_from_c4c6, prod, copy
from red_gens import reduce_tgens, reduce_gens
from trace_hash import TraceHashClass
from files import make_line, datafile_columns, MATSCHKE_DIR, parse_line_label_cols, split
from codec import parse_int_list, point_to_weighted_proj, proj_to_point

from sage.databases.cremona import parse_cremona_label, class_to_int, cremona_letter_code

mwrank_saturation_precision = 300
mwrank_saturation_maxprime = 200000
GP = '/usr/local/bin/gp'

def print_data(outfile, code, ainvs, r, t):
    print("Code = {}".format(code))
    print("Curve = {}".format(EllipticCurve(ainvs)))
    print("rank = {}".format(r))
    print("torsion = {}".format(t))

def put_allcurves_line(outfile, N, cl, num, ainvs, r, t):
    line = ' '.join([str(N),cl,str(num),str(ainvs).replace(' ',''),str(r),str(t)])
    outfile.write(line+'\n')

def make_allcurves_lines(outfile, code, ainvs, r, t):
    E = EllipticCurve(ainvs)
    N, cl, n = parse_cremona_label(code)
    for i, F in enumerate(E.isogeny_class().curves):
        put_allcurves_line(outfile,N,cl,str(i+1),list(F.ainvs()),r,F.torsion_order())
    outfile.flush()

def process_curve_file(infilename, outfilename, use):
    infile = open(infilename)
    outfile = open(outfilename, mode='a')
    for L in infile.readlines():
        N, iso, num, ainvs, r, tor, d = L.split()
        code = N+iso+num
        N = int(N)
        num = int(num)
        r = int(r)
        tor = int(tor)
        ainvs = eval(ainvs)
        use(outfile, code, ainvs, r, tor)
    infile.close()
    outfile.close()

def liststr(l):
    return str(l).replace(' ','')

def shortstr(E):
    return liststr(list(E.ainvs()))

def shortstrlist(Elist):
    return str([list(F.ainvs()) for F in Elist]).replace(' ','')

def pointstr(P):
    P = list(P)
    z = P[1].denominator()
    P = [z*c for c in P]
    return '['+':'.join([str(c) for c in P])+']'
#    return str(P).replace('(','[').replace(')',']').replace(' ','')

# convert '[x:y:z]' to '[x/z,y/z]'
def pointPtoA(P):
    x,y,z = [Integer(c) for c in P[1:-1].split(":")]
    return [x/z,y/z]


def matstr(m):
    return str(list(m)).replace('(','[').replace(')',']').replace(' ','')

def mat_to_list_list(M):
    m,n = M.dimensions()
    return [[M[i][j] for j in range(n)] for i in range(m)]

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
    MS = mE.MordellWeilShaInformation(nvals=3)
    rank_bounds = [r.sage() for r in MS[0]]
    gens = [E(P.Eltseq().sage()) for P in MS[1]]
    return rank_bounds, gens

def get_gens_mwrank(E):
    return E.gens(algorithm='mwrank_lib', descent_second_limit=15, sat_bound=2)

def get_rank1_gens(E, mE, verbose=False):
    if verbose:
        print(" - trying a point search...")
    gens = E.point_search(15)
    if gens:
        if verbose:
            print("--success: P = {}".format(gens[0]))
        return gens
    if verbose:
        print("--failed.  Trying a pari's ellheegner...")
    gens = [pari_rank1_gen(E)]
    if gens:
        if verbose:
            print("--success: P = {}".format(gens[0]))
        return gens
    if verbose:
        print("--failed.  Trying a Magma's HeegnerPoint...")
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

def get_gens(E, ar, verbose=False):
    if ar==0:
        return []
    mag = Magma()
    mE = mag(E)
    if ar==1:
        if verbose:
            print("a.r.=1, finding a generator")
        gens = get_rank1_gens(E,mE, verbose)
    else: # ar >=2
        if verbose:
            print("a.r.={}, finding generators using Magma".format(ar))
        rb, gens = get_magma_gens(E, mE)
    mag.quit()
    return gens

# Given a matrix of isogenies and a list of points on the initial
# curve returns a# list of their images on each other curve.  The
# complication is that the isogenies will only exist when they have
# prime degree.
def map_points(maps, Plist):
    ncurves = len(maps)
    if len(Plist)==0:
        return [[] for _ in range(ncurves)]
    if ncurves==1:
        return [Plist]
    Qlists = [Plist] + [[]]*(ncurves-1)
    nfill = 1
    for i in range(ncurves):
        if nfill==ncurves:
            return Qlists
        for j in range(1,ncurves):
            if not (maps[i][j] == 0) and Qlists[j]==[]:
                Qlists[j] = [maps[i][j](P) for P in Qlists[i]]
                nfill += 1


# Find integral points in a fail-safe way uing both Sage and Magma,
# comparing, returning the union in all cases and outputting a warning
# message if they disagree.
def get_integral_points_with_sage(E, gens):
    return [P[0] for P in E.integral_points(mw_base=gens)]

def get_integral_points_with_magma(E, gens):
    mag = Magma()
    mE = mag(E)
    xs = [E(P.Eltseq().sage())[0] for P in mE.IntegralPoints(FBasis=[mE(list(P)) for P in gens])]
    mag.quit()
    return xs

def get_integral_points(E, gens, verbose=True):
    x_list_magma = get_integral_points_with_magma(E, gens)
    x_list_sage = get_integral_points_with_sage(E, gens)
    if x_list_magma != x_list_sage:
        if verbose:
            print("Curve {} = {}: \n".format(E.ainvs))
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


# Given a filename like curves.000-999, read the data in the file,
# compute the isogeny class for each curve, and output (1)
# allcurves.000-999, (2) allisog.000-999 (with the same suffix).  Also
# compute the bsd data for each curve and output (3) allbsd.000-999,
# (4) allgens.000-999 (with the same suffix), (5) degphi.000-999, (6)
# intpts.000-999, (7) alldegphi.000-999

# Version to compute gens & torsion gens too
def make_datafiles(infilename, mode='w', verbose=False, prefix="t"):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    allcurvefile = open(prefix+"allcurves."+suf, mode=mode)
    allisogfile = open(prefix+"allisog."+suf, mode=mode)
    allbsdfile = open(prefix+"allbsd."+suf, mode=mode)
    allgensfile = open(prefix+"allgens."+suf, mode=mode)
    degphifile = open(prefix+"degphi."+suf, mode=mode)
    alldegphifile = open(prefix+"alldegphi."+suf, mode=mode)
    apfile = open(prefix+"aplist."+suf, mode=mode)
    intptsfile = open(prefix+"intpts."+suf, mode=mode)
    for L in infile.readlines():
        if verbose: print("="*72)
        N, cl, num, ainvs, r, tor, d = L.split()
        E = EllipticCurve(eval(ainvs))
        r=int(r)
        # Compute the isogeny class
        Cl = E.isogeny_class()
        Elist = Cl.curves
        mat = Cl.matrix()
        maps = Cl.isogenies()
        ncurves = len(Elist)
        print("class {} (rank {}) has {} curve(s)".format(N+cl,r,ncurves))
        line = ' '.join([str(N),cl,str(1),ainvs,shortstrlist(Elist),matstr(mat)])
        allisogfile.write(line+'\n')
        if verbose: print("allisogfile: {}".format(line))
        # compute BSD data for each curve
        torgroups = [F.torsion_subgroup() for F in Elist]
        torlist = [G.order() for G in torgroups]
        torstruct = [list(G.invariants()) for G in torgroups]
        torgens = [[P.element() for P in G.gens()] for G in torgroups]
        cplist  = [F.tamagawa_product() for F in Elist]
        omlist  = [F.real_components()*F.period_lattice().real_period() for F in Elist]

        Lr1 = E.pari_curve().ellanalyticrank()[1].sage() / factorial(r)
        if r==0:
            genlist = [[] for F in Elist]
            reglist = [1 for F in Elist]
        else:
            Plist = get_gens(E, r, verbose)
            genlist = map_points(maps,Plist)
            prec0=mwrank_get_precision()
            mwrank_set_precision(mwrank_saturation_precision)
#            genlist = [Elist[i].saturation(genlist[i], max_prime=mwrank_saturation_maxprime)[0] for i in range(ncurves)]
            if verbose: print("genlist (before saturation) = {}".format(genlist))
            genlist = [Elist[i].saturation(genlist[i])[0] for i in range(ncurves)]
            if verbose: print("genlist (before reduction) = {}".format(genlist))
            genlist = [Elist[i].lll_reduce(genlist[i])[0] for i in range(ncurves)]
            mwrank_set_precision(prec0)
            if verbose: print("genlist  (after reduction)= {}".format(genlist))
            reglist = [Elist[i].regulator_of_points(genlist[i]) for i in range(ncurves)]
        shalist = [Lr1*torlist[i]**2/(cplist[i]*omlist[i]*reglist[i]) for i in range(ncurves)]
        squares = [n*n for n in srange(1,100)]
        for i,s in enumerate(shalist):
            if not round(s) in squares:
                print("bad sha value %s for %s" % (s,str(N)+cl+str(i+1)))
                print("Lr1 = %s" % Lr1)
                print("#t  = %s" % torlist[i])
                print("cp  = %s" % cplist[i])
                print("om  = %s" % omlist[i])
                print("reg = %s" % reglist[i])
                return Elist[i]
        if verbose: print("shalist = {}".format(shalist))

        # compute modular degrees
        # degphilist = [e.modular_degree(algorithm='magma') for e in Elist]
        # degphi = degphilist[0]
        # NB The correctness of the following relies on E being optimal!
        try:
            degphi = E.modular_degree(algorithm='magma')
        except RuntimeError:
            degphi = 0
        degphilist1 = [degphi*mat[0,j] for j in range(ncurves)]
        degphilist = degphilist1
        if verbose: print("degphilist = {}".format(degphilist))

        # compute aplist for optimal curve only
        aplist = my_aplist(E)
        if verbose: print("aplist = {}".format(aplist))

        # Compute integral points (x-coordinates)
        intpts = [get_integral_points(Elist[i],genlist[i]) for i in range(ncurves)]
        #if verbose: print("intpts = {}".format(intpts))
        for i, Ei, xs in zip(range(ncurves),Elist,intpts):
            print("{}{}{} = {}: intpts = {}".format(N,cl,(i+1),Ei.ainvs(),xs))

        # Output data for optimal curves

        # aplist
        line = ' '.join([N,cl,aplist])
        if verbose: print("aplist: {}".format(line))
        apfile.write(line+'\n')

        # degphi
        line = ' '.join([N,cl,'1',str(degphi),str(Set(degphi.prime_factors())).replace(' ',''),shortstr(E)])
        if verbose: print("degphifile: {}".format(line))
        degphifile.write(line+'\n')

        # Output data for each curve
        for i in range(ncurves):
            # allcurves
            line = ' '.join([N,cl,str(i+1),shortstr(Elist[i]),str(r),str(torlist[i])])
            allcurvefile.write(line+'\n')
            if verbose: print("allcurvefile: {}".format(line))

            # allbsd
            line = ' '.join([N,cl,str(i+1),shortstr(Elist[i]),str(r),str(torlist[i]),str(cplist[i]),str(omlist[i]),str(Lr1),str(reglist[i]),str(shalist[i])])
            allbsdfile.write(line+'\n')
            if verbose: print("allbsdfile: {}".format(line))

            # allgens (including torsion gens, listed last)
            line = ' '.join([str(N),cl,str(i+1),shortstr(Elist[i]),str(r)]
                            + [liststr(torstruct[i])]
                            + [pointstr(P) for P in genlist[i]]
                            + [pointstr(P) for P in torgens[i]]
                            )
            allgensfile.write(line+'\n')
            if verbose:
                print("allgensfile: {}".format(line))
            # intpts
            line = ''.join([str(N),cl,str(i+1)]) + ' ' + shortstr(Elist[i]) + ' ' + liststr(intpts[i])
            intptsfile.write(line+'\n')
            if verbose: print("intptsfile: {}".format(line))

            # alldegphi
            line = ' '.join([str(N),cl,str(i+1),shortstr(Elist[i]),liststr(degphilist[i])])
            alldegphifile.write(line+'\n')
            if verbose: print("alldegphifile: {}".format(line))

    infile.close()
    allbsdfile.close()
    allgensfile.close()
    allcurvefile.close()
    allisogfile.close()
    degphifile.close()
    intptsfile.close()
    apfile.close()


# Read allgens file (with torsion) and output paricurves file
#
def make_paricurves(infilename, mode='w', verbose=False, prefix="t"):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    paricurvesfile = open(prefix+"paricurves."+suf, mode=mode)
    for L in infile.readlines():
        N, cl, num, ainvs, r, gens = L.split(' ',5)
        if int(r)==0:
            gens = []
        else:
            gens = gens.split()[1:1+int(r)] # ignore torsion struct and gens
            gens = [pointPtoA(P) for P in gens]

        line = '[' + ', '.join(['"'+N+cl+num+'"',ainvs,str(gens).replace(' ','')]) + ']'
        paricurvesfile.write(line+'\n')
    infile.close()
    paricurvesfile.close()

# Create alldegphi files from allcurves files:

def make_alldegphi(infilename, mode='w', verbose=False, prefix="t"):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    alldegphifile = open(prefix+"alldegphi."+suf, mode=mode)
    for L in infile.readlines():
        N, cl, num, ainvs, rest = L.split(' ',4)
        label = "".join([N,cl,num])
        E = EllipticCurve(eval(ainvs))
        degphi = get_modular_degree(E, label)
        line = ' '.join([str(N),cl,str(num),shortstr(E),liststr(degphi)])
        alldegphifile.write(line+'\n')
        if verbose: print("alldegphifile: {}".format(line))
    infile.close()
    alldegphifile.close()

# Create manin constant files from allcurves files:
#
# NB We assume that curve 1 in each class is optimal with constant 1
#
# Use this on a range where we have established that the optimal curve
# is #1.  Otherwise the C++ program h1pperiods outputs a file
# opt_man.<range> which includes what can be deduced about optimality
# (without a full check) and also outputs Manin constants conditional
# on the optimal curve being #1.

# The infilename should be an allcurves file, e.g. run in an ecdata
# directory and use infilename=allcurves/allcurves.250000-259999.  The
# part after the "." (here 250000-259999) will be used as suffix to
# the output file.

# Some classes may be output as "special": two curves in the class
# linked by a 2-isogeny with the first lattice having positive
# discriminant and the second Manin constant=2.  In several cases this
# has been an indication that the wrong curve has been tagged as
# optimal.

def make_manin(infilename, mode='w', verbose=False, prefix=""):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    allmaninfile = open(prefix+"opt_man."+suf, mode=mode)
    allisogfile = open("allisog/allisog."+suf)
    last_class = ""
    manins = []
    area0 = 0    # else pyflakes objects
    degrees = [] # else pyflakes objects
    for L in infile.readlines():
        N, cl, num, ainvs, rest = L.split(' ',4)
        this_class = N+cl
        num = int(num)
        E = EllipticCurve(eval(ainvs))
        lattice1 = None
        if this_class == last_class:
            deg = degrees[int(num)-1]
            #assert ideg==deg
            area = E.period_lattice().complex_area()
            mc = round((deg*area/area0).sqrt())
            # print("{}{}{}: ".format(N,cl,num))
            # print("degree =     {}".format(deg))
            # print("area   =     {}".format(area))
            # print("area ratio = {}".format(area/area0))
            # print("c^2 =        {}".format(deg*area/area0))
            manins.append(mc)
            if num==3:
                lattice1 = None
            elif not (num==2 and deg==2 and mc==2):
                lattice1 = None
            if lattice1:
                print("Class {} is special".format(lattice1))
        else: # start a new class
            if manins and verbose: # else we're at the start
                print("class {} has Manin constants {}".format(last_class,manins))
            isogmat = allisogfile.readline().split()[-1]
            isogmat = eval(isogmat)
            degrees = isogmat[0]
            if verbose:
                print("class {}".format(this_class))
                print("isogmat: {}".format(isogmat))
                print("degrees: {}".format(isogmat))
            area0 = E.period_lattice().complex_area()
            manins = [1]
            if num==1  and len(isogmat)==2 and isogmat[0][1]==2 and E.discriminant()>0:
                lattice1 = this_class
            last_class = this_class

        # construct output line for this curve
        #
        # SPECIAL CASE 990h
        #
        if N==90 and cl=='h':
            opt = str(int(num==3))
            manins = [1,1,1,1]
        else:
            opt = str(int(num==1))
        mc = str(manins[-1])
        line = ' '.join([str(N),cl,str(num),shortstr(E),opt,mc])

        allmaninfile.write(line+'\n')
        if verbose: print("allmaninfile: {}".format(line))
        # if int(N)>100:
        #     break
    if manins and verbose: # else we're at the start
        # This output line is the only reason we keep the list of all m.c.s
        print("class {} has Manin constants {}".format(last_class,manins))
    infile.close()
    allmaninfile.close()

def make_opt_input(N):
    """Parse a file in optimality/ and produce an input file for
    runh1firstx1 which runs h1first1.

    Find lines containing "possible optimal curves", extract the
    isogeny class label, convert iso class code to a number, and
    output.  One line is output for each level N, of the form

    N i_1 i_2 ... i_k 0

    where N is the level and i_1,...,i_k are the numbers (from 1) of
    the newforms / isogeny classes where we do not know which curve is
    optimal.

    For example the input lines starting with 250010 are

    250010b: c=1; optimal curve is E1
    250010d: c=1; 3 possible optimal curves: E1 E2 E3

    and produce the output line

    250010 2 4 0

    while the lines starting with 250020 are

    250020b: c=1; optimal curve is E1
    250020d: c=1; optimal curve is E1

    and produce no output
    """
    outfile = "optx.{0:02d}".format(N)
    o = open(outfile, 'w')
    n = 0
    lastN = 0
    for L in open("optimality/optimality.{0:02d}".format(N)):
        if "possible" in L:
            N, c, i = parse_cremona_label(L.split()[0][:-1])
            if N==lastN:
                o.write(" {}".format(1+class_to_int(c)))
            else:
                if lastN!=0:
                    o.write(" 0\n")
                o.write("{} {}".format(N, 1+class_to_int(c)))
                lastN = N
                n += 1
    o.write(" 0\n")
    n += 1
    o.close()
    print("wrote {} lines to {}".format(n,outfile))

def check_sagedb(N1, N2, a4a6bound=100):
    """Sanity check that Sage's elliptic curve database contains all
    curves [a1,a2,a3,a4,a6] with a1,a3 in [0,1], a2 in [-1,0,1], a4,a6
    in [-100..100] and conductor in the given range.

    Borrowed from Bill Allombert (who found a missing curve of
    conductor 406598 this way).
    """
    from sage.all import CremonaDatabase
    CDB = CremonaDatabase()
    def CDB_curves(N):
        return [c[0] for c in CDB.allcurves(N).values()]
    Nrange = srange(N1,N2+1)
    ncurves = 12*(2*a4a6bound+1)**2
    print("Testing {} curves".format(ncurves))
    n=0
    for a1 in range(2):
        for a2 in range(-1,2):
            for a3 in range(2):
                for a4 in range(-a4a6bound,a4a6bound+1):
                    for a6 in range(-a4a6bound,a4a6bound+1):
                        ai = [a1,a2,a3,a4,a6]
                        n += 1
                        if n%1000==0:
                            print("test #{}/{}".format(n,ncurves))
                        try:
                            E = EllipticCurve(ai).minimal_model()
                        except ArithmeticError: #singular curve
                            continue
                        N = E.conductor()
                        if not N in Nrange:
                            continue
                        if not list(E.ainvs()) in CDB_curves(N):
                            print("Missing N={}, ai={}".format(N,ai))

# From here on, functions to process a set of "raw" curves, to create
# all the data files needed for LMFDB upload

# Assume that in directory CURVE_DIR there is a file 'curves.NN'
# containing one curve per line (valid formats: as in
# process_raw_curves() below).  For correct isogeny class labelling it
# is necessary that for every conductor present, at least one curve
# from each isogeny class is present.  The input list need not be
# closed under isogeny: if it is not then the 'allcurves.NN' file
# will contain more curves than the input file, otherwise they will
# contain the same curves. The file suffix 'NN' can be anything
# (nonempty) to identify the dataset, e.g. it could be a single
# conductor, or a range of conductors.
#
# Step 1: sort and label.  Run sage from ecdata/scripts:
#
# sage: %runfile ecdb.py
# sage: NN = '2357' # (say)
# sage: curves_file = 'curves.{}'.format(NN)
# sage: allcurves_file = 'allcurves.{}'.format(NN)
# sage: process_raw_curves(curves_file, allcurves_file, base_dir=CURVE_DIR)
#
# Step 2: compute most of the data, run also in ecdb/scripts after reading ecdb.py:
#
# sage: read_write_data(allcurves_file, CURVE_DIR)
#
# We now have files CURVE_DIR/ft/ft.NN for ft in ['curvedata', 'classdata', 'intpts', 'alldegphi'].
#
# Step 3: run magma scripts to compute galrep and 2adic data. From ecdata/scripts:
# ./make_galrep.sh NN CURVE_DIR
#
# This creates CURVE_DIR/ft/ft.NN for ft in ['galrep', '2adic'].
#
# Step 4: create 6 files in UPLOAD_DIR (default ${HOME}/ecq-upload):
#
# ec_curvedata.NN, ec_classdata.NN, ec_localdata.NN,
# ec_mwbsd.NN, ec_galrep.NN, ec_2adic.NN
#
# In ecdata/scripts run sage:
# sage: %runfile files.py
# sage: data = read_data(CURVE_DIR, file_types=main_file_types, ranges=[NN])
# sage: make_all_upload_files(data, tables = main_tables, NN=NN)
#
# ready for upload to LMFDB's db.ec_curvedata (etc) using
#
# db.ec_curvedata.update_from_file()


def curve_from_inv_string(s):
    invs = parse_int_list(s)
    if len(invs)==5:
        E = EllipticCurve(invs).minimal_model()
    elif len(invs)==2:
        E = EllipticCurve_from_c4c6(*invs).minimal_model()
    else:
        raise ValueError("{}: invariant list must have length 2 or 5".format(s))
    return E

def process_raw_curves(infilename, outfilename, base_dir='.', split_by_N=False, verbose=1):
    """File infilename should contain one curve per line, with
    a-ainvariants or c-invariants as a list (with no internal spaces),
    optionally preceded by the conductor.

    Sample lines (all defining the same curve):

    [0,-1,1,-10,-20]
    11 [0,-1,1,-10,-20]
    [496,20008]
    11 [496,20008]

    We do not assume minimal or reduced models, and if the conductor
    is given it is checked for consistency.

    For each input curve we check whether an isogenous curve has
    already been processed; if not, we compute its isogeny class and
    process all of them.

    OUTPUT: one line per curve, sorted by conductor and then by
    isogeny class.  Each line has the format

    N class_code number ainvs lmfdb_number

    where "number" is the index of the curve in its class (counting
    from 1), sorted by Faltings heights, with the lex. order of
    a-invariants as a tie-breaker.

    If split_by_N is False, the output will be all in one file.  If
    True, there will be one output file per conductor N, whose name is
    allcurves file with suffix ".{N}".

    """
    # allcurves will have conductors as keys, values lists of lists of
    # ainvs, subdivided by isogeny class
    allcurves = {}
    ncurves = ncurves_complete = 0

    with open(os.path.join(base_dir,infilename)) as infile:
        for L in infile:
            data = L.split()
            assert len(data) in [1,2]
            invs = data[-1]
            E = curve_from_inv_string(invs)
            ainvs = E.ainvs()
            N = E.conductor()
            if len(data)==2:
                N1 = ZZ(data[0])
                if N1!=N:
                    raise ValueError("curve with invariants {} has conductor {}, not {} as input".format(invs, N, N1))
            ncurves += 1
            if ncurves%1000==0 and verbose:
                print("{} curves read from {} so far...".format(ncurves, infilename))

            # now we have the curve E, check if we have seen it (or an isogenous curve) before
            repeat = False

            if N in allcurves:
                if any(ainvs in c for c in allcurves[N]):
                    repeat = True
                    if verbose>1:
                        print("Skipping {} of conductor {}: repeat".format(ainvs,N))
                else:
                    if verbose>1:
                        print("New curve {} of conductor {}".format(ainvs,N))
            else:
                allcurves[N] = []
            if repeat:
                continue # to next input line
            newcurves = [E2.ainvs() for E2 in E.isogeny_class().curves]
            ncurves_complete += len(newcurves)
            allcurves[N].append(newcurves)

    print("{} curves read from {}".format(ncurves, infilename))
    if ncurves != ncurves_complete:
        print("input curves not closed under isogeny! completed list contains {} curves".format(ncurves_complete))

    conductors = list(allcurves.keys())
    conductors.sort()

    for N in conductors:
        C = allcurves[N]
        ncl = len(C)
        nc = sum([len(cl) for cl in C])
        print("Conductor {}:\t{} isogeny classes,\t{} curves".format(N,ncl,nc))

    # Now we work out labels and output a new file

    def curve_key_LMFDB(E): # lex order of ainvs
        return E.ainvs()

    def curve_key_Faltings(E): # Faltings height, with tie-break
        return [-E.period_lattice().complex_area(), E.ainvs()]

    def output_one_conductor(N, allcurves_N, outfile):
        nNcl = nNcu = 0
        sN = str(N)
        # construct the curves from their ainvs:
        CC = [[EllipticCurve(ai) for ai in cl] for cl in allcurves_N]
        nap=100
        ok = False
        while not ok:
            aplists = dict([(cl[0],cl[0].aplist(100,python_ints=True)) for cl in CC])
            aps = list(aplists.values())
            ok = (len(list(Set(aps))) == len(aps))
            if not ok:
                print("{} ap not enough for conductor {}, increasing...".format(nap,N))
            nap += 100
        # sort the isogeny classes:
        CC.sort(key=lambda cl: aplists[cl[0]])
        for ncl, cl in enumerate(CC):
            nNcl += 1
            class_code = cremona_letter_code(ncl)
            # sort the curves in two ways (LMFDB, Faltings)
            cl.sort(key = curve_key_LMFDB)
            cl_Faltings = copy(cl)
            cl_Faltings.sort(key=curve_key_Faltings)
            class_size = len(cl)
            for nE_F, E in enumerate(cl_Faltings):
                nNcu += 1
                ainvs = shortstr(E)
                nE_L = cl.index(E)
                line = " ".join([sN,class_code,str(class_size),str(nE_F+1),str(nE_L+1),ainvs])
                #print(line)
                outfile.write(line+"\n")
        return nNcu, nNcl

    if split_by_N:
        for N in conductors:
            outfilename_N = ".".join([outfilename,str(N)])
            with open(os.path.join(base_dir,outfilename_N), 'w') as outfile:
                nNcu, nNcl = output_one_conductor(N, allcurves[N], outfile)
                print("N={}: {} curves in {} classes output to {}".format(N,nNcu,nNcl,outfilename_N))

    else:
        with open(os.path.join(base_dir,outfilename), 'w') as outfile:
            for N in conductors:
                nNcu, nNcl = output_one_conductor(N, allcurves[N], outfile)
                print("N={}: {} curves in {} classes output to {}".format(N,nNcu,nNcl,outfilename))

def parse_allgens_line_simple(line):
    r"""
    Parse one line from an allgens file

    Lines contain 6+t+r fields (columns)

    conductor iso number ainvs r torsion_structure <tgens> <gens>

    where:

    torsion_structure is a list of t = 0,1,2 ints
    <tgens> is t fields containing torsion generators
    <gens> is r fields containing generators mod torsion

    """
    label, record = parse_line_label_cols(line, 3, True)
    E = EllipticCurve(record['ainvs'])
    data = split(line)
    rank = int(data[4])
    record['gens'] = [proj_to_point(gen, E) for gen in data[6:6 + rank]]
    return label,  record

def make_new_data(infilename, base_dir, PRECISION=100, verbose=1, allgensfilename=None):
    alldata = {}
    nc = 0
    with open(os.path.join(base_dir, infilename)) as infile:
        for L in infile:
            nc += 1
            sN, isoclass, class_size, number, lmfdb_number, ainvs = L.split()
            iso = ''.join([sN,isoclass])
            label = ''.join([iso,number])
            lmfdb_number = int(lmfdb_number)
            lmfdb_isoclass = isoclass
            lmfdb_iso = '.'.join([sN,isoclass])
            lmfdb_label = ''.join([lmfdb_iso,number])
            iso_nlabel = class_to_int(isoclass)
            number = int(number)
            class_size = int(class_size)
            ainvs = parse_int_list(ainvs)
            N = ZZ(sN)
            bad_p = N.prime_factors() # will be sorted

            record = {
                'label': label,
                'isoclass': isoclass,
                'iso': iso,
                'number': number,
                'iso_nlabel': iso_nlabel,
                'lmfdb_number': lmfdb_number,
                'lmfdb_isoclass': lmfdb_isoclass,
                'lmfdb_iso': lmfdb_iso,
                'lmfdb_label':lmfdb_label,
                'faltings_index': number,
                'class_size': class_size,
                'ainvs': ainvs,
                'conductor': N,
                'bad_primes': bad_p,
                'num_bad_primes': len(bad_p),
                }
            alldata[label] = record

    if verbose:
        print("{} curves read from {}".format(nc, infilename))

    if allgensfilename:
        print("Reading from {}".format(allgensfilename))
        n = 0
        with open(os.path.join(base_dir, allgensfilename)) as allgensfile:
            for L in allgensfile:
                n+=1
                label, record = parse_allgens_line_simple(L)
                if label in alldata:
                    alldata[label].update(record)
                else:
                    print("ignoring allgens data for {}".format(label))
                if n%1000==0:
                    print("Read {} curves from {}".format(n,allgensfilename))


    for label, record in alldata.items():
        if verbose:
            print("Processing {}".format(label))
        N = record['conductor']
        iso = record['iso']
        number = record['number']
        first = (number==1) # tags first curve in each isogeny class
        ncurves = record['class_size']
        if first:
            alllabels = [iso+str(n+1) for n in range(ncurves)]
            allcurves = [EllipticCurve(alldata[lab]['ainvs']) for lab in alllabels]
            E = allcurves[0]
        else:
            record1 = alldata[iso+'1']
            E = allcurves[number-1]
        assert N==E.conductor()
        if verbose>1:
            print("E = {}".format(E))

        record['jinv'] = j = E.j_invariant()
        record['potential_good_reduction'] = (j.denominator()==1)
        record['signD'] = int(E.discriminant().sign())
        record['cm'] = int(E.cm_discriminant()) if E.has_cm() else 0

        if first:
            record['aplist'] = E.aplist(100,python_ints=True)
            record['anlist'] = E.anlist(20,python_ints=True)
            record['trace_hash'] = TraceHashClass(record['iso'], E)
        else:
            record['aplist'] = record1['aplist']
            record['anlist'] = record1['anlist']
            record['trace_hash'] = record1['trace_hash']

        if verbose>1:
            print("aplist done")
        local_data = [{'p': int(ld.prime().gen()),
                       'ord_cond':int(ld.conductor_valuation()),
                       'ord_disc':int(ld.discriminant_valuation()),
                       'ord_den_j':int(max(0,-(E.j_invariant().valuation(ld.prime().gen())))),
                       'red':int(ld.bad_reduction_type()),
                       'rootno':int(E.root_number(ld.prime().gen())),
                       'kod':ld.kodaira_symbol()._pari_code(),
                       'cp':int(ld.tamagawa_number())}
                      for ld in E.local_data()]

        record['tamagawa_numbers']    = cps = [ld['cp'] for ld in local_data]
        record['kodaira_symbols']           = [ld['kod'] for ld in local_data]
        record['reduction_types']           = [ld['red'] for ld in local_data]
        record['root_numbers']              = [ld['rootno'] for ld in local_data]
        record['conductor_valuations'] = cv = [ld['ord_cond'] for ld in local_data]
        record['discriminant_valuations']   = [ld['ord_disc'] for ld in local_data]
        record['j_denominator_valuations']  = [ld['ord_den_j'] for ld in local_data]

        record['semistable'] = all([v==1 for v in cv])
        record['tamagawa_product'] = tamprod = prod(cps)
        if verbose>1:
            print("local data done")

        T = E.torsion_subgroup()
        tgens = [P.element() for P in T.gens()]
        tgens.sort(key=lambda P:P.order())
        tgens = reduce_tgens(tgens)
        tor_struct = [P.order() for P in tgens]
        record['torsion_generators'] = [point_to_weighted_proj(gen) for gen in tgens]
        record['torsion_structure'] = tor_struct
        record['torsion'] = torsion = prod(tor_struct)
        record['torsion_primes'] = [int(p) for p in Integer(torsion).support()]
        if verbose>1:
            print("torsion done")

        if first: # else add later to avoid recomputing a.r.

            # Analytic rank and special L-value
            ar,sv = E.pari_curve().ellanalyticrank(precision=PRECISION)
            record['analytic_rank'] = ar = ar.sage()
            record['special_value'] = sv = sv.sage()/factorial(ar)

            # compute isogenies so we can map points from #1 to the rest:
            cl = E.isogeny_class(order=tuple(allcurves))
            record['isogeny_matrix'] = mat = mat_to_list_list(cl.matrix())
            record['class_deg'] = max(max(r) for r in mat)
            record['isogeny_degrees'] = mat[0]
            isogenies = cl.isogenies() # for mapping points later
            record['degree'] = degphi = get_modular_degree(E, label)
            record['degphilist'] = degphilist = [degphi*mat[0][j] for j in range(ncurves)]
            if verbose>1: print("degphilist = {}".format(degphilist))

        else:

            record['analytic_rank'] = ar = record1['analytic_rank']
            record['special_value'] = sv = record1['special_value']
            record['class_deg'] = record1['class_deg']
            record['isogeny_degrees'] = record1['isogeny_matrix'][number-1]
            record['degree'] = record1['degphilist'][number-1]

        if verbose>1:
            print("analytic rank done: {}".format(ar))
        if ar==0:
            record['gens'] = gens = []
            record['regulator'] = reg = 1
            record['ngens'] = 0
            record['heights'] = []
            record['rank'] = 0
            record['rank_bounds'] = [0,0]

        else: # positive rank
            if first:
                if verbose>1:
                    print("{}: an.rk.={}, finding generators".format(label, ar))
                if 'gens' in record:
                    gens = record['gens']
                    if verbose>1:
                        print("..already have gens {}".format(gens))
                else:
                    gens = get_gens(E, ar, verbose)
                ngens = len(gens)
                if ngens <ar:
                    print("{}: analytic rank = {} but we only found {} generators".format(label,ar,ngens))
                else:
                    if verbose>1:
                        print("...done, generators {}".format(gens))
                record['rank_bounds'] = [ngens, ar]
                record['rank'] = ngens if ngens == ar else None
                # so the other curves in the class know their gens:
                record['allgens'] = map_points(isogenies, gens)
            else:
                gens = record1['allgens'][number-1]
                record['rank'] = record1['rank']
                record['rank_bounds'] = record1['rank_bounds']

            record['gens'] = [point_to_weighted_proj(gen) for gen in gens]
            prec0=mwrank_get_precision()
            mwrank_set_precision(mwrank_saturation_precision)
            maxp = 0 if first else max(record['class_deg'].support())
            if verbose>1:
                print("Starting saturation (maxp={})...".format(maxp))
            gens, index, reg = E.saturation(gens, max_prime=maxp)
            mwrank_set_precision(prec0)
            if verbose>1:
                print("... finished saturation (index {}, new reg={}), reducing...".format(index, reg))
            gens, tgens = reduce_gens(gens, tgens, False, label)
            record['heights'] = heights = [P.height(precision=PRECISION) for P in gens]
            reg = E.regulator_of_points(gens, PRECISION)
            record['ngens'] = len(gens)
            record['regulator'] = reg
            if verbose>1:
                print("... finished reduction, gens are now {}, heights {}, reg={}".format(gens,heights,reg))

        L = E.period_lattice()
        record['real_period'] = om = L.omega(prec=PRECISION) # includes #R-components factor
        record['area'] = A = L.complex_area(prec=PRECISION)
        record['faltings_height'] = -A.log()/2

        # Analytic Sha
        sha_an = sv*torsion**2 / (tamprod*reg*om)
        sha = sha_an.round()
        assert sha>0
        assert sha.is_square()
        assert ((sha-sha_an).abs() < 1e-10)
        record['sha_an'] = sha_an
        record['sha'] = int(sha)
        record['faltings_ratio'] = 1 if first else (record1['area']/A).round()

        if verbose>1:
            print(" -- getting integral points...")
        record['intpts'] = get_integral_points(E, gens)
        if verbose>1:
            print(" ...done: {}".format(record['intpts']))


        Etw, Dtw = E.minimal_quadratic_twist()
        if Etw.conductor()==N:
            record['min_quad_twist_ainvs'] = record['ainvs']
            record['min_quad_twist_disc']  = 1
        else:
            record['min_quad_twist_ainvs'] = [int(a) for a in Etw.ainvs()]
            record['min_quad_twist_disc']  = int(Dtw)

        if verbose:
            print("Finished processing {}".format(label))
            if verbose>1:
                print("data for {}:\n{}".format(label,record))

    return alldata


def write_curvedata(data, r, base_dir=MATSCHKE_DIR):
    r"""
    Write file base_dir/curvedata/curvedata.<r>

    """
    cols = datafile_columns['curvedata']
    filename = os.path.join(base_dir, 'curvedata', 'curvedata.{}'.format(r))
    print("Writing data to {}".format(filename))
    n = 0
    with open(filename, 'w') as outfile:
        for label, record in data.items():
            line = make_line(record, cols)
            outfile.write(line +"\n")
            n += 1
    print("{} lines written to {}".format(n, filename))

def write_classdata(data, r, base_dir=MATSCHKE_DIR):
    r"""
    Write file base_dir/classdata/classdata.<r>

    """
    cols = datafile_columns['classdata']
    filename = os.path.join(base_dir, 'classdata', 'classdata.{}'.format(r))
    print("Writing data to {}".format(filename))
    n = 0
    with open(filename, 'w') as outfile:
        for label, record in data.items():
            if record['number']==1:
                line = make_line(record, cols)
                outfile.write(line +"\n")
                n += 1
    print("{} lines written to {}".format(n, filename))

def write_intpts(data, r, base_dir=MATSCHKE_DIR):
    r"""
    Write file base_dir/intpts/intpts.<r>

    """
    cols = ['label', 'ainvs', 'intpts']
    filename = os.path.join(base_dir, 'intpts', 'intpts.{}'.format(r))
    print("Writing data to {}".format(filename))
    n = 0
    with open(filename, 'w') as outfile:
        for label, record in data.items():
            line = make_line(record, cols)
            outfile.write(line +"\n")
            n += 1
    print("{} lines written to {}".format(n, filename))

def write_degphi(data, r, base_dir=MATSCHKE_DIR):
    r"""
    Write file base_dir/alldegphi/alldegphi.<r>

    """
    cols = ['conductor', 'isoclass', 'number', 'ainvs', 'degree']
    filename = os.path.join(base_dir, 'alldegphi', 'alldegphi.{}'.format(r))
    print("Writing data to {}".format(filename))
    n = 0
    with open(filename, 'w') as outfile:
        for label, record in data.items():
            line = make_line(record, cols)
            outfile.write(line +"\n")
            n += 1
    print("{} lines written to {}".format(n, filename))

def write_datafiles(data, r, base_dir=MATSCHKE_DIR):
    r"""Write file base_dir/<ft>/<ft>.<r> for ft in ['curvedata',
    'classdata', 'intpts', 'alldegphi']
    """
    for writer in [write_curvedata, write_classdata, write_intpts, write_degphi]:
        writer(data, r, base_dir)

def read_write_data(infilename, base_dir=MATSCHKE_DIR, verbose=1):
    print("Reading from {}".format(infilename))
    N = infilename.split(".")[-1]
    data = make_new_data(infilename, base_dir=base_dir, verbose=verbose)
    write_datafiles(data, N, base_dir)

# How to make a 2adic file: run 2adic.m on a file containing one line
# per curve, where each line has the form
#
# N c i [a1,a2,a3,a4,a6] *
#
# The magma script only read the a-invariants from what it finds
# between the first "[" and the first "]", and also copies anything
# which precedes the first "[" verbatim into the output line (one line
# per curve).
#
# Similarly, Sutherland's script for mod p Galois representations runs
# the function ComputeQGaloisImages(infilename, outfilename) in
# ~/galrep/nfgalrep.m, where the format of infilename is
#
# label:[a1,a2,a3,a4,a6]
#
# To run both on the range NN, given a file  curvedata/ curvedata.NN in directory D:
#
# ~/ecdata/scripts/make_galrep.sh NN D
#
# from a directory containing curvedata/curvedata.${r}, which will
# create (or overwrite) the files 2adic/2adic.${r} and
# galrep/galrep.${r}.


