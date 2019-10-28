from sage.all import EllipticCurve, Integer, QQ, Set, magma, prime_range, factorial, mwrank_get_precision, mwrank_set_precision, srange, pari

from sage.databases.cremona import parse_cremona_label, class_to_int
try:
    from sage.databases.cremona import cmp_code
except:
    pass

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

# Assuming that E is known to have rank 1, returns a point on E
# computed by Magma's HeegnerPoint command
def magma_rank1_gen(E):
    mP = magma(E).HeegnerPoint(nvals=2)[1]
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
    magma.eval('E:=EllipticCurve({});'.format(list(E.ainvs())))
    magma.eval('pts:=[];')
    for P in gens:
        magma.eval('Append(~pts,E!{});'.format(list(P)))
    res = magma.eval('IntegralPoints(E : FBasis:=pts);')
    return [p[0] for p in eval(res.split("\n")[0].replace(":",","))]

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

# Sage's E.aplist(100) returns a list of the Fourier coefficients for
# p<100.  We want to replace the coefficient for p|N with the
# W-eigenvalue (the root number) and append the W-eigenvalues for p|N,
# p>100.

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
# allcurves.000-999, (2) allisog.000-999 (with the same suffix).

def make_allcurves_and_allisog(infilename, mode='w'):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    allcurvefile = open("tallcurves."+suf, mode=mode)
    allisogfile = open("tallisog."+suf, mode=mode)
    for L in infile.readlines():
        N, cl, num, ainvs, r, tor, d = L.split()
        E = EllipticCurve(eval(ainvs))
        Cl = E.isogeny_class()
        Elist = Cl.curves
        mat = Cl.matrix()
        torlist = [F.torsion_order() for F in Elist]
        for i in range(len(Elist)):
            line = ' '.join([N,cl,str(i+1),shortstr(Elist[i]),r,str(torlist[i])])
            allcurvefile.write(line+'\n')
            print("allcurvefile: {}".format(line))
        line = ' '.join([str(N),cl,str(1),ainvs,shortstrlist(Elist),matstr(mat)])
        allisogfile.write(line+'\n')
        print("allisogfile:  {}".format(line))
    infile.close()
    allcurvefile.close()
    allisogfile.close()
    
# Version using David Roe's new Isogeny Class class (trac #12768)
def make_allcurves_and_allisog_new(infilename, mode='w', verbose=False):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    allcurvefile = open("tallcurves."+suf, mode=mode)
    allisogfile = open("tallisog."+suf, mode=mode)
    count=0
    for L in infile.readlines():
        count +=1
        if count%1000==0:
            print(L)
        N, cl, num, ainvs, r, tor, d = L.split()
        E = EllipticCurve(eval(ainvs))
        Cl = E.isogeny_class(order="database")
        Elist = Cl.curves
        torlist = [F.torsion_order() for F in Elist]
        for i in range(len(Elist)):
            line = ' '.join([N,cl,str(i+1),shortstr(Elist[i]),r,str(torlist[i])])
            allcurvefile.write(line+'\n')
            if verbose:
                print("allcurvefile: {}".format(line))
        mat = Cl.matrix()
        line = ' '.join([str(N),cl,str(1),ainvs,shortstrlist(Elist),matstr(mat)])
        allisogfile.write(line+'\n')
        if verbose:
            print("allisogfile: {}".format(line))
    infile.close()
    allcurvefile.close()
    allisogfile.close()
    
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
        # LE = E.lseries()
        # LEdok = LE.dokchitser(100) # bits precision (default is 53)
        if r==0:
            genlist = [[] for F in Elist]
            reglist = [1 for F in Elist]
            # Lr1 = LEdok(1)
        else:
            # Lr1 = LEdok.derivative(1,r) / factorial(r)
            # #Lr1 = E.lseries().dokchitser().derivative(1,r)/factorial(r)
            if r==1:
                Plist = E.point_search(15)
                if len(Plist)==0:
                    try:
                        #Plist = [magma_rank1_gen(E)]
                        print("using GP's ellheegner() to find generator")
                        Plist = [pari_rank1_gen(E)]
                        print("P = {}".format(Plist[0]))
                    except:  # Magma/pari bug or something
                        Plist = E.gens()
            else:
                if torlist[0]%2==1:
                    if N+cl=="322074i":
                        P1 = E(QQ(95209997)/361, QQ(-796563345544)/6859)
                        P2 = E(QQ(67511363092960062552491477869533612821)/167548532744324594465910917052304,
                               QQ(-546962755962107290021339666753477014846325372323086316509)/2168757247628325524167944948382918905481652710592)
                        Plist = [P1,P2]
                        print("Special case gens for {}{}: {}".format(N,cl,Plist))
                    else:
                        try:
                            s = E.simon_two_descent(lim3=5000)
                            Plist = E.gens()
                        except:
                            print("Simon failed, using mwrank: ")
                            Plist = E.gens()
                else:
                    Plist = E.gens()
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

#
# Compute torsion gens only from allcurves file
#
def make_rank0_torsion(infilename, mode='w', verbose=False, prefix="t"):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    allgensfile = open(prefix+"allgens0."+suf, mode=mode)
    for L in infile.readlines():
        N, cl, num, ainvs, r, tor = L.split()
        if int(r)==0:
            E = EllipticCurve(eval(ainvs))

        # compute torsion data
            T = E.torsion_subgroup()
            torstruct = list(T.invariants())
            torgens = [P.element() for P in T.gens()]
            gens = []

        # Output data
            line = ' '.join([str(N),cl,num,ainvs,r]
                            + [liststr(torstruct)]
                            + [pointstr(P) for P in gens]
                            + [pointstr(P) for P in torgens]
                            )
            allgensfile.write(line+'\n')
            if verbose: 
                print("allgensfile: {}".format(line))

    infile.close()
    allgensfile.close()

# Read allgens file without torsion and output allgens file with torsion
#

def add_torsion(infilename, mode='w', verbose=False, prefix="t"):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    allgensfile = open(prefix+"allgens."+suf, mode=mode)
    for L in infile.readlines():
        N, cl, num, ainvs, r, gens = L.split(' ',5)
        gens = gens.split()
        E = EllipticCurve(eval(ainvs))
        T = E.torsion_subgroup()
        torstruct = list(T.invariants())
        torgens = [P.element() for P in T.smith_gens()]
        
        # allgens (including torsion gens, listed last)
        line = ' '.join([N,cl,num,ainvs,r]
                        + [liststr(torstruct)]
                        + gens #[pointstr(P) for P in gens]
                        + [pointstr(P) for P in torgens]
                        )
        if verbose: print(line)
        allgensfile.write(line+'\n')

    infile.close()
    allgensfile.close()

# Read allgens file and for curves with non-cyclic torsion, make sure
# that the gens are in the same order as the group structure
# invariants:

def fix_torsion(infilename, mode='w', verbose=False, prefix="t"):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    allgensfile = open(prefix+"allgens."+suf, mode=mode)
    for L in infile.readlines():
        if verbose: 
            print("old line")
            print(L)
        N, cl, num, ainvs, r, gens = L.split(' ',5)
        gens = gens.split()
        tor_invs = gens[0]
        inf_gens = gens[1:int(r)+1]
        tor_gens = gens[int(r)+1:]
        if verbose: 
            print("old line rank = %s, gens=%s"%(r,gens))
            print(tor_invs, inf_gens, tor_gens)
        if len(tor_gens)<2:
            allgensfile.write(L)
        else:
            if verbose: 
                print("old line")
                print(L)
            E = EllipticCurve(eval(ainvs))
            T = E.torsion_subgroup()
            tor_struct = list(T.invariants())
            tor_gens = [P.element() for P in T.smith_form_gens()]
            assert all([P.order()==n for P,n in zip(tor_gens,tor_struct)])

        # allgens (including torsion gens, listed last)
            line = ' '.join([N,cl,num,ainvs,r]
                        + [liststr(tor_struct)]
                        + inf_gens #[pointstr(P) for P in gens]
                        + [pointstr(P) for P in tor_gens]
                        )
            if verbose: 
                print("new line")
                print(line)
            allgensfile.write(line+'\n')

    infile.close()
    allgensfile.close()

def fix_all_torsion():
    for n in range(23):
        ns=str(n)
        filename = "allgens."+ns+"0000-"+ns+"9999"
        print(filename)
        fix_torsion(filename)


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

def compare(Ncc1,Ncc2):
    d = Integer(Ncc1[0])-Integer(Ncc2[0])
    if d!=0:
        return d
    code1 = Ncc1[1]+Ncc1[2]
    code2 = Ncc2[1]+Ncc2[2]
    d = cmp_code(code1,code2)
    return d

def merge_gens(infile1, infile2):
    pre, suf = infile1.split(".")
    infile1 = open(infile1)
    infile2 = open(infile2)
    allgensfile = open("mallgens."+suf, mode='w')
    L1 = infile1.readline()
    L2 = infile2.readline()
    while len(L1)>0 and len(L2)>0:
        N1,cl1,cu1,rest = L1.split(' ',3)
        N2,cl2,cu2,rest = L2.split(' ',3)
        if compare([N1,cl1,cu1],[N2,cl2,cu2])<0:
            allgensfile.write(L1)
            L1 = infile1.readline()
        else:
            allgensfile.write(L2)
            L2 = infile2.readline()
    while len(L1)>0:
        allgensfile.write(L1)
        L1 = infile1.readline()
    while len(L2)>0:
        allgensfile.write(L2)
        L2 = infile2.readline()

    infile1.close()
    infile2.close()
    allgensfile.close()

# Create alldegphi files from allcurves files:

def make_alldegphi(infilename, mode='w', verbose=False, prefix="t"):
    infile = open(infilename)
    pre, suf = infilename.split(".")
    alldegphifile = open(prefix+"alldegphi."+suf, mode=mode)
    for L in infile.readlines():
        if verbose: print(L)
        N, cl, num, ainvs, rest = L.split(' ',4)
        if verbose: print(ainvs)
        E = EllipticCurve(eval(ainvs))
        try:
            degphi = E.modular_degree()
        except RuntimeError:
            degphi = 0

        # alldegphi
        line = ' '.join([str(N),cl,str(num),shortstr(E),liststr(degphi)])
        alldegphifile.write(line+'\n')
        if verbose: print("alldegphifile: {}".format(line))
    infile.close()
    alldegphifile.close()

def check_degphi(infilename):
    infile = open(infilename)
    for L in infile.readlines():
        N, cl, num, ainvs, d = L.split()
        if int(N)==990 and cl=='h':
            continue
        d=int(d)
        if int(num)==1:
            d1=d
        else:
            if d1<d:
                pass
            else:
                print("%s: d=%s but d1=%s (ratio %s)"%(N+cl+str(num),d,d1,d1//d))
    infile.close()

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
    CDB = CremonaDatabase()
    def CDB_curves(N):
        return [c[0] for c in CDB.allcurves(N).values()]
    Nrange = srange(N1,N2+1)
    ncurves = 12*(2*a4a6bound+1)**2
    print("Testing {} curves".format(ncurves))
    n=0
    for a1 in xrange(2):
        for a2 in xrange(-1,2):
            for a3 in xrange(2):
                for a4 in xrange(-a4a6bound,a4a6bound+1):
                    for a6 in xrange(-a4a6bound,a4a6bound+1):
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
