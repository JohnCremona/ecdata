from sage.databases.cremona import parse_cremona_label, cmp_code

mwrank_saturation_precision = 300
mwrank_saturation_maxprime = 200000
GP = '/home/jec/bin/gp'

def print_data(outfile, code, ainvs, r, t):
    print "Code = ",code
    print "Curve = ",EllipticCurve(ainvs)
    print "rank = ",r
    print "torsion = ",t

def put_allcurves_line(outfile, N, cl, num, ainvs, r, t):
    line = ' '.join([str(N),cl,str(num),str(ainvs).replace(' ',''),str(r),str(t)])
#    print line
#    print N, cl, num, ainvs, r, t
    outfile.write(line+'\n')

def make_allcurves_lines(outfile, code, ainvs, r, t):
    E = EllipticCurve(ainvs)
    N, cl, n = parse_cremona_label(code)
    for i, F in enumerate(E.isogeny_class().curves):
        put_allcurves_line(outfile,N,cl,str(i+1),list(F.ainvs()),r,F.torsion_order())
    outfile.flush()

def process_curve_file(infilename, outfilename, use):
    infile = file(infilename)
    outfile = file(outfilename, mode='a')
    for L in file(infilename).readlines():
        #print L
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
    #print P
    x,y,z = [Integer(c) for c in P[1:-1].split(":")]
    return [x/z,y/z]


def matstr(m):
    return str(list(m)).replace('(','[').replace(')',']').replace(' ','')

# Assuming that E is known to have rank 1, returns a point on E
# computed by Magma's HeegnerPoint command
def magma_rank1_gen(E):
    mP = magma(E).HeegnerPoint(nvals=2)[1]
#    print "Magma computes P = ", mP
    P = E([mP[i].sage() for i in [1,2,3]])
    return P

# Assuming that E is known to have rank 1, returns a point on E
# computed by GP's ellheegner() command
def pari_rank1_gen(E, stacksize=1024000000):
    from os import system, getpid
    f = 'tempfile-'+str(getpid())
    comm = "LD_LIBRARY_PATH=/usr/local/lib; echo `echo 'ellheegner(ellinit("+str(list(E.ainvs()))+"))' | %s -q -f -s %s` > %s;" % (GP,stacksize,f)
    system(comm)
    P = file(f).read().partition("[")[2].partition("]")[0]
    os.unlink(f)
    #print("PARI computes P = [%s]" % P)
    P = E([QQ(c) for c in P.split(',')])
    #print(" ---> P = %s" % P)
    return P

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
    infile = file(infilename)
    pre, suf = infilename.split(".")
    allcurvefile = file("tallcurves."+suf, mode=mode)
    allisogfile = file("tallisog."+suf, mode=mode)
    for L in file(infilename).readlines():
        #print L
        N, cl, num, ainvs, r, tor, d = L.split()
        E = EllipticCurve(eval(ainvs))
        Cl = E.isogeny_class()
        Elist = Cl.curves
        mat = Cl.matrix()
        torlist = [F.torsion_order() for F in Elist]
        for i in range(len(Elist)):
            line = ' '.join([N,cl,str(i+1),shortstr(Elist[i]),r,str(torlist[i])])
            allcurvefile.write(line+'\n')
            print "allcurvefile: ",line
        line = ' '.join([str(N),cl,str(1),ainvs,shortstrlist(Elist),matstr(mat)])
        allisogfile.write(line+'\n')
        print "allisogfile:  ",line
    infile.close()
    allcurvefile.close()
    allisogfile.close()
    
# Version using David Roe's new Isogeny Class class (trac #12768)
def make_allcurves_and_allisog_new(infilename, mode='w', verbose=False):
    infile = file(infilename)
    pre, suf = infilename.split(".")
    allcurvefile = file("tallcurves."+suf, mode=mode)
    allisogfile = file("tallisog."+suf, mode=mode)
    count=0
    for L in file(infilename).readlines():
        count +=1
        if count%1000==0:
            print L
        N, cl, num, ainvs, r, tor, d = L.split()
        E = EllipticCurve(eval(ainvs))
        Cl = E.isogeny_class(order="database")
        Elist = Cl.curves
        torlist = [F.torsion_order() for F in Elist]
        for i in range(len(Elist)):
            line = ' '.join([N,cl,str(i+1),shortstr(Elist[i]),r,str(torlist[i])])
            allcurvefile.write(line+'\n')
            if verbose:
                print "allcurvefile: ",line
        mat = Cl.matrix()
        line = ' '.join([str(N),cl,str(1),ainvs,shortstrlist(Elist),matstr(mat)])
        allisogfile.write(line+'\n')
        if verbose:
            print "allisogfile:  ",line
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
    infile = file(infilename)
    pre, suf = infilename.split(".")
    allcurvefile = file(prefix+"allcurves."+suf, mode=mode)
    allisogfile = file(prefix+"allisog."+suf, mode=mode)
    allbsdfile = file(prefix+"allbsd."+suf, mode=mode)
    allgensfile = file(prefix+"allgens."+suf, mode=mode)
    degphifile = file(prefix+"degphi."+suf, mode=mode)
    alldegphifile = file(prefix+"alldegphi."+suf, mode=mode)
    apfile = file(prefix+"aplist."+suf, mode=mode)
    intptsfile = file(prefix+"intpts."+suf, mode=mode)
    for L in file(infilename).readlines():
        if verbose: print "="*72
        #print L
        N, cl, num, ainvs, r, tor, d = L.split()
        E = EllipticCurve(eval(ainvs))
        r=int(r)
        # Compute the isogeny class
        Cl = E.isogeny_class()
        Elist = Cl.curves
        mat = Cl.matrix()
        maps = Cl.isogenies()
        ncurves = len(Elist)
        print "class ",N+cl," (rank ",r, ") has ",ncurves," curve(s)"
        line = ' '.join([str(N),cl,str(1),ainvs,shortstrlist(Elist),matstr(mat)])
        allisogfile.write(line+'\n')
        if verbose: print "allisogfile:  ",line
        # compute BSD data for each curve
        torgroups = [F.torsion_subgroup() for F in Elist]
        torlist = [G.order() for G in torgroups]
        torstruct = [list(G.invariants()) for G in torgroups]
        torgens = [[P.element() for P in G.gens()] for G in torgroups]
        cplist  = [F.tamagawa_product() for F in Elist]
        omlist  = [F.real_components()*F.period_lattice().real_period() for F in Elist]

        Lr1 = E.pari_curve().ellanalyticrank()[1].sage()
        #print("computed L^(%s)(E,1) = %s" % (r,Lr1))
        # LE = E.lseries()
        # print("constructed LE")
        # LEdok = LE.dokchitser(100) # bits precision (default is 53)
        # print("constructed LEdok")
        if r==0:
            genlist = [[] for F in Elist]
            reglist = [1 for F in Elist]
            # Lr1 = LEdok(1)
            # print("computed L(E,1) = %s" % Lr1)
        else:
            # Lr1 = LEdok.derivative(1,r) / factorial(r)
            # print("computed L^(%s)(E,1) = %s" % (r,Lr1))
            # #Lr1 = E.lseries().dokchitser().derivative(1,r)/factorial(r)
            if r==1:
                Plist = E.point_search(15)
                if len(Plist)==0:
                    try:
                        #print "using Magma's HeegnerPoint to find generator"
                        #Plist = [magma_rank1_gen(E)]
                        #print "P = ",Plist[0]
                        print "using GP's ellheegner() to find generator"
                        Plist = [pari_rank1_gen(E)]
                        print "P = ",Plist[0]
                    except:  # Magma/pari bug or something
                        Plist = E.gens()
            else:
                #print "torsion order = ",torlist[0]
                if torlist[0]%2==1:
                    if N+cl=="322074i":
                        P1 = E(QQ(95209997)/361, QQ(-796563345544)/6859)
                        P2 = E(QQ(67511363092960062552491477869533612821)/167548532744324594465910917052304,
                               QQ(-546962755962107290021339666753477014846325372323086316509)/2168757247628325524167944948382918905481652710592)
                        Plist = [P1,P2]
                        print "Special case gens for "+N+cl+": ",Plist
                    else:
                        try:
                            s = E.simon_two_descent(lim3=5000)
                            #print "using simon_two_descent: ",s
                            Plist = E.gens()
                            #print "gens: ",Plist
                        except:
                            print "Simon failed, using mwrank: "
                            Plist = E.gens()
                            #print "gens: ",Plist
                else:
                    #print "using mwrank: "
                    Plist = E.gens()
                    #print "gens: ",Plist
 #           print "gens: ",Plist
            genlist = map_points(maps,Plist)
            prec0=mwrank_get_precision()
            mwrank_set_precision(mwrank_saturation_precision)
#            genlist = [Elist[i].saturation(genlist[i], max_prime=mwrank_saturation_maxprime)[0] for i in range(ncurves)]
            genlist = [Elist[i].saturation(genlist[i], verbose=False)[0] for i in range(ncurves)]
            mwrank_set_precision(prec0)
            if verbose: print "genlist = ",genlist
            reglist = [Elist[i].regulator_of_points(genlist[i]) for i in range(ncurves)]
        shalist = [Lr1*torlist[i]**2/(cplist[i]*omlist[i]*reglist[i]) for i in range(ncurves)]
        if verbose: print "shalist = ",shalist

        # compute modular degrees
        # degphilist = [e.modular_degree(algorithm='magma') for e in Elist]
        # degphi = degphilist[0]
        # NB The correctness of the following relies on E being optimal!
        try:
            degphi = E.modular_degree()
        except RuntimeError:
            degphi = 0
        degphilist1 = [degphi*mat[0,j] for j in range(ncurves)]
        degphilist = degphilist1
        # if degphilist!=degphilist1:
        #     print "Problem with modular degrees for class %s: perhaps curve #1 is not optimal?"%(N+cl)
        if verbose: print "degphilist = ",degphilist

        # compute aplist for optimal curve only
        aplist = my_aplist(E)
        if verbose: print "aplist = ",aplist

        # Compute integral points (x-coordinates)
        intpts = [[P[0] for P in Elist[i].integral_points(mw_base=genlist[i])] for i in range(ncurves)]
        if verbose: print "intpts = ",intpts

        # Output data for optimal curves

        # aplist
        line = ' '.join([N,cl,aplist])
        if verbose: print "aplist: ",line
        apfile.write(line+'\n')

        # degphi
        line = ' '.join([N,cl,'1',str(degphi),str(Set(degphi.prime_factors())).replace(' ',''),shortstr(E)])
        if verbose: print "degphifile: ",line
        degphifile.write(line+'\n')

        # Output data for each curve
        for i in range(ncurves):
            # allcurves
            line = ' '.join([N,cl,str(i+1),shortstr(Elist[i]),str(r),str(torlist[i])])
            allcurvefile.write(line+'\n')
            if verbose: print "allcurvefile: ",line

            # allbsd
            line = ' '.join([N,cl,str(i+1),shortstr(Elist[i]),str(r),str(torlist[i]),str(cplist[i]),str(omlist[i]),str(Lr1),str(reglist[i]),str(shalist[i])])
            allbsdfile.write(line+'\n')
            if verbose: print "allbsdfile: ",line

            # allgens (including torsion gens, listed last)
            line = ' '.join([str(N),cl,str(i+1),shortstr(Elist[i]),str(r)]
                            + [liststr(torstruct[i])]
                            + [pointstr(P) for P in genlist[i]]
                            + [pointstr(P) for P in torgens[i]]
                            )
            allgensfile.write(line+'\n')
            if verbose: 
                print "allgensfile:  ",line
            # intpts
            line = ''.join([str(N),cl,str(i+1)]) + ' ' + shortstr(Elist[i]) + ' ' + liststr(intpts[i])
            intptsfile.write(line+'\n')
            if verbose: print "intptsfile: ",line

            # alldegphi
            line = ' '.join([str(N),cl,str(i+1),shortstr(Elist[i]),liststr(degphilist[i])])
            alldegphifile.write(line+'\n')
            if verbose: print "alldegphifile: ",line

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
    infile = file(infilename)
    pre, suf = infilename.split(".")
    allgensfile = file(prefix+"allgens0."+suf, mode=mode)
    for L in file(infilename).readlines():
        #print L
        N, cl, num, ainvs, r, tor = L.split()
        if int(r)==0:
            E = EllipticCurve(eval(ainvs))

        # compute torsion data
            T = E.torsion_subgroup()
            ntor = T.order()
            torstruct = list(T.invariants())
            torgens = [P.element() for P in T.gens()]
            gens = []
            reg = 1

        # Output data
            line = ' '.join([str(N),cl,num,ainvs,r]
                            + [liststr(torstruct)]
                            + [pointstr(P) for P in gens]
                            + [pointstr(P) for P in torgens]
                            )
            allgensfile.write(line+'\n')
            if verbose: 
                print "allgensfile:  ",line

    infile.close()
    allgensfile.close()

# Read allgens file without torsion and output allgens file with torsion
#

def add_torsion(infilename, mode='w', verbose=False, prefix="t"):
    infile = file(infilename)
    pre, suf = infilename.split(".")
    allgensfile = file(prefix+"allgens."+suf, mode=mode)
    for L in file(infilename).readlines():
        #print L
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
        if verbose: print line
        allgensfile.write(line+'\n')

# Read allgens file and for curves with non-cyclic torsion, make sure
# that the gens are in the same order as the group structure
# invariants:

def fix_torsion(infilename, mode='w', verbose=False, prefix="t"):
    infile = file(infilename)
    pre, suf = infilename.split(".")
    allgensfile = file(prefix+"allgens."+suf, mode=mode)
    for L in file(infilename).readlines():
        if verbose: 
            print "old line"
            print L
        N, cl, num, ainvs, r, gens = L.split(' ',5)
        gens = gens.split()
        tor_invs = gens[0]
        inf_gens = gens[1:int(r)+1]
        tor_gens = gens[int(r)+1:]
        if verbose: 
            print "old line rank = %s, gens=%s"%(r,gens)
            print tor_invs, inf_gens, tor_gens
        if len(tor_gens)<2:
            allgensfile.write(L)
        else:
            if verbose: 
                print "old line"
                print L
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
                print "new line"
                print line
            allgensfile.write(line+'\n')

def fix_all_torsion():
    for n in range(23):
        ns=str(n)
        filename = "allgens."+ns+"0000-"+ns+"9999"
        print filename
        fix_torsion(filename)


# Read allgens file (with torsion) and output paricurves file
#
def make_paricurves(infilename, mode='w', verbose=False, prefix="t"):
    infile = file(infilename)
    pre, suf = infilename.split(".")
    paricurvesfile = file(prefix+"paricurves."+suf, mode=mode)
    for L in file(infilename).readlines():
        #print L
        N, cl, num, ainvs, r, gens = L.split(' ',5)
        if int(r)==0:
            gens = []
        else:            
            gens = gens.split()[1:1+int(r)] # ignore torsion struct and gens
            gens = [pointPtoA(P) for P in gens]
        
        line = '[' + ', '.join(['"'+N+cl+num+'"',ainvs,str(gens).replace(' ','')]) + ']'
        #print line
        paricurvesfile.write(line+'\n')

def compare(Ncc1,Ncc2):
#    print "Comparing %s and %s"%(Ncc1,Ncc2)
    d = Integer(Ncc1[0])-Integer(Ncc2[0])
    if d!=0: 
#        print "d=%s"%d
        return d
    code1 = Ncc1[1]+Ncc1[2]
    code2 = Ncc2[1]+Ncc2[2]
#    print "code1=%s, code2=%s"%(code1,code2)
    d = cmp_code(code1,code2)
#    print "d=%s"%d
    return d

def merge_gens(infile1, infile2):
    pre, suf = infile1.split(".")
    infile1 = file(infile1)
    infile2 = file(infile2)
    allgensfile = file("mallgens."+suf, mode='w')
    L1 = infile1.readline()
    L2 = infile2.readline()
    while len(L1)>0 and len(L2)>0:
        #print "Line 1: %s"%L1
        #print "Line 2: %s"%L2
        N1,cl1,cu1,rest = L1.split(' ',3)
        N2,cl2,cu2,rest = L2.split(' ',3)
        if compare([N1,cl1,cu1],[N2,cl2,cu2])<0:
            #print "Outputting Line 1"
            allgensfile.write(L1)
            L1 = infile1.readline()
        else:
            #print "Outputting Line 2"
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
    infile = file(infilename)
    pre, suf = infilename.split(".")
    alldegphifile = file(prefix+"alldegphi."+suf, mode=mode)
    for L in file(infilename).readlines():
        if verbose: print L
        N, cl, num, ainvs, rest = L.split(' ',4)
        if verbose: print ainvs
        E = EllipticCurve(eval(ainvs))
        try:
            degphi = E.modular_degree()
        except RuntimeError:
            degphi = 0

        # alldegphi
        line = ' '.join([str(N),cl,str(num),shortstr(E),liststr(degphi)])
        alldegphifile.write(line+'\n')
        if verbose: print "alldegphifile: ",line
    infile.close()
    alldegphifile.close()

def check_degphi(infilename):
    infile = file(infilename)
    for L in file(infilename).readlines():
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
                print "%s: d=%s but d1=%s (ratio %s)"%(N+cl+str(num),d,d1,d1//d)
