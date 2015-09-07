# given Elist, a list of elliptic curves over Q which is
#  "conductor-complete", i.e. contains every curve of conductor N for
#  every N which occurs.

# 1. Create a list of triples [N, reversed-aplist, ainvs] which when
#  sorted (by default, i.e. lexicographically) put the curves into
#  lmfdb order

def make_sorted_Edata(Elist):
    Edata = [[E.conductor(), tuple(E.aplist(100)), E.ainvs()]
             for E in Elist]
    return(sorted(Edata))

# 2. Create from the above a (sorted) list of triples [label,anvs]
#  with first entry the lmfdb label

from sage.databases.cremona import cremona_letter_code as code

def make_labels(Edata):
    data = {}
    for E in Edata:
        N, ap, ai = E
        if not N in data:
           data[N] = {}
        if not ap in data[N]:
           data[N][ap] = {}
        data[N][ap][ai] = "?"
    for N in sorted(data.keys()):
        ncl = len(data[N])
        ncu = sum([len(data[N][C]) for C in data[N]])
        if N<300000: # check counts against database
           assert ncl==len(list(cremona_optimal_curves([N])))
           assert ncu==len(list(cremona_curves([N])))
        print "Conductor %s=%s: %s curves in %s isogeny classes" % (N,N.factor(),ncu,ncl)
    	for iso, ap in enumerate(sorted(data[N].keys())):
	    for cu, ai in enumerate(sorted(data[N][ap].keys())):
	        data[N][ap][ai] = ".".join([str(N),code(iso)+str(cu+1)])
    return data

# 3. Make a table from previous

def make_table_from_data(Edata, plist, filename=None):
    if filename:
       f = open(filename,'w')
    else:
       f = sys.stdout

    f.write("label\t[a1,a2,a3,a4,a6]\tord_p(N) for p in %s\n"%plist)
    for N in sorted(Edata.keys()):
        vals = " ".join([str(N.valuation(p)) for p in plist])
        Ndata=Edata[N]
        for isodata in sorted(Ndata.keys()):
             curdata=Ndata[isodata]
             for C in sorted(curdata.keys()):
                 f.write("%s\t%s\t%s\n" % (curdata[C],list(C),vals))
    if filename:
        f.close()

# 4. Put it all together

def make_table(Elist,plist,filename=None):
    make_table_from_data(make_labels(make_sorted_Edata(Elist)),plist,filename)

#########################################################
# processing a Stein-Watkins data file to make labels etc
#########################################################

def class_size_from_width(w):
    if w==1: return 1
    if w in [12,16]: return 8
    if w in [18,8]: return 6
    if w in [4,6,10,14,15,21,27]: return 4
    if w in [9,25]: return 3
    return 2

def liststr(l):
    return str(list(l)).replace(' ','')

def torsion_structure(tx):
    if 'x' in tx:
        return [2,2*ZZ(tx.replace('x',''))]
    return [ZZ(tx)]

def write_file(ofile,data,N):
    for iso, ap in enumerate(sorted(data[N].keys())):
        for cu, ai in enumerate(sorted(data[N][ap].keys())):
            C = data[N][ap][ai]
            lab = [str(N),code(iso),str(cu+1)]
            C['label'] = ".".join([str(N),code(iso)+str(cu+1)])
            aistr = liststr(ai)
            ofile.write(" ".join(lab+[aistr,str(C['rank']),str(C['ntorsion'])]))
            ofile.write("\n")


def read_sw(filename, verbose=False):
    f = open(filename)
    o = open(filename+".allcurves",'w')
    data = {}
    N = 0
    for L in f.readlines():
        if L[0]=='[':  # then it's a curve line
           ai, ordD, Sha, ntors = L.split()
           ai = tuple(eval(ai))
           E = EllipticCurve(ai)
           gtors = torsion_structure(ntors)
           ntors = prod(gtors,1)
           assert E.conductor() == N
           if verbose: print E.ainvs()
           aplist = tuple(E.aplist(100))
           if not aplist in data[N]:
               data[N][aplist] = {}
           data[N][aplist][ai] = {}
           data[N][aplist][ai]['label'] = '?'
           data[N][aplist][ai]['rank'] = r
           data[N][aplist][ai]['ntorsion'] = ntors
        else: # it's a (new) isogeny class line
           newN, Nfact, r, Lr1, class_degree, degphi = L.split()
           newN = ZZ(newN) # conductor
           if newN!=N:
               if N:
                   write_file(o,data,N)
               N=newN
           if not N in data:
               data[N] = {}
           r = ZZ(r) # analytic rank
           Lr1 = float(Lr1) # leading L-series coeff
           class_degree = int(class_degree) # width of isogeny class
           class_size = class_size_from_width(class_degree)
           if verbose: print "Isogeny class N=%s, r=%s, size=%s" % (N,r,class_size)
    f.close()
    write_file(o,data,N)
    o.close()
    return data
