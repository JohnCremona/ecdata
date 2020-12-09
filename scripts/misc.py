import os
from sage.all import ZZ, QQ, Integer, EllipticCurve, class_to_int

try:
    from sage.databases.cremona import cmp_code
except:
    pass

from files import read_data, MATSCHKE_DIR
from codec import parse_int_list, point_to_weighted_proj
from ecdb import write_curvedata, get_modular_degree, liststr, shortstr, shortstrlist, matstr, pointstr, mat_to_list_list

# one-off function to fix curvedata encoding of gens and torsion_generators
# from e.g. [(-5:625:1),(580/9:1250/27:1)]

def fix_gens(r, base_dir=MATSCHKE_DIR):
    data = read_data(base_dir, ['curvedata'], [r], True)

    def parse_points(s, E):
        if s=='[]':
            return []
        else:
            return [E([QQ(c) for c in pt.split(":")]) for pt in s[2:-2].split("),(")]

    for label, record in data.items():
        E = EllipticCurve(parse_int_list(record['ainvs']))
        for t in ['gens', 'torsion_generators']:
            record[t] = [point_to_weighted_proj(gen) for gen in parse_points(record[t], E)]

    write_curvedata(data, r+'.new', base_dir=base_dir)

# one-off function to compute modular degrees when we originally
# forgot to include them in the main script.

def make_degrees(infilename, base_dir, verbose=1):
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

    nc=0
    for label, record in alldata.items():
        nc += 1
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
            cl = E.isogeny_class(order=tuple(allcurves))
            record['isogeny_matrix'] = mat = mat_to_list_list(cl.matrix())
            record['class_deg'] = max(max(r) for r in mat)
            record['isogeny_degrees'] = mat[0]
            record['degree'] = degphi = get_modular_degree(E, label)
            record['degphilist'] = degphilist = [degphi*mat[0][j] for j in range(ncurves)]
            if verbose: print("degphilist = {}".format(degphilist))
        else:
            record1 = alldata[iso+'1']
            E = allcurves[number-1]
            record['class_deg'] = record1['class_deg']
            record['isogeny_degrees'] = record1['isogeny_matrix'][number-1]
            record['degree'] = record1['degphilist'][number-1]

        if verbose:
            print("Modular degree = {}".format(record['degree']))
        if nc%1000==0:
            print("{} curves processed...".format(nc))
    print("Finished processing {} curves".format(nc))
    return alldata

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

def compare(Ncc1,Ncc2):
    d = Integer(Ncc1[0])-Integer(Ncc2[0])
    if d!=0:
        return d
    code1 = Ncc1[1]+Ncc1[2]
    code2 = Ncc2[1]+Ncc2[2]
    d = cmp_code(code1,code2)
    return d

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

