import os
import sys
sys.path.insert(0, '/home/jec/ecdata/scripts/')

from sage.all import (EllipticCurve, Integer, ZZ, Set, factorial,
                      mwrank_get_precision, mwrank_set_precision, srange, prod, copy, gcd)
from magma import get_magma

from red_gens import reduce_tgens, reduce_gens
from trace_hash import TraceHashClass

from files import (parse_line_label_cols, parse_curvedata_line,
                   parse_allgens_line_simple, parse_extra_gens_line, write_datafiles,
                   make_paricurves)

from codec import (parse_int_list, parse_int_list_list, point_to_weighted_proj,
                   weighted_proj_to_affine_point, split_galois_image_code,
                   curve_from_inv_string, shortstr, liststr, matstr, shortstrlist,
                   point_to_proj, mat_to_list_list)

from twoadic import get_2adic_data
from galrep import get_galrep_data
from intpts import get_integral_points
from aplist import my_aplist
from moddeg import get_modular_degree

from ec_utils import (mwrank_saturation_precision,
                      mwrank_saturation_maxprime, get_gens,
                      map_points)


from sage.databases.cremona import parse_cremona_label, class_to_int, cremona_letter_code

HOME = os.getenv("HOME")
sys.path.append(os.path.join(HOME, 'ecdata', 'scripts'))

FR_values = srange(1, 20) + [21, 25, 37, 43, 67, 163] # possible values of Faltings ratio

modular_degree_bound = 1000000 # do not compute modular degree if conductor greater than this

# Given a filename like curves.000-999, read the data in the file,
# compute the isogeny class for each curve, and output (1)
# allcurves.000-999, (2) allisog.000-999 (with the same suffix).  Also
# compute the bsd data for each curve and output (3) allbsd.000-999,
# (4) allgens.000-999 (with the same suffix), (5) degphi.000-999, (6)
# intpts.000-999, (7) alldegphi.000-999

# Version to compute gens & torsion gens too

def make_datafiles(infilename, mode='w', verbose=False, prefix="t"):
    infile = open(infilename)
    _, suf = infilename.split(".")
    allcurvefile = open(prefix+"allcurves."+suf, mode=mode)
    allisogfile = open(prefix+"allisog."+suf, mode=mode)
    allbsdfile = open(prefix+"allbsd."+suf, mode=mode)
    allgensfile = open(prefix+"allgens."+suf, mode=mode)
    degphifile = open(prefix+"degphi."+suf, mode=mode)
    alldegphifile = open(prefix+"alldegphi."+suf, mode=mode)
    apfile = open(prefix+"aplist."+suf, mode=mode)
    intptsfile = open(prefix+"intpts."+suf, mode=mode)
    for L in infile.readlines():
        if verbose:
            print("="*72)
        N, cl, _, ainvs, r, _, _ = L.split()
        E = EllipticCurve(parse_int_list(ainvs))
        r = int(r)
        # Compute the isogeny class
        Cl = E.isogeny_class()
        Elist = Cl.curves
        mat = Cl.matrix()
        maps = Cl.isogenies()
        ncurves = len(Elist)
        print("class {} (rank {}) has {} curve(s)".format(N+cl, r, ncurves))
        line = ' '.join([str(N), cl, str(1), ainvs, shortstrlist(Elist), matstr(mat)])
        allisogfile.write(line+'\n')
        if verbose:
            print("allisogfile: {}".format(line))
        # compute BSD data for each curve
        torgroups = [F.torsion_subgroup() for F in Elist]
        torlist = [G.order() for G in torgroups]
        torstruct = [list(G.invariants()) for G in torgroups]
        torgens = [[P.element() for P in G.gens()] for G in torgroups]
        cplist = [F.tamagawa_product() for F in Elist]
        omlist = [F.real_components()*F.period_lattice().real_period() for F in Elist]

        Lr1 = E.pari_curve().ellanalyticrank()[1].sage() / factorial(r)
        if r == 0:
            genlist = [[] for F in Elist]
            reglist = [1 for F in Elist]
        else:
            Plist = get_gens(E, r, verbose)
            genlist = map_points(maps, Plist)
            prec0 = mwrank_get_precision()
            mwrank_set_precision(mwrank_saturation_precision)
            if verbose:
                print("genlist (before saturation) = {}".format(genlist))
            genlist = [Elist[i].saturation(genlist[i], max_prime=mwrank_saturation_maxprime)[0] for i in range(ncurves)]
            if verbose:
                print("genlist (before reduction) = {}".format(genlist))
            genlist = [Elist[i].lll_reduce(genlist[i])[0] for i in range(ncurves)]
            mwrank_set_precision(prec0)
            if verbose:
                print("genlist  (after reduction)= {}".format(genlist))
            reglist = [Elist[i].regulator_of_points(genlist[i]) for i in range(ncurves)]
        shalist = [Lr1*torlist[i]**2/(cplist[i]*omlist[i]*reglist[i]) for i in range(ncurves)]
        squares = [n*n for n in srange(1, 100)]
        for i, s in enumerate(shalist):
            if round(s) not in squares:
                print("bad sha value %s for %s" % (s, str(N)+cl+str(i+1)))
                print("Lr1 = %s" % Lr1)
                print("#t  = %s" % torlist[i])
                print("cp  = %s" % cplist[i])
                print("om  = %s" % omlist[i])
                print("reg = %s" % reglist[i])
                return Elist[i]
        if verbose:
            print("shalist = {}".format(shalist))

        # compute modular degrees
        # degphilist = [e.modular_degree(algorithm='magma') for e in Elist]
        # degphi = degphilist[0]
        # NB The correctness of the following relies on E being optimal!
        try:
            degphi = E.modular_degree(algorithm='magma')
        except RuntimeError:
            degphi = ZZ(0)
        degphilist1 = [degphi*mat[0, j] for j in range(ncurves)]
        degphilist = degphilist1
        if verbose:
            print("degphilist = {}".format(degphilist))

        # compute aplist for optimal curve only
        aplist = my_aplist(E)
        if verbose:
            print("aplist = {}".format(aplist))

        # Compute integral points (x-coordinates)
        intpts = [get_integral_points(Elist[i], genlist[i]) for i in range(ncurves)]
        if verbose:
            for i, Ei, xs in zip(range(ncurves), Elist, intpts):
                print("{}{}{} = {}: intpts = {}".format(N, cl, (i+1), Ei.ainvs(), xs))

        # Output data for optimal curves

        # aplist
        line = ' '.join([N, cl, aplist])
        if verbose:
            print("aplist: {}".format(line))
        apfile.write(line+'\n')

        # degphi
        if degphi:
            line = ' '.join([N, cl, '1', str(degphi), str(Set(degphi.prime_factors())).replace(' ', ''), shortstr(E)])
        else:
            line = ' '.join([N, cl, '1', str(degphi), str(Set([]).replace(' ', '')), shortstr(E)])
        if verbose:
            print("degphifile: {}".format(line))
        degphifile.write(line+'\n')

        # Output data for each curve
        for i in range(ncurves):
            # allcurves
            line = ' '.join([N, cl, str(i+1), shortstr(Elist[i]), str(r), str(torlist[i])])
            allcurvefile.write(line+'\n')
            if verbose:
                print("allcurvefile: {}".format(line))

            # allbsd
            line = ' '.join([N, cl, str(i+1), shortstr(Elist[i]), str(r), str(torlist[i]), str(cplist[i]), str(omlist[i]), str(Lr1), str(reglist[i]), str(shalist[i])])
            allbsdfile.write(line+'\n')
            if verbose:
                print("allbsdfile: {}".format(line))

            # allgens (including torsion gens, listed last)
            line = ' '.join([str(N), cl, str(i+1), shortstr(Elist[i]), str(r)]
                            + [liststr(torstruct[i])]
                            + [point_to_proj(P) for P in genlist[i]]
                            + [point_to_proj(P) for P in torgens[i]]
                           )
            allgensfile.write(line+'\n')
            if verbose:
                print("allgensfile: {}".format(line))
            # intpts
            line = ''.join([str(N), cl, str(i+1)]) + ' ' + shortstr(Elist[i]) + ' ' + liststr(intpts[i])
            intptsfile.write(line+'\n')
            if verbose:
                print("intptsfile: {}".format(line))

            # alldegphi
            line = ' '.join([str(N), cl, str(i+1), shortstr(Elist[i]), liststr(degphilist[i])])
            alldegphifile.write(line+'\n')
            if verbose:
                print("alldegphifile: {}".format(line))

    infile.close()
    allbsdfile.close()
    allgensfile.close()
    allcurvefile.close()
    allisogfile.close()
    degphifile.close()
    intptsfile.close()
    apfile.close()


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
    _, suf = infilename.split(".")
    allmaninfile = open(prefix+"opt_man."+suf, mode=mode)
    allisogfile = open("allisog/allisog."+suf)
    last_class = ""
    manins = []
    area0 = 0    # else pyflakes objects
    degrees = [] # else pyflakes objects
    for L in infile.readlines():
        N, cl, num, ainvs, _ = L.split(' ', 4)
        this_class = N+cl
        num = int(num)
        E = EllipticCurve(parse_int_list(ainvs))
        lattice1 = None
        if this_class == last_class:
            deg = degrees[int(num)-1]
            #assert ideg == deg
            area = E.period_lattice().complex_area()
            mc = round((deg*area/area0).sqrt())
            # print("{}{}{}: ".format(N, cl, num))
            # print("degree =     {}".format(deg))
            # print("area   =     {}".format(area))
            # print("area ratio = {}".format(area/area0))
            # print("c^2 =        {}".format(deg*area/area0))
            manins.append(mc)
            if num == 3:
                lattice1 = None
            elif not (num == 2 and deg == 2 and mc == 2):
                lattice1 = None
            if lattice1:
                print("Class {} is special".format(lattice1))
        else: # start a new class
            if manins and verbose: # else we're at the start
                print("class {} has Manin constants {}".format(last_class, manins))
            isogmat = allisogfile.readline().split()[-1]
            isogmat = parse_int_list_list(isogmat)
            degrees = isogmat[0]
            if verbose:
                print("class {}".format(this_class))
                print("isogmat: {}".format(isogmat))
                print("degrees: {}".format(isogmat))
            area0 = E.period_lattice().complex_area()
            manins = [1]
            if num == 1  and len(isogmat) == 2 and isogmat[0][1] == 2 and E.discriminant() > 0:
                lattice1 = this_class
            last_class = this_class

        # construct output line for this curve
        #
        # SPECIAL CASE 990h
        #
        if N == 90 and cl == 'h':
            opt = str(int(num == 3))
            manins = [1, 1, 1, 1]
        else:
            opt = str(int(num == 1))
        mc = str(manins[-1])
        line = ' '.join([str(N), cl, str(num), shortstr(E), opt, mc])

        allmaninfile.write(line+'\n')
        if verbose:
            print("allmaninfile: {}".format(line))
        # if int(N) > 100:
        #     break
    if manins and verbose: # else we're at the start
        # This output line is the only reason we keep the list of all m.c.s
        print("class {} has Manin constants {}".format(last_class, manins))
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
            N, c, _ = parse_cremona_label(L.split()[0][:-1])
            if N == lastN:
                o.write(" {}".format(1+class_to_int(c)))
            else:
                if lastN != 0:
                    o.write(" 0\n")
                o.write("{} {}".format(N, 1+class_to_int(c)))
                lastN = N
                n += 1
    o.write(" 0\n")
    n += 1
    o.close()
    print("wrote {} lines to {}".format(n, outfile))

def check_sagedb(N1, N2, a4a6bound=100):
    """Sanity check that Sage's elliptic curve database contains all
    curves [a1, a2, a3, a4, a6] with a1, a3 in [0, 1], a2 in [-1, 0, 1], a4,a6
    in [-100..100] and conductor in the given range.

    Borrowed from Bill Allombert (who found a missing curve of
    conductor 406598 this way).
    """
    from sage.all import CremonaDatabase
    CDB = CremonaDatabase()
    def CDB_curves(N):
        return [c[0] for c in CDB.allcurves(N).values()]
    Nrange = srange(N1, N2+1)
    ncurves = 12*(2*a4a6bound+1)**2
    print("Testing {} curves".format(ncurves))
    n = 0
    for a1 in range(2):
        for a2 in range(-1, 2):
            for a3 in range(2):
                for a4 in range(-a4a6bound, a4a6bound+1):
                    for a6 in range(-a4a6bound, a4a6bound+1):
                        ai = [a1, a2, a3, a4, a6]
                        n += 1
                        if n%1000 == 0:
                            print("test #{}/{}".format(n, ncurves))
                        try:
                            E = EllipticCurve(ai).minimal_model()
                        except ArithmeticError: #singular curve
                            continue
                        N = E.conductor()
                        if N not in Nrange:
                            continue
                        if list(E.ainvs()) not in CDB_curves(N):
                            print("Missing N={}, ai={}".format(N, ai))

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
# If you add split_by_N=True to process_raw_curves() then instead of
# one allcurves file you will get one per conductor N, each of the
# form allcurves.<N>.  This can be useful for using parallel
# processing in the next step, after which you have to concatenate all
# the output files.  In this case allcurves_file will be ignored, and
# a file <curves_file>.conductors will also be written containing the
# distinct conductors in the input file.
#
#
# Step 2: compute most of the data, run also in ecdb/scripts after reading ecdb.py:
#
# sage: read_write_data(allcurves_file, CURVE_DIR)
#
# We now have files CURVE_DIR/F/F.NN for F in ['curvedata', 'classdata', 'intpts', 'alldegphi'].
#
# Step 2a: if we have processed allcurves.RANGE with more than one
#  conductor, then for each conductor N we will have files F/F.N for F
#  in ['curvedata', 'classdata', 'intpts', 'alldegphi'].  We now check
#  that these are complete and merge them into just 4 files F/F.RANGE.
#
# (o) RANGE=...
#     CURVE_DIR=...
#     cd ${CURVE_DIR}
#
# (i) extract conductors
# awk '{print $1;}' allcurves/allcurves.${RANGE} | sort -n | uniq > conductors.RANGE
#
# (ii) check files are complete (this is quick)

# for F in curvedata classdata intpts alldegphi; do echo $F; for N in `cat conductors.${RANGE}`; do if ! [ -e ${CURVE_DIR}/${F}/${F}.${N} ]; then echo ${F}"."${N}" missing"; fi; done; done;
#
# (iii) merge files (this takes ages for ~200000 conductors)
# Here we don't bother with alldegphi files as they will be empty for large N.

# for F in curvedata classdata intpts; do echo $F; for N in `cat conductors.${RANGE}`; do cat ${CURVE_DIR}/${F}/${F}.${N} >> ${CURVE_DIR}/${F}/${F}.${RANGE}; done; done
#
# Better way if you can find all the files easily with a wildcard:
# e.g. for all conductors in range 1e7-1e8:
# echo curvedata.???????? | xargs cat > curvedata.${RANGE}
#
# or
# pushd curvedata
# echo $(for N in $(cat ../conductors.${RANGE}); do echo curvedata\.${N}; done) | xargs cat > curvedata.${RANGE}
# popd
#
# Step 3: run magma scripts to compute galrep and 2adic data. From ecdata/scripts:
# ./make_galrep.sh NN CURVE_DIR
#
# This uses CURVE_DIR/allcurves/allcurves.NN and
# creates   CURVE_DIR/ft/ft.NN for ft in ['galrep', '2adic'].
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

    N class_code class_size number lmfdb_number ainvs

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

    with open(os.path.join(base_dir, infilename)) as infile:
        for L in infile:
            data = L.split()
            assert len(data) in [1, 2]
            invs = data[-1]
            E = curve_from_inv_string(invs)
            ainvs = E.ainvs()
            N = E.conductor()
            if len(data) == 2:
                N1 = ZZ(data[0])
                if N1 != N:
                    raise ValueError("curve with invariants {} has conductor {}, not {} as input".format(invs, N, N1))
            ncurves += 1
            if ncurves%1000 == 0 and verbose:
                print("{} curves read from {} so far...".format(ncurves, infilename))

            # now we have the curve E, check if we have seen it (or an isogenous curve) before
            repeat = False

            if N in allcurves:
                if any(ainvs in c for c in allcurves[N]):
                    repeat = True
                    if verbose > 1:
                        print("Skipping {} of conductor {}: repeat".format(ainvs, N))
                else:
                    if verbose > 1:
                        print("New curve {} of conductor {}".format(ainvs, N))
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
        print("Conductor {}:\t{} isogeny classes,\t{} curves".format(N, ncl, nc))

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
        nap = 100
        ok = False
        while not ok:
            aplists = dict([(cl[0], cl[0].aplist(100, python_ints=True)) for cl in CC])
            aps = list(aplists.values())
            ok = (len(list(Set(aps))) == len(aps))
            if not ok:
                print("{} ap not enough for conductor {}, increasing...".format(nap, N))
            nap += 100
        # sort the isogeny classes:
        CC.sort(key=lambda cl: aplists[cl[0]])
        for ncl, cl in enumerate(CC):
            nNcl += 1
            class_code = cremona_letter_code(ncl)
            # sort the curves in two ways (LMFDB, Faltings)
            cl.sort(key=curve_key_LMFDB)
            cl_Faltings = copy(cl)
            cl_Faltings.sort(key=curve_key_Faltings)
            class_size = len(cl)
            for nE_F, E in enumerate(cl_Faltings):
                nNcu += 1
                ainvs = shortstr(E)
                nE_L = cl.index(E)
                line = " ".join([sN, class_code, str(class_size), str(nE_F+1), str(nE_L+1), ainvs])
                #print(line)
                outfile.write(line+"\n")
        return nNcu, nNcl

    if split_by_N:
        with open(os.path.join(base_dir, infilename+".conductors"), 'w') as Nfile:
            for N in conductors:
                Nfile.write("{}\n".format(N))
                outfilename_N = "allcurves/allcurves.{}".format(N)
                with open(os.path.join(base_dir, outfilename_N), 'w') as outfile:
                    nNcu, nNcl = output_one_conductor(N, allcurves[N], outfile)
                    print("N={}: {} curves in {} classes output to {}".format(N, nNcu, nNcl, outfilename_N))

    else:
        with open(os.path.join(base_dir, outfilename), 'w') as outfile:
            for N in conductors:
                nNcu, nNcl = output_one_conductor(N, allcurves[N], outfile)
                print("N={}: {} curves in {} classes output to {}".format(N, nNcu, nNcl, outfilename))

def make_new_data(infilename, base_dir, Nmin=None, Nmax=None, PRECISION=128, verbose=1,
                  allgensfilename=None, oldcurvedatafile=None, extragensfilename=None):
    alldata = {}
    nc = 0
    labels_by_conductor = {}
    tempdir = os.path.join(base_dir, "temp")

    with open(os.path.join(base_dir, "allcurves", infilename)) as infile:
        for L in infile:
            sN, isoclass, class_size, number, lmfdb_number, ainvs = L.split()
            N = ZZ(sN)
            if Nmin and N < Nmin:
                continue
            if Nmax and N > Nmax:
                continue
            nc += 1
            iso = ''.join([sN, isoclass])
            label = ''.join([iso, number])
            lmfdb_isoclass = isoclass
            lmfdb_iso = '.'.join([sN, isoclass])
            lmfdb_label = ''.join([lmfdb_iso, lmfdb_number])
            iso_nlabel = class_to_int(isoclass)
            number = int(number)
            lmfdb_number = int(lmfdb_number)
            class_size = int(class_size)
            ainvs = parse_int_list(ainvs)
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
            if label in alldata:
                raise RuntimeError("duplicate label {} in {}".format(label, infilename))
            alldata[label] = record
            if N in labels_by_conductor:
                labels_by_conductor[N].append(label)
            else:
                labels_by_conductor[N] = [label]

    allN = list(labels_by_conductor.keys())
    allN.sort()
    firstN = min(allN)
    lastN = max(allN)
    if verbose:
        print("{} curves read from {}, with {} distinct conductors from {} to {}".format(nc, infilename, len(labels_by_conductor), firstN, lastN))

    gens_dict = {}
    if allgensfilename:
        print("Reading from {}".format(allgensfilename))
        n = 0
        with open(os.path.join(base_dir, allgensfilename)) as allgensfile:
            for L in allgensfile:
                if Nmin or Nmax:
                    label, record = parse_line_label_cols(L)
                    N = record['conductor']
                    if Nmin and N < Nmin:
                        continue
                    if Nmax and N > Nmax:
                        continue
                n += 1
                label, record = parse_allgens_line_simple(L)
                gens_dict[tuple(record['ainvs'])] = record['gens']
                if n%10000 == 0:
                    print("Read {} curves from {}".format(n, allgensfilename))

    if oldcurvedatafile:
        print("Reading from {}".format(oldcurvedatafile))
        n = 0
        with open(os.path.join(base_dir, "curvedata", oldcurvedatafile)) as oldfile:
            for L in oldfile:
                label, record = parse_curvedata_line(L)
                N = int(record['conductor'])
                if Nmin and N < Nmin:
                    continue
                if Nmax and N > Nmax:
                    continue
                n += 1
                gens = [weighted_proj_to_affine_point(P) for P in record['gens']]
                gens_dict[tuple(record['ainvs'])] = gens
                if n%10000 == 0:
                    print("Read {} curves from {}".format(n, oldcurvedatafile))

    if extragensfilename:
        print("Reading from {}".format(extragensfilename))
        n = 0
        with open(os.path.join(base_dir, "allgens", extragensfilename)) as genfile:
            for L in genfile:
                N = ZZ(L.split()[0])
                if Nmin and N < Nmin:
                    continue
                if Nmax and N > Nmax:
                    continue
                N, ainvs, gens = parse_extra_gens_line(L)
                n += 1
                if gens:
                    gens_dict[ainvs] = gens
                if n%10000 == 0:
                    print("Read {} curves from {}".format(n, extragensfilename))

    N0 = 0
    for N, labels in labels_by_conductor.items():
        if N != N0 and N0:
            # extract from alldata the labels for conductor N0, and output it
            if verbose:
                print("Finished conductor {}".format(N0))
                print("Writing data files for conductor {}".format(N0))
            write_datafiles(dict([(lab, alldata[lab]) for lab in labels_by_conductor[N0]]),
                            N0, tempdir)
        N0 = N
        if verbose:
            print("Processing conductor {}".format(N))
        for label in labels:
            record = alldata[label]
            if verbose:
                print("Processing curve {}".format(label))
            iso = record['iso']
            number = record['number']
            first = (number == 1) # tags first curve in each isogeny class
            ncurves = record['class_size']
            if first:
                alllabels = [iso+str(k+1) for k in range(ncurves)]
                allcurves = [EllipticCurve(alldata[lab]['ainvs']) for lab in alllabels]
                E = allcurves[0]
            else:
                record1 = alldata[iso+'1']
                E = allcurves[number-1]
            assert N == E.conductor()
            if verbose > 1:
                print("E = {}".format(E.ainvs()))

            record['jinv'] = Ej = E.j_invariant()
            record['potential_good_reduction'] = (Ej.denominator() == 1)
            D = E.discriminant()
            record['absD'] = ZZ(D).abs()
            record['signD'] = int(D.sign())
            record['cm'] = int(E.cm_discriminant()) if E.has_cm() else 0

            if first:
                record['aplist'] = E.aplist(100, python_ints=True)
                record['anlist'] = E.anlist(20, python_ints=True)
                record['trace_hash'] = TraceHashClass(record['iso'], E)
            else:
                record['aplist'] = record1['aplist']
                record['anlist'] = record1['anlist']
                record['trace_hash'] = record1['trace_hash']

            if verbose > 1:
                print("aplist done")
            local_data = [{'p': int(ld.prime().gen()),
                           'ord_cond':int(ld.conductor_valuation()),
                           'ord_disc':int(ld.discriminant_valuation()),
                           'ord_den_j':int(max(0, -(Ej.valuation(ld.prime().gen())))),
                           'red':int(ld.bad_reduction_type()),
                           'rootno':int(E.root_number(ld.prime().gen())),
                           'kod':ld.kodaira_symbol()._pari_code(),
                           'cp':int(ld.tamagawa_number())}
                          for ld in E.local_data()]

            record['tamagawa_numbers'] = cps = [ld['cp'] for ld in local_data]
            record['kodaira_symbols'] = [ld['kod'] for ld in local_data]
            record['reduction_types'] = [ld['red'] for ld in local_data]
            record['root_numbers'] = [ld['rootno'] for ld in local_data]
            record['conductor_valuations'] = cv = [ld['ord_cond'] for ld in local_data]
            record['discriminant_valuations'] = [ld['ord_disc'] for ld in local_data]
            record['j_denominator_valuations'] = [ld['ord_den_j'] for ld in local_data]

            record['semistable'] = all([v == 1 for v in cv])
            record['tamagawa_product'] = tamprod = prod(cps)
            if verbose > 1:
                print("local data done")

            twoadic_data = get_2adic_data(E, get_magma())
            record.update(twoadic_data)
            if verbose > 1:
                print("2-adic data done")

            record['modp_images'] = image_codes = get_galrep_data(E, get_magma())
            record['nonmax_primes'] = pr = [int(split_galois_image_code(s)[0]) for s in image_codes]
            record['nonmax_rad'] = prod(pr)
            if verbose > 1:
                print("galrep data done")

            T = E.torsion_subgroup()
            tgens = [P.element() for P in T.gens()]
            tgens.sort(key=lambda P: P.order())
            tgens = reduce_tgens(tgens)
            tor_struct = [P.order() for P in tgens]
            record['torsion_generators'] = [point_to_weighted_proj(gen) for gen in tgens]
            record['torsion_structure'] = tor_struct
            record['torsion'] = torsion = prod(tor_struct)
            record['torsion_primes'] = [int(p) for p in Integer(torsion).support()]
            if verbose > 1:
                print("torsion done")

            if first: # else add later to avoid recomputing a.r.

                # Analytic rank and special L-value
                ar, sv = E.pari_curve().ellanalyticrank(precision=PRECISION)
                record['analytic_rank'] = ar = ar.sage()
                record['special_value'] = sv = sv.sage()/factorial(ar)

                # compute isogenies so we can map points from #1 to the rest:
                cl = E.isogeny_class(order=tuple(allcurves))
                record['isogeny_matrix'] = mat = mat_to_list_list(cl.matrix())
                record['class_deg'] = max(max(r) for r in mat)
                record['isogeny_degrees'] = mat[0]
                isogenies = cl.isogenies() # for mapping points later
                if N <= modular_degree_bound:
                    record['degree'] = degphi = get_modular_degree(E, label)
                    if verbose > 1:
                        print("degphi = {}".format(degphi))
                else:
                    record['degree'] = degphi = 0
                    if verbose > 1:
                        print("degphi not computed as conductor > {}".format(modular_degree_bound))
                record['degphilist'] = [degphi*mat[0][j] for j in range(ncurves)]
            else:

                record['analytic_rank'] = ar = record1['analytic_rank']
                record['special_value'] = sv = record1['special_value']
                record['class_deg'] = record1['class_deg']
                record['isogeny_degrees'] = record1['isogeny_matrix'][number-1]
                record['degree'] = record1['degphilist'][number-1]

            if verbose > 1:
                print("analytic rank done: {}".format(ar))

            gens_missing = False

            if ar == 0:
                record['gens'] = gens = []
                record['regulator'] = reg = 1
                record['ngens'] = 0
                record['heights'] = []
                record['rank'] = 0
                record['rank_bounds'] = [0, 0]

            else: # positive rank
                if first:
                    if verbose > 1:
                        print("{}: an.rk.={}, finding generators".format(label, ar))
                    #if 'gens' in record:
                    ainvs = tuple(record['ainvs'])
                    if ainvs in gens_dict:
                        gens = [E(P) for P in gens_dict[ainvs]]
                        if verbose > 1:
                            print("..already have gens {}, just saturating (p<{})...".format(gens, mwrank_saturation_maxprime))
                        gens, n, _ = E.saturation(gens, max_prime=mwrank_saturation_maxprime)
                        ngens = len(gens)
                        if ngens < ar:
                            if verbose:
                                print("Warning: a.r.={} but we only have {} gens".format(ar, ngens))
                        if verbose > 1:
                            if n and n > 1:
                                print("..saturation index was {}, new gens: {}".format(n, gens))
                            else:
                                print("..saturated already")
                    else:
                        gens = get_gens(E, ar, verbose) # this returns saturated points
                    ngens = len(gens)
                    if ngens < ar:
                        gens_missing = True
                        print("{}: analytic rank = {} but we only found {} generators".format(label, ar, ngens))
                    else:
                        if verbose > 1:
                            print("...done, generators {}".format(gens))
                    record['rank_bounds'] = [ngens, ar]
                    record['rank'] = None if gens_missing else ngens
                    # so the other curves in the class know their gens:
                    record['allgens'] = map_points(isogenies, gens, verbose) # returns saturated sets of gens
                else:
                    gens = record1['allgens'][number-1]
                    record['rank'] = record1['rank']
                    record['rank_bounds'] = record1['rank_bounds']
                    gens_missing = (record['rank'] is None)

                gens, tgens = reduce_gens(gens, tgens, False, label)
                record['gens'] = [point_to_weighted_proj(gen) for gen in gens]
                record['heights'] = heights = [P.height(precision=PRECISION) for P in gens]
                reg = None if gens_missing else E.regulator_of_points(gens, PRECISION)
                record['ngens'] = len(gens)
                record['regulator'] = reg
                if verbose > 1:
                    print("... finished reduction, gens are now {}, heights {}, reg={}".format(gens, heights, reg))

            L = E.period_lattice()
            record['real_period'] = om = L.omega(prec=PRECISION) # includes #R-components factor
            record['area'] = A = L.complex_area(prec=PRECISION)
            F_ht = -A.log()/2
            R = om.parent()
            record['faltings_height'] = F_ht
            record['stable_faltings_height'] = F_ht - R(gcd(D, E.c4()**3)).log()/12

            # Analytic Sha
            if gens_missing:
                if verbose:
                    print("Unable to compute analytic Sha since #gens < analytic rank")
                record['sha_an'] = None
                record['sha'] = None
            else:
                sha_an = sv*torsion**2 / (tamprod*reg*om)
                sha = sha_an.round()
                warn = "sha_an = {}, rounds to {}".format(sha_an, sha)
                assert sha > 0, warn
                assert sha.is_square(), warn
                assert ((sha-sha_an).abs() < 1e-10), warn
                record['sha_an'] = sha_an
                record['sha'] = int(sha)

            # Faltings ratio
            FR = 1 if first else (record1['area']/A).round()
            assert FR in FR_values, "F ratio = {}/{} = {}".format(record1['area'], A, FR)
            record['faltings_ratio'] = FR

            if verbose > 1:
                print(" -- getting integral points...")
            record['xcoord_integral_points'] = get_integral_points(E, gens)
            if verbose > 1:
                print(" ...done: {}".format(record['xcoord_integral_points']))


            Etw, Dtw = E.minimal_quadratic_twist()
            if Etw.conductor() == N:
                record['min_quad_twist_ainvs'] = record['ainvs']
                record['min_quad_twist_disc'] = 1
            else:
                record['min_quad_twist_ainvs'] = [int(a) for a in Etw.ainvs()]
                record['min_quad_twist_disc'] = int(Dtw)

            if verbose:
                print("Finished processing {}".format(label))
                if verbose > 1:
                    print("data for {}:\n{}".format(label, record))

    # don't forget to output last conductor
    if verbose:
        print("Finished conductor {}".format(N0))
        print("Writing data files for conductor {}".format(N0))
    write_datafiles(dict([(lab, alldata[lab]) for lab in labels_by_conductor[N0]]),
                    N0, tempdir)
    outfilename = infilename.replace("allcurves.", "")
    if Nmin or Nmax:
        outfilename = "{}.{}-{}".format(outfilename, firstN, lastN)
    if verbose:
        print("Finished all, writing data files for {}".format(outfilename))

    write_datafiles(alldata, outfilename, base_dir)

    return alldata

def read_write_data(infilename, base_dir, verbose=1):
    print("Reading from {}".format(infilename))
    N = ".".join(infilename.split(".")[1:])
    data = make_new_data(infilename, base_dir=base_dir, verbose=verbose)
    write_datafiles(data, N, base_dir)

################################################################################
