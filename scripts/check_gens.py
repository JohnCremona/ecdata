######################################################################
#
# Functions to check Minkowksi-reduction of generators, and compare
#
######################################################################

import os
from sage.all import EllipticCurve
from codec import split, parse_int_list, proj_to_point, point_to_proj
from red_gens import reduce_gens, check_minkowski

HOME = os.getenv("HOME")
ECDATA_DIR = os.path.join(HOME, "ecdata")
BASE_DIR = ECDATA_DIR

def check_mink(filename, base_dir=BASE_DIR):
    infilename = os.path.join(base_dir, filename)
    with open(infilename) as infile:
        n = 0
        nbad = 0
        for line in infile:
            n += 1
            data = split(line)
            label = "".join(data[:3])
            ainvs = parse_int_list(data[3])
            E = EllipticCurve(ainvs)
            rank = int(data[4])
            gens = [proj_to_point(gen, E) for gen in data[6:6 + rank]]
            if not check_minkowski(gens):
                nbad += 1
                print("{}: gens {} are not Minkowski-reduced".format(label, gens))
                print("heights: {}".format([P.height() for P in gens]))
    print("Out of {} curves, {} were reduced and {} not".format(n, n-nbad, nbad))

def rewrite_gens(filename, base_dir=BASE_DIR):
    oldfile = os.path.join(base_dir, filename)
    newfile = ".".join([oldfile, 'new'])
    with open(oldfile) as infile, open(newfile, 'w') as outfile:
        n = 0
        for line in infile:
            n += 1
            data = split(line)
            label = "".join(data[:3])
            ainvs = parse_int_list(data[3])
            E = EllipticCurve(ainvs)
            rank = int(data[4])
            gens = [proj_to_point(gen, E) for gen in data[6:6 + rank]]
            tgens = [proj_to_point(gen, E) for gen in data[6 + rank:]]
            if len(tgens) == 2 and tgens[0].order() > tgens[1].order():
                print("torsion gens for {} in wrong order".format(label))
            newgens, tgens = reduce_gens(gens, tgens)
            newgens = [point_to_proj(P) for P in newgens]
            tgens = [point_to_proj(P) for P in tgens]
            data[6:6+rank] = newgens
            data[6+rank:] = tgens
            outfile.write(" ".join(data) + "\n")
            if n%1000 == 0:
                print("{} lines output...".format(n))
    print("Done. {} lines output".format(n))

def compare_gens(filename1, filename2, base_dir=BASE_DIR):
    file1 = os.path.join(base_dir, filename1)
    file2 = os.path.join(base_dir, filename2)
    Egens = {}
    with open(file1) as infile1, open(file2) as infile2:
        for infile, infilename in zip([infile1, infile2], [file1, file2]):
            for line in infile:
                data = split(line)
                label = "".join(data[:3])
                ainvs = parse_int_list(data[3])
                E = EllipticCurve(ainvs)
                rank = int(data[4])
                gens = [proj_to_point(gen, E) for gen in data[6:6 + rank]]
                if label not in Egens:
                    Egens[label] = {}
                Egens[label][infilename] = gens
            print("Finished reading {}".format(infilename))
    nsame = 0
    ndiff = 0
    nsame_uptosign = 0
    nsame_uptotorsion = 0
    ndiff_really = 0
    with open(file1 + '.uptosign.txt', 'w') as logfile2, open(file1 + '.uptotorsion.txt', 'w') as logfile3:
        for label in Egens:
            bothgens = Egens[label]
            gens1 = bothgens[file1]
            gens2 = bothgens[file2]
            if gens1 == gens2:
                nsame += 1
            else:
                ndiff += 1
                if all([(P == Q) or (P == -Q) for P, Q in zip(gens1, gens2)]):
                    nsame_uptosign += 1
                    logfile2.write("{} gens1: {}\n{} gens2: {}\n".format(label, gens1, label, gens2))
                else:
                    if all([(P - Q).has_finite_order() or (P + Q).has_finite_order() for P, Q in zip(gens1, gens2)]):
                        nsame_uptotorsion += 1
                        logfile3.write("{} gens1: {}\n{} gens2: {}\n".format(label, gens1, label, gens2))
                    else:
                        ndiff_really += 1
    print("gens equal in {} cases, different in {} cases".format(nsame, ndiff))
    print("gens equal up to sign in {} more cases".format(nsame_uptosign))
    print("gens equal up to torsion in {} more cases".format(nsame_uptotorsion))
    if ndiff_really:
        print("gens different, even up to sign and torsion in {} cases".format(ndiff_really))
    else:
        print("In all cases the gens are the same up to sign and torsion")
