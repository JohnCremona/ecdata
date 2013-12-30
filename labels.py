# Function to create alllabels file mapping Cremona labels to LMFDB labels
#
# Input: curves file, e.g. curves.230000-239999
#
# Output: alllabels file, e.g. alllabels.230000-239999
#
# Format: conductor iso number conductor lmfdb_iso lmfdb_number
#
# NB we do this by computing the isogeny class and (re)sorting it in
# each case.  This is ONLY intended for ranges where the classes
# themselves are in the correct order.
#
def make_alllabels(infilename, mode='w', pref='t', verbose=False):
    infile = file(infilename)
    pre, suf = infilename.split(".")
    alllabelsfile = file(pref+"alllabels."+suf, mode=mode)
    count=0
    for L in file(infilename).readlines():
        count +=1
        if count%1000==0:
            print L
        N, cl, num, ainvs, r, tor, d = L.split()
        E = EllipticCurve(eval(ainvs))
        curves = E.isogeny_class(use_tuple=False, order="sage")
        reordered_curves = sorted(curves, key = lambda E: E.a_invariants())
        lab1 = range(1,len(curves)+1)
        lab2 = [1+reordered_curves.index(EE) for EE in curves]
        for j in range(len(curves)):
            line = ' '.join([N,cl,str(lab1[j]),N,cl,str(lab2[j])])
            alllabelsfile.write(line+'\n')
            if verbose:
                print "alllabelsfile:  ",line
    infile.close()
    alllabelsfile.close()
