# Sage's E.aplist(100) returns a list of the Fourier coefficients for
# p<100.  For the aplist files, we want to replace the coefficient for
# p|N with the W-eigenvalue (the root number) and append the
# W-eigenvalues for p|N, p>100.  Not relevant for making LMFDBupload
# files.

from sage.all import prime_range

def wstr(n, w):  # str(n) with enough spaces prepended to give width w
    a = str(n)
    if len(a) < w:
        a = ' ' * (w - len(a)) + a
    return a

def my_ap(E, D, p):
    if p.divides(D):
        return E.root_number(p)
    return E.ap(p)

def my_ap_str(E, D, p):
    if p.divides(D):
        a = E.root_number(p)
        if a == 1:
            if p > 23:
                return '  +'
            return ' +'
        if p > 23:
            return '  -'
        return ' -'
    if p > 23:
        return wstr(E.ap(p), 3)
    return wstr(E.ap(p), 2)

def my_aplist(E):
    D = E.discriminant()
    ap = [my_ap_str(E, D, p) for p in prime_range(100)]
    qlist = D.support()
    for q in qlist:
        if q > 100:
            if E.root_number(q) == 1:
                ap.append('+('+str(q)+')')
            else:
                ap.append('-('+str(q)+')')
    return ' '.join(ap)
