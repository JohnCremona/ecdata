# Find integral points in a fail-safe way uing both Sage and Magma,
# comparing, returning the union in all cases and outputting a warning
# message if they disagree.

from sage.all import Set
from magma import get_magma

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

