# Custom function to make a latex 'equation' string from a-invariants
#
# assume that [a1,a2,a3] are one of the 12 reduced triples.

# This is the same as latex(EllipticCurve(ainvs)).replace("
# ","").replace("{3}","3").replace("{2}","2"), i.e. the only
# difference is that we have x^3 instead of x^{3} (and x^2 instead of
# x^{2} when a2!=0), and have no spaces.  I checked this on all curves
# of conductor <10000!

def latex_equation(ainvs):
    a1,a2,a3,a4,a6 = [int(a) for a in ainvs]
    return ''.join(['\(y^2',
                    '+xy' if a1 else '',
                    '+y' if a3 else '',
                    '=x^3',
                    '+x^2' if a2==1 else '-x^2' if a2==-1 else '',
                    '{:+}x'.format(a4) if abs(a4)>1 else '+x' if a4==1 else '-x' if a4==-1 else '',
                    '{:+}'.format(a6) if a6 else '',
                    '\)'])

