######################################################################
#
#   Utility and coding/decoding functions
#
######################################################################

from sage.all import ZZ, QQ, RR, sage_eval
import re

whitespace = re.compile(r'\s+')

def split(line):
    return whitespace.split(line.strip())

def parse_int_list(s, delims=True):
    r"""
    Given a string like '[a1,a2,a3,a4,a6]' returns the list of integers [a1,a2,a3,a4,a6]
    """
    ss = s[1:-1] if delims else s
    return [] if ss=='' else [int(a) for a in ss.split(',')]

def parse_int_list_list(s, delims=True):
    r"""
    Given a string like '[[1,2,3],[4,5,6]]' returns the list of lists of integers [[1,2,3],[4,5,6]]
    """
    ss = s.replace(" ","")
    return [] if ss=='[]' else [parse_int_list(a, False) for a in ss[2:-2].split('],[')]

def proj_to_aff(s):
    r"""
    Converts projective coordinate string '[x:y:z]' to affine coordinate string '(x/z,y/z)'
    """
    x, y, z = [ZZ(c) for c in s[1:-1].split(":")]
    return "({},{})".format(x/z,y/z)

def proj_to_weighted_proj(s):
    r"""Converts projective coordinate string '[x:y:z]' to list [a,b,c]
    where [x,y,z]=[ac,b,c^3] and [x/z,y/z]=[a/c^2,b/c^3]
    """
    x, b, z = [ZZ(c) for c in s[1:-1].split(":")]
    c = x.gcd(z)
    a = x//c
    return [a,b,c]

def point_to_weighted_proj(P):
    r"""Converts rational point P=(x,y) to weighted projective coordinates [a,b,c]
    where x=a/c^2, y=b/c^3
    """
    x, y, _ = list(P)
    a = x.numerator()
    b = y.numerator()
    c = y.denominator() // x.denominator()
    return [a,b,c]

def point_to_proj(P):
    r"""Converts rational point P=(x,y) to projective coordinates [a,b,c]
    where x=a/c, y=b/c
    """
    x, y, _ = list(P)
    c = y.denominator()
    a = ZZ(c*x)
    b = ZZ(c*y)
    return "[" + ":".join([str(x) for x in [a,b,c]]) + "]"

def proj_to_point(s, E):
    r"""
    Converts projective coordinate string '[x:y:z]' to a point on E
    """
    return E.point([ZZ(c) for c in s[1:-1].split(":")])

def split_galois_image_code(s):
    """Each code starts with a prime (1-3 digits but we allow for more)
    followed by an image code for that prime.  This function returns
    two substrings, the prefix number and the rest.
    """
    p = re.findall(r'\d+', s)[0]
    return p, s[len(p):]

######################################################################
#
# Coding and decoding functions
#

str_type = type('abc')
bool_type = type(True)
list_type = type([1,2,3])
int_type = type(int(1))
ZZ_type = type(ZZ(1))
QQ_type = type(QQ(1))
RR_type = type(RR(1))
number_types = [int_type, ZZ_type, RR_type]

def encode(col):
    t = type(col)
    if t is str_type:
        return col
    if t is bool_type:
        return str(int(col))
    if t in number_types:
        return str(col)
    if t is QQ_type:
        return str([col.numer(),col.denom()]).replace(" ","")
    if t is list_type:
        return str(col).replace(" ","")
    print("no encoding for type {}".format(t))
    return col

str_cols = ['label', 'iso', 'isoclass', 'lmfdb_label', 'lmfdb_isoclass', 'lmfdb_iso']
int_cols = ['number', 'lmfdb_number', 'iso_nlabel', 'faltings_index',
            'faltings_ratio', 'conductor', 'cm', 'signD',
            'min_quad_twist_disc', 'rank', 'analytic_rank', 'ngens',
            'torsion', 'tamagawa_product', 'sha', 'class_size', 'class_deg']
bigint_cols = ['trace_hash']
int_list_cols = ['ainvs', 'isogeny_degrees', 'min_quad_twist_ainvs',
                 'bad_primes', 'tamagawa_numbers', 'kodaira_symbols',
                 'reduction_types', 'root_numbers', 'conductor_valuations',
                 'discriminant_valuations',
                 'j_denominator_valuations', 'rank_bounds',
                 'torsion_structure',
                 'aplist', 'anlist']
int_list_list_cols = ['isogeny_matrix', 'gens', 'torsion_generators']
bool_cols = ['semistable']
QQ_cols = ['jinv']
RR_cols = ['regulator', 'real_period', 'area', 'faltings_height', 'special_value', 'sha_an']
RR_list_cols = ['heights']

def decode(colname, data):
    if colname in str_cols:
        return data
    elif colname in bigint_cols:
        return ZZ(data)
    elif colname in int_cols:
        return int(data)
    elif colname in int_list_cols:
        return parse_int_list(data)
    elif colname in bool_cols:
        return bool(int(data))
    elif colname in int_list_list_cols:
        return parse_int_list_list(data)
    elif colname in RR_cols:
        return sage_eval(data)
    elif colname in QQ_cols:
        return QQ(tuple(parse_int_list(data)))
    elif colname in RR_list_cols:
        return [] if data=='[]' else [sage_eval(d) for d in data[1:-1].split(",")]
    #
    print("No decoder set for column {} (data = {})".format(colname, data))
    return data

