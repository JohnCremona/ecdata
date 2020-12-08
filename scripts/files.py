# Functions to read data from the files:
#
# alllabels.*, allgens.*, alldegphi.*, allisog.*,
# intpts.*, opt_man.*, 2adic.*, galrep.*
#
# and also torsion growth and Iwasawa data files. The latter used to
# be arranged differently; now they are not, but only exist in the
# ranges up to 50000.

import os
from sage.all import ZZ, QQ, EllipticCurve, Integer, prod, factorial, primes
from sage.databases.cremona import class_to_int, parse_cremona_label
from trace_hash import TraceHashClass
from codec import split, parse_int_list, proj_to_point, point_to_weighted_proj, decode, encode, split_galois_image_code
from red_gens import reduce_gens

HOME = os.getenv("HOME")

# Most data from John Cremona (https://github.com/JohnCremona/ecdata)
# which we assume cloned in the home directory
ECDATA_DIR = os.path.join(HOME, "ecdata")

UPLOAD_DIR = os.path.join(HOME, "ecq-upload")

# Iwasawa data from Rob Pollack (https://github.com/rpollack9974/Iwasawa-invariants) but reorganised here:
IWASAWA_DATA_DIR = ECDATA_DIR

DEFAULT_PRECISION=53
PRECISION=53 # precision of heights, regulator, real period and
             # special value, and hence analytic sha, when these are
             # computed and not just read from the files.

######################################################################
#
# Functions to parse single lines from each of the file types
#
# In each case the function returns a full label and a dict whose keys
# are (exactly?) the relevant table columns
#
######################################################################
#
# Parsing common label/ainvs columns:
#
# allgens, allisog, alldegphi, opt_man, 2adic: first 4 cols are N,iso,number,ainvs
# alllabels:                                   first 3 cols are N,iso,number
# intpts:                                      first 2 cols are label, ainvs
# galrep:                                      first 1 col  is  label
#
# so there are 1 or 3 label columns, and there may or may not be an ainvs column

def parse_line_label_cols(L, label_cols=3, ainvs=True, raw=False):
    r"""
    Parse the first columns of one line to extract label and/or ainvs

    If label_cols is 3, the first 3 columns are conductor, iso, number, else the first column is label.
    If ainvs is True, the next columnm is the ainvs.

    If raw is False, 'conductor' is an int and 'ainvs' a list of ints, otherwise they stay as strings.

    Cols filled: 'label', 'conductor', 'iso', 'number', and optionally 'ainvs'.
    """
    data = L.split()
    record = {}

    if label_cols == 1:
        record['label'] = label = data[0]
        N, isoclass, num = parse_cremona_label(label)
        sN = str(N)
        record['conductor'] = sN if raw else N
        record['isoclass'] = isoclass
        record['iso'] = ''.join([sN,isoclass])
        record['number'] = str(num) if raw else num
    else:
        record['conductor'] = data[0] if raw else int(data[0])
        record['isoclass']  = data[1]
        record['iso']       = ''.join(data[:2])
        record['number']    = data[2] if raw else int(data[2])
        record['label']     = ''.join(data[:3])

    if ainvs:
        record['ainvs'] = data[label_cols] if raw else parse_int_list(data[label_cols])

    return record['label'], record

######################################################################
#
# allgens parser
#
# This is the biggest since it's the only one where curves have to be
# constructed and nontrivial dependent column data computed.  After
# running this on all files and outputting new the computed data to a
# new set of files, this will no longer be needed.

def parse_allgens_line(line):
    r"""
    Parse one line from an allgens file

    Lines contain 6+t+r fields (columns)

    conductor iso number ainvs r torsion_structure <tgens> <gens>

    where:

    torsion_structure is a list of t = 0,1,2 ints
    <tgens> is t fields containing torsion generators
    <gens> is r fields containing generators mod torsion

    """
    label, record = parse_line_label_cols(line, 3, True)
    first = record['number'] # tags first curve in each isogeny class
    data = split(line)

    ainvs = record['ainvs']
    E = EllipticCurve(ainvs)
    N = E.conductor()
    assert N==record['conductor']
    record['bad_primes'] = bad_p = N.prime_factors() # will be sorted
    record['num_bad_primes'] = len(bad_p)
    record['jinv'] = E.j_invariant()
    record['signD'] = int(E.discriminant().sign())
    record['cm'] = int(E.cm_discriminant()) if E.has_cm() else 0

    if first:
        record['aplist'] = E.aplist(100,python_ints=True)
        record['anlist'] = E.anlist(20,python_ints=True)
        # passing the iso means that we'll only do the computation once per isogeny class
        record['trace_hash'] = TraceHashClass(record['iso'], E)

    local_data = [{'p': int(ld.prime().gen()),
                   'ord_cond':int(ld.conductor_valuation()),
                   'ord_disc':int(ld.discriminant_valuation()),
                   'ord_den_j':int(max(0,-(E.j_invariant().valuation(ld.prime().gen())))),
                   'red':int(ld.bad_reduction_type()),
                   'rootno':int(E.root_number(ld.prime().gen())),
                   'kod':ld.kodaira_symbol()._pari_code(),
                   'cp':int(ld.tamagawa_number())}
                  for ld in E.local_data()]

    record['tamagawa_numbers']    = cps = [ld['cp'] for ld in local_data]
    record['kodaira_symbols']           = [ld['kod'] for ld in local_data]
    record['reduction_types']           = [ld['red'] for ld in local_data]
    record['root_numbers']              = [ld['rootno'] for ld in local_data]
    record['conductor_valuations'] = cv = [ld['ord_cond'] for ld in local_data]
    record['discriminant_valuations']   = [ld['ord_disc'] for ld in local_data]
    record['j_denominator_valuations']  = [ld['ord_den_j'] for ld in local_data]

    record['semistable'] = all([v==1 for v in cv])
    record['tamagawa_product'] = tamprod = prod(cps)

    # NB in the allgens file all points are stored in projective
    # coordinates as [x:y:z].  In the database (as of 2020.11.18) we
    # store both sets of generators as lists of strings, using
    # projective coordinates '(x:y:z)' for gens of infinite order but
    # affine coordinates '(x/z,y/z)' for torsion gens.  Then in the
    # code (web_ec.py) we convert projective coords to affine anyway,
    # using code which can handle either, but has less to do if the
    # point is already affine.  So let's do the conversion here.

    record['rank'] = rank = int(data[4])
    record['rank_bounds'] = [rank,rank]
    record['ngens'] = rank
    tor_struct = parse_int_list(data[5])
    record['torsion_structure'] = tor_struct
    record['torsion'] = torsion = prod(tor_struct)
    record['torsion_primes'] = [int(p) for p in Integer(torsion).support()]

    # read and reduce generators and torsion generators

    gens = [proj_to_point(gen, E) for gen in data[6:6 + rank]]
    tgens = [proj_to_point(gen, E) for gen in data[6 + rank:]]
    gens, tgens = reduce_gens(gens, tgens, False, label)

    record['gens']               = [point_to_weighted_proj(gen) for gen in gens]
    record['torsion_generators'] = [point_to_weighted_proj(gen) for gen in tgens]
    record['heights'] = [P.height(precision=PRECISION) for P in gens]
    reg = E.regulator_of_points(gens, precision=PRECISION) if gens else 1
    record['regulator'] = reg

    L = E.period_lattice()
    record['real_period'] = om = L.omega(prec=PRECISION) # includes #R-components factor
    record['area'] = A = L.complex_area(prec=PRECISION)
    record['faltings_height'] = -A.log()/2

    if first: # else add later to avoid recomputing a.r.

        # Analytic rank and special L-value
        ar,sv = E.pari_curve().ellanalyticrank(precision=PRECISION)
        record['analytic_rank'] = ar = ar.sage()
        record['special_value'] = sv = sv.sage()/factorial(ar)

        # Analytic Sha
        sha_an = sv*torsion**2 / (tamprod*reg*om)
        sha = sha_an.round()
        assert sha>0
        assert sha.is_square()
        assert ((sha-sha_an).abs() < 1e-10)
        record['sha_an'] = sha_an
        record['sha'] = int(sha)
        record['sha_primes'] = [int(p) for p in sha.prime_divisors()]

    Etw, Dtw = E.minimal_quadratic_twist()
    if Etw.conductor()==N:
        record['min_quad_twist_ainvs'] = ainvs
        record['min_quad_twist_disc']  = 1
    else:
        record['min_quad_twist_ainvs'] = [int(a) for a in Etw.ainvs()]
        record['min_quad_twist_disc']  = int(Dtw)

    return label,  record

######################################################################
#
# alllabels parser
#

def parse_alllabels_line(line):
    r""" Parses one line from an alllabels file.  Returns the label
    and a dict containing seven fields, 'conductor', 'iso', 'number',
    'lmfdb_label', 'lmfdb_iso', 'iso_nlabel', 'lmfdb_number', being strings or ints.

    Cols filled: 'label', 'conductor', 'iso', 'number', 'lmfdb_label', 'lmfdb_iso', 'lmfdb_number', 'iso_nlabel'.

    [NO] Also populates two global dictionaries lmfdb_label_to_label and
    label_to_lmfdb_label, allowing other upload functions to look
    these up.

    Input line fields:

    conductor iso number conductor lmfdb_iso lmfdb_number

    Sample input line:

    57 c 2 57 b 1

    """
    data = split(line)
    if data[0] != data[3]:
        raise ValueError("Inconsistent conductors in alllabels file: %s" % line)

    label, record = parse_line_label_cols(line, 3, False)

    record['lmfdb_isoclass'] = data[4]
    record['lmfdb_iso'] = lmfdb_iso = ''.join([data[3], '.', data[4]])
    record['lmfdb_label'] = ''.join([lmfdb_iso, data[5]])
    record['lmfdb_number'] = int(data[5])
    record['iso_nlabel'] = class_to_int(data[4])

    return label, record

######################################################################
#
# allisog parser
#

def parse_allisog_line(line):
    r"""
    Parse one line from an allisog file

    Input line fields:

    conductor iso number ainvs all_ainvs isogeny_matrix

    Sample input line:

    11 a 1 [0,-1,1,-10,-20] [[0,-1,1,-10,-20],[0,-1,1,-7820,-263580],[0,-1,1,0,0]] [[1,5,5],[5,1,25],[5,25,1]]

    """
    label, record = parse_line_label_cols(line, 3, False)
    assert record['number'] == 1

    isomat = split(line)[5][2:-2].split("],[")
    record['isogeny_matrix'] = mat = [[int(a) for a in r.split(",")] for r in isomat]
    record['class_size'] = len(mat)
    record['class_deg'] = max(max(r) for r in mat)
    record['all_iso_degs'] = dict([[n+1,sorted(list(set(row)))] for n,row in enumerate(mat)])
    record['isogeny_degrees'] = record['all_iso_degs'][1]

    # NB Every curve in the class has the same 'isogeny_matrix',
    # 'class_size', 'class_deg', and the for the i'th curve in the
    # class (for i=1,2,3,...) its 'isogeny_degrees' column is
    # all_iso_degs[i].

    return label, record

######################################################################
#
# alldegphi parser
#

def parse_alldegphi_line(line, raw=False):
    r""" Parses one line from an alldegphi file.

    Input line fields:

    conductor iso number ainvs degree

    Sample input line:

    11 a 1 [0,-1,1,-10,-20] 1
    """
    label, record = parse_line_label_cols(line, 3, False, raw=raw)
    deg = split(line)[4]
    record['degree'] = deg if raw else int(deg)
    return label, record

######################################################################
#
# intpts parser
#

def make_y_coords(ainvs,x):
    a1, a2, a3, a4, a6 = ainvs
    f = ((x + a2) * x + a4) * x + a6
    b = (a1*x + a3)
    d = (ZZ(b*b + 4*f)).isqrt()
    y = (-b+d)//2
    return [y, -b-y] if d else [y]

def count_integral_points(ainvs, xs):
    return sum([len(make_y_coords(ainvs,x)) for x in xs])

def parse_intpts_line(line, raw=False):
    r""" Parses one line from an intpts file.

    Input line fields:

    label ainvs x-coordinates_of_integral_points

    Sample input line:

    11a1 [0,-1,1,-10,-20] [5,16]
    """
    label, record = parse_line_label_cols(line, 1, True, raw=raw)
    rxs = split(line)[2]
    xs = parse_int_list(rxs)
    ainvs = record['ainvs']
    if raw:
        ainvs = parse_int_list(ainvs)
    nip = count_integral_points(ainvs, xs)
    record['xcoord_integral_points'] = rxs if raw else xs
    record['num_int_pts'] = str(nip) if raw else nip

    return label, record

######################################################################
#
# opt_man parser
#

def parse_opt_man_line(line, raw=False):
    r"""Parses one line from an opt_man file, giving optimality and Manin
    constant data.

    Input line fields:

    N iso num ainvs opt mc

    where opt = (0 if not optimal, 1 if optimal, n>1 if one of n
    possibly optimal curves in the isogeny class), and mc = Manin
    constant *conditional* on curve #1 in the class being the optimal
    one.

    Sample input lines with comments added:

    11 a 1 [0,-1,1,-10,-20] 1 1       # optimal, mc=1
    11 a 2 [0,-1,1,-7820,-263580] 0 1 # not optimal, mc=1
    11 a 3 [0,-1,1,0,0] 0 5           # not optimal, mc=5
    499992 a 1 [0,-1,0,4481,148204] 3 1       # one of 3 possible optimal curves in class g, mc=1 for all whichever is optimal
    499992 a 2 [0,-1,0,-29964,1526004] 3 1    # one of 3 possible optimal curves in class g, mc=1 for all whichever is optimal
    499992 a 3 [0,-1,0,-446624,115024188] 3 1 # one of 3 possible optimal curves in class g, mc=1 for all whichever is optimal
    499992 a 4 [0,-1,0,-164424,-24344100] 0 1 # not optimal, mc=1

    """
    label, record = parse_line_label_cols(line, 3, False)
    opt, mc = split(line)[4:]
    record['optimality'] = opt if raw else int(opt)
    record['manin_constant'] = mc if raw else int(mc)
    return label, record

######################################################################
#
# 2adic parser
#

def parse_twoadic_line(line, raw=False):
    r""" Parses one line from a 2adic file.

    Input line fields:

    conductor iso number ainvs index level gens label

    Sample input lines:

    110005 a 2 [1,-1,1,-185793,29503856] 12 4 [[3,0,0,1],[3,2,2,3],[3,0,0,3]] X24
    27 a 1 [0,0,1,0,-7] inf inf [] CM
    """
    label, record = parse_line_label_cols(line, 3, False, raw=raw)
    data = split(line)
    assert len(data)==8
    model = data[7]
    if model == 'CM':
        record['twoadic_index'] = '0' if raw else int(0)
        record['twoadic_log_level'] = None
        record['twoadic_gens'] = None
        record['twoadic_label'] = None
        return label, record

    record['twoadic_label'] = model
    record['twoadic_index'] = data[4] if raw else int(data[4])
    log_level = ZZ(data[5]).valuation(2)
    record['twoadic_log_level'] = str(log_level) if raw else int(log_level)

    rgens = data[6]
    if raw:
        record['twoadic_gens'] = rgens
    else:
        if rgens=='[]':
            record['twoadic_gens'] = []
        else:
            gens = rgens[1:-1].replace('],[','];[').split(';')
            record['twoadic_gens'] = [[int(c) for c in g[1:-1].split(',')] for g in gens]
    return label, record

######################################################################
#
# galrep parser
#

def parse_galrep_line(line, raw = False):
    r"""Parses one line from a galrep file.

    Codes follow Sutherland's coding scheme for subgroups of GL(2,p).
    Note that these codes start with a 1 or 2 digit prime followed a
    letter in ['B','C','N','S'].

    Input line fields:

    label codes

    Sample input line:

    66c3 2B 5B.1.2

    """
    label, record = parse_line_label_cols(line, 1, False, raw=raw)
    record['modp_images'] = image_codes = split(line)[1:]
    record['nonmax_primes'] = pr = [ int(split_galois_image_code(s)[0]) for s in image_codes]
    record['nonmax_rad'] = prod(pr)
    return label, record

######################################################################
#
# iwasawa parser
#

def parse_iwasawa_line(line, debug=0, raw=False):
    r"""Parses one line from an Iwasawa data input file.

    Sample line: 11 a 1 0,-1,1,-10,-20 7 1,0 0,1,0 0,0 0,1

    Fields: label (3 fields)
            a-invariants (1 field but no brackets)
            p0
            For each bad  prime:  'a'                if additive
                                  lambda,mu          if multiplicative (or 'o?' if unknown)
            For each good prime:  lambda,mu          if ordinary (or 'o?' if unknown)
                                  lambda+,lambda-,mu if supersingular (or 's?' if unknown)

    """
    if debug: print("Parsing input line {}".format(line[:-1]))
    label, record = parse_line_label_cols(line, 3, False, raw=raw)
    badp = Integer(record['conductor']).support()
    nbadp = len(badp)

    data = split(line)
    rp0 = data[4]
    p0 = int(rp0)
    record['iwp0'] = rp0 if raw else p0
    if debug: print("p0={}".format(p0))

    iwdata = {}

    # read data for bad primes

    for p,pdat in zip(badp,data[5:5+nbadp]):
        p = str(p)
        if debug>1: print("p={}, pdat={}".format(p,pdat))
        if pdat in ['o?','a']:
            iwdata[p]=pdat
        else:
            iwdata[p]=[int(x) for x in pdat.split(",")]

    # read data for all primes

    # NB Current data has p<50: if this increases to over 1000, change the next line.

    for p,pdat in zip(primes(1000),data[5+nbadp:]):
        p = str(p)
        if debug>1: print("p={}, pdat={}".format(p,pdat))
        if pdat in ['s?','o?','a']:
            iwdata[p] = pdat
        else:
            iwdata[p] = parse_int_list(pdat, delims=False)

    record['iwdata'] = iwdata
    if debug: print("label {}, data {}".format(label,record))
    return label, record

######################################################################
#
# growth parser
#

def parse_growth_line(line, raw=False):
    r"""Parses one line from a torsion growth file.

    Sample line: 14a1 [3,6][1,1,1] [2,6][2,-1,1]

    Fields: label (single field, Cremona label)

            1 or more items of the form TF (with no space between)
            with T =[n] or [m,n] and F a list of integers of length
            d+1>=3 containing the coefficients of a monic polynomial
            of degree d defining a number field (constant coeff
            first).

    Notes: (1) in each file d is fixed and contained in the filename
    (e.g. growth2.000000-399999) but can be recovered from any line
    from the length of the coefficient lists.

    (2) The files for degree d only have lines for curves where there
    is growth in degree d, so each line has at least 2 fields in it.

    (3) The returned record (dict) has one relevant key
    'torsion_growth' which is a dict with a unique key (the degree)
    and value a list of pairs [F,T] where both F and T are lists of
    ints.  It is up to the calling function to merge these for
    different degrees.

    """
    label, record = parse_line_label_cols(line, 1, False, raw=raw)
    data = [[parse_int_list(F,delims=False),parse_int_list(T,delims=False)]
            for T,F in [s[1:-1].split("][") for s in split(line)[1:]]]
    degree = len(data[0][0])-1
    record['torsion_growth'] = {degree: data}
    return label, record

iwasawa_ranges = ["{}0000-{}9999".format(n,n) for n in range(15)]

######################################################################
#
# Function to read all growth data (or a subset)
#
# (special treatment needed because of the nonstandard filenames, in directories by degree)
#

growth_degrees = [2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,21]
#NB the 0'th range is usually '00000-09999' but for growth files it's just '0-9999'
growth_ranges = ["0-9999"] + ["{}0000-{}9999".format(k,k) for k in range(1,40)]

def read_all_growth_data(base_dir=ECDATA_DIR, degrees=growth_degrees, ranges=growth_ranges, raw=False):
    r"""Read all the data in files base_dir/growth/<d>/growth.<r> where d
    is a list of degrees and r is a range.

    Return a single dict with keys labels and values curve records
    with label keys and one 'torsion_growth' key.

    """
    all_data = {}
    for r in ranges:
        if r=='00000-09999':
            r = '0-9999'
        if not r in growth_ranges:
            continue
        for d in degrees:
            if not d in growth_degrees:
                continue
            data_filename = os.path.join(base_dir, 'growth/{}/growth{}.{}'.format(d,d,r))
            n = 0
            with open(data_filename) as data:
                for L in data:
                    label, record = parse_growth_line(L, raw=raw)
                    n+=1
                    if label in all_data:
                        all_data[label]['torsion_growth'].update(record['torsion_growth'])
                    else:
                        all_data[label] = record
    return all_data

datafile_columns = {
    'curvedata': ['label', 'isoclass', 'number', 'lmfdb_label', 'lmfdb_isoclass',
                  'lmfdb_number', 'iso_nlabel', 'faltings_index', 'faltings_ratio',
                  'conductor', 'ainvs', 'jinv', 'cm',
                  'isogeny_degrees', 'semistable', 'signD',
                  'min_quad_twist_ainvs', 'min_quad_twist_disc',
                  'bad_primes', 'tamagawa_numbers', 'kodaira_symbols',
                  'reduction_types', 'root_numbers', 'conductor_valuations',
                  'discriminant_valuations', 'j_denominator_valuations',
                  'rank', 'rank_bounds', 'analytic_rank', 'ngens', 'gens',
                  'heights', 'regulator', 'torsion', 'torsion_structure',
                  'torsion_generators', 'tamagawa_product', 'real_period',
                  'area', 'faltings_height', 'special_value', 'sha_an', 'sha'],

    'classdata': ['iso', 'lmfdb_iso', 'trace_hash',
                  'class_size', 'class_deg', 'isogeny_matrix', 'aplist', 'anlist'],
    }

def parse_curvedata_line(line, raw=False):
    """
    """
    data = split(line)
    if raw:
        record = dict([(col, data[n]) for n, col in enumerate(datafile_columns['curvedata']) ])
    else:
        record = dict([(col, decode(col, data[n])) for n, col in enumerate(datafile_columns['curvedata']) ])
    record['sha_primes'] = [int(p) for p in Integer(record['sha']).prime_divisors()]
    badp = record['bad_primes']
    record['num_bad_primes'] = str(1+badp.count(",")) if raw else len(badp)
    record['torsion_primes'] = [int(p) for p in Integer(record['torsion']).prime_divisors()]
    record['lmfdb_iso'] = ".".join([str(record['conductor']),record['lmfdb_isoclass']])
    isodegs = record['isogeny_degrees']
    record['class_size'] = str(1+isodegs.count(",")) if raw else len(isodegs)
    jinv = parse_int_list(record['jinv']) if raw else record['jinv']
    record['potential_good_reduction'] = (jinv[1]==1)

    return record['label'], record

def parse_classdata_line(line, raw=False):
    """
    """
    data = split(line)
    if raw:
        record = dict([(col, data[n]) for n, col in enumerate(datafile_columns['classdata']) ])
    else:
        record = dict([(col, decode(col, data[n])) for n, col in enumerate(datafile_columns['classdata']) ])
    return record['iso']+'1', record

######################################################################
#

parsers = {'allgens': parse_allgens_line,
           'alllabels': parse_alllabels_line,
           'allisog': parse_allisog_line,
           'alldegphi': parse_alldegphi_line,
           'intpts': parse_intpts_line,
           'opt_man': parse_opt_man_line,
           '2adic': parse_twoadic_line,
           'galrep': parse_galrep_line,
           'curvedata': parse_curvedata_line,
           'classdata': parse_classdata_line,
           'growth': parse_growth_line,
           'iwasawa': parse_iwasawa_line,
}

all_file_types = list(parsers.keys())
old_file_types = ['alllabels', 'allgens', 'allisog']
new_file_types = [ft for ft in all_file_types if not ft in old_file_types]

all_ranges = ["{}0000-{}9999".format(n,n) for n in range(50)]
iwasawa_ranges = ["{}0000-{}9999".format(n,n) for n in range(15)]

######################################################################
#
# Function to read data from ['allgens', 'alllabels', 'allisog'] in
# one or more ranges and fill in additional data required for
# curvedata and classdata files.
#
# This is used in the (one-off) function make_curvedata() which makes
# curvedata and classdata files for each range.

def read_old_data(base_dir=ECDATA_DIR, ranges=all_ranges):
    r"""Read all the data in files base_dir/<ft>/<ft>.<r> where ft is each
    of ['allgens', 'alllabels', 'allisog'] and r is a range.

    Return a single dict with keys labels and values complete
    curve records.
    """
    all_data = {}

    for r in ranges:
        for ft in ['allgens', 'alllabels', 'allisog']:
            data_filename = os.path.join(base_dir, '{}/{}.{}'.format(ft,ft,r))
            parser = parsers[ft]
            n = 0
            with open(data_filename) as data:
                for L in data:
                    label, record = parser(L)
                    first = (record['number']==1)
                    # if n>100 and first:
                    #     break
                    if label:
                        if first:
                            n += 1
                        if label in all_data:
                            all_data[label].update(record)
                        else:
                            all_data[label] = record
                    if n%1000==0 and first:
                        print("Read {} classes from {}".format(n,data_filename))
            print("Read {} lines from {}".format(n,data_filename))

    # Fill in isogeny data for all curves in each class:
    for label, record in all_data.items():
        n = record['number']
        if n>1:
            record1 = all_data[label[:-1]+'1']
            record['isogeny_degrees'] = record1['all_iso_degs'][n]
            for col in ['isogeny_matrix', 'class_size', 'class_deg']:
                record[col]  = record1[col]

    # Fill in MW & BSD data for all curves in each class:
    FRout = open("FRlists.txt", 'w')
    for label, record in all_data.items():
        n = record['number']
        if n==1:
            # We sort the curves in each class by Faltings height.
            # Note that while there is always a unique curve with
            # minimal height (whose period lattice is a sublattice of
            # all the others), other heights may appear multiple
            # times.  The possible ordered lists of ratios which have
            # repeats include: [1,2,2,2], [1,2,4,4], [1,2,4,4,8,8],
            # [1,2,4,4,8,8,16,16], [1,3,3], [1,3,3,9],
            # [1,2,3,4,4,6,12,12].  As a tie-breaker we use the LMFDB
            # ordering.
            sort_key = lambda lab: [all_data[lab]['faltings_height'],all_data[lab]['lmfdb_number']]
            class_size = record['class_size']
            if class_size==1:
                record['faltings_index'] = 0
                record['faltings_ratio'] = 1
            else:
                class_labels = [label[:-1]+str(n+1) for n in range(class_size)]
                class_labels.sort(key = sort_key)
                base_label = class_labels[0]
                base_record = all_data[base_label]
                area = base_record['area']
                base_record['faltings_index'] = 0
                base_record['faltings_ratio'] = 1
                for i, lab in enumerate(class_labels):
                    if i==0:
                        continue
                    rec = all_data[lab]
                    area_ratio = area/rec['area'] # real, should be an integer
                    rec['faltings_index'] = i
                    rec['faltings_ratio'] = area_ratio.round()
                FRout.write("{}: {}\n".format(label,[all_data[lab]['faltings_ratio'] for lab in class_labels]))
        else:
            record1 = all_data[label[:-1]+'1']
            for col in ['analytic_rank', 'special_value', 'aplist', 'anlist', 'trace_hash']:
                record[col]  = record1[col]


            # Analytic Sha
            sha_an = record['special_value']*record['torsion']**2 / (record['tamagawa_product']*record['regulator']*record['real_period'])
            sha = sha_an.round()
            assert sha>0
            assert sha.is_square()
            assert ((sha-sha_an).abs() < 1e-10)
            record['sha_an'] = sha_an
            record['sha'] = int(sha)
            record['sha_primes'] = [int(p) for p in sha.prime_divisors()]

    FRout.close()
    return all_data


######################################################################
#
# Output functions
#
######################################################################
#
# make one line

def make_line(E, columns):
    """
    Given a curve record E, return a string of the selected columns.
    """
    return ' '.join([encode(E[col]) for col in columns])

######################################################################
#
# one-off function to read {allgens, alllabels, allisog} files for one
# or more ranges, and write {curvedata, classdata} files for
# the same ranges.

def make_curvedata(base_dir=ECDATA_DIR, ranges=all_ranges, prec=DEFAULT_PRECISION):
    r"""Read all the data in files base_dir/<ft>/<ft>.<r> for ft in
    ['allgens', 'alllabels', 'allisog'] and r in ranges.

    Write files base_dir/<f>/<f>.<r> for the same r and for f in
    ['curvedata', 'classdata'].

    """
    global PRECISION
    PRECISION=prec
    for r in ranges:
        print("Reading data for range {}".format(r))
        all_data = read_old_data(base_dir=base_dir, ranges=[r])

        # write out data

        for ft in ['curvedata', 'classdata']:
            cols = datafile_columns[ft]
            filename = os.path.join(base_dir, '{}/{}.{}'.format(ft,ft,r))
            print("Writing data to {}".format(filename))
            n = 0
            with open(filename, 'w') as outfile:
                for label, record in all_data.items():
                    if ft=='curvedata' or record['number']==1:
                        line = make_line(record, cols)
                        outfile.write(line +"\n")
                        n += 1
                print("{} lines written to {}".format(n, filename))

def read_data(base_dir=ECDATA_DIR, file_types=new_file_types, ranges=all_ranges, raw=True):
    r"""Read all the data in files base_dir/<ft>/<ft>.<r> where ft is a file type
    and r is a range.

    Return a single dict with keys labels and values complete
    curve records.
    """
    all_data = {}

    for r in ranges:
        for ft in file_types:
            if ft == 'growth': # special case below
                continue
            if ft == 'iwasawa' and not r in iwasawa_ranges + ['0-999']:
                continue
            data_filename = os.path.join(base_dir, '{}/{}.{}'.format(ft,ft,r))
            print("Starting to read from {}".format(data_filename))
            parser = parsers[ft]
            n = 0
            with open(data_filename) as data:
                for L in data:
                    label, record = parser(L, raw=raw)
                    first = (ft=='classdata') or (int(record['number'])==1)
                    if label:
                        if first:
                            n += 1
                        if label in all_data:
                            all_data[label].update(record)
                        else:
                            all_data[label] = record
                    if n%10000==0 and first:
                        print("Read {} classes so far from {}".format(n,data_filename))
            print("Finished reading {} classes from {}".format(n,data_filename))

    if 'curvedata' in file_types and 'classdata' in file_types:
        print("filling in class_deg from class to curve")
        for label, record in all_data.items():
            if int(record['number']) > 1:
                record['class_deg'] = all_data[label[:-1]+'1']['class_deg']

    if 'growth' in file_types:
        print("reading growth data")
        growth_data = read_all_growth_data(ranges=ranges)
        for label, record in all_data.items():
            if label in growth_data:
                record.update(growth_data[label])

    return all_data

######################################################################
#
# Function to output files which can be uploaded to the database using copy_from() or update_from_file()
#
# NB postgresql has various integer types of different *fixed*
# bit-lengths, of which te largest if 'bigint' but even that is too
# big for a 20-digit integer, so quite a few of the columns have to
# use the 'numeric' type.  The website code will cast to integers
# where necessary.
#

schemas = { 'ec_curvedata': {'label': 'text', 'lmfdb_label': 'text', 'iso': 'text', 'lmfdb_iso': 'text',
                               'iso_nlabel': 'smallint', 'number': 'smallint', 'lmfdb_number': 'smallint',
                               'ainvs': 'numeric[]', 'jinv': 'numeric[]', 'conductor': 'integer',
                               'cm': 'smallint', 'isogeny_degrees': 'smallint[]',
                               'nonmax_primes': 'smallint[]', 'nonmax_rad': 'integer',
                               'bad_primes': 'integer[]', 'num_bad_primes': 'smallint',
                               'semistable': 'boolean', 'potential_good_reduction': 'boolean',
                               'optimality': 'smallint', 'manin_constant': 'smallint',
                               'num_int_pts': 'integer', 'torsion': 'smallint',
                               'torsion_structure': 'smallint[]', 'torsion_primes': 'smallint[]',
                               'rank': 'smallint', 'analytic_rank': 'smallint',
                               'sha': 'integer',  'sha_primes': 'smallint[]', 'regulator': 'numeric',
                               'signD': 'smallint', 'degree': 'bigint', 'class_deg': 'smallint', 'class_size': 'smallint',
                               'min_quad_twist_ainvs': 'numeric[]', 'min_quad_twist_disc': 'smallint',
                               'faltings_index': 'smallint', 'faltings_ratio': 'smallint'},

            # local data: one row per (curve, bad prime)
            'ec_localdata': {'label': 'text', 'lmfdb_label': 'text',
                               'prime': 'integer', 'tamagawa_number': 'smallint', 'kodaira_symbol': 'smallint',
                               'reduction_type': 'smallint', 'root_number': 'smallint',
                               'conductor_valuation': 'smallint', 'discriminant_valuation': 'smallint',
                               'j_denominator_valuation': 'smallint'},

            'ec_mwbsd': {'label': 'text', 'lmfdb_label': 'text',
                                'torsion_generators': 'numeric[]', 'xcoord_integral_points': 'numeric[]',
                                'special_value': 'numeric', 'real_period': 'numeric', 'area': 'numeric',
                                'tamagawa_product': 'integer', 'sha_an': 'numeric', 'rank_bounds': 'smallint[]',
                                'ngens': 'smallint', 'gens': 'numeric[]', 'heights': 'numeric[]'},

            # class data: one row per isogeny class
            'ec_classdata': {'iso': 'text', 'lmfdb_iso': 'text',
                               'trace_hash': 'bigint', 'class_size': 'smallint', 'class_deg': 'smallint',
                               'isogeny_matrix': 'smallint[]',
                               'aplist': 'smallint[]', 'anlist': 'smallint[]'},

            'ec_2adic': {'label': 'text', 'lmfdb_label': 'text',
                         'twoadic_label': 'text', 'twoadic_index': 'smallint',
                         'twoadic_log_level': 'smallint', 'twoadic_gens': 'smallint[]'},

            # galrep data: one row per (curve, non-maximal prime)
            'ec_galrep': {'label': 'text', 'lmfdb_label': 'text',
                            'prime': 'smallint', 'image': 'text'},

            # torsion growth data: one row per (curve, extension field)
            'ec_torsion_growth': {'label': 'text', 'lmfdb_label': 'text',
                                  'degree': 'smallint', 'field': 'numeric[]', 'torsion': 'smallint[]'},

            'ec_iwasawa': {'label': 'text', 'lmfdb_label': 'text',
                           'iwdata': 'jsonb',  'iwp0': 'smallint'}
}

######################################################################

def postgres_encode(col, coltype):
    """
    Encoding of column data into a string for output to an upload file.

    NB A list stored in the database as a postgres array (e.g. int[] or
    numeric[]) must appear as (e.g.) {1,2,3} not [1,2,3].
    """
    if col is None:
            return "\\N"
    if coltype == "boolean":
        return "t" if col else "f"
    if type(col) == type(QQ(1)): # to handle the j-invariant
        col = [col.numer(), col.denom()]
    col = str(col).replace(" ","")
    if coltype == 'jsonb':
        col = col.replace("'",'"')
    if '[]' in coltype:
        col = col.replace("[","{").replace("]","}")
    return col

def table_cols(table, include_id=False):
    """
    Get the list of column names for a table, sorted for consistency,
    with 'label' (or 'iso' for the classdata table) moved to the
    front, and 'id' at the very front if wanted.
    """
    if table == 'ec_galrep':
        return ['label', 'lmfdb_label', 'prime', 'image']

    if table == 'ec_torsion_growth':
        return ['label', 'lmfdb_label', 'degree', 'field', 'torsion']

    cols = sorted(list(schemas[table].keys()))

    # We want the first two columns to be 'id', 'label' or 'id', 'iso' if present
    if table == 'ec_classdata':
        cols.remove('iso')
        cols.remove('lmfdb_iso')
        cols = ['iso', 'lmfdb_iso'] + cols
    else:
        cols.remove('label')
        cols.remove('lmfdb_label')
        cols = ['label', 'lmfdb_label'] + cols
    if 'id' in cols:
        cols.remove('id')
    if include_id:
        cols = ['id'] + cols
    return cols

def data_to_string(table, cols, record):
    """
    table:  a table name: one of schemas.keys()
    cols:  list of columns to output
    record: a complete curve or class record
    """
    schema = schemas[table]
    if 'id' in cols:
        schema['id'] = 'bigint'

    return "|".join([postgres_encode(record.get(col, None), schema[col]) for col in cols])

tables1 = ['ec_curvedata', 'ec_mwbsd', 'ec_2adic', 'ec_iwasawa'] # tables with one row per curve
tables2 = ['ec_classdata']                                       # table with one row per isogeny class
tables3 = ['ec_localdata',      # one row per bad prime
           'ec_galrep',         # one row per non-maximal prime
           'ec_torsion_growth', # one row per extension degree
]

all_tables = tables1 + tables2 + tables3

def make_table_upload_file(data, table, rows=None, include_id=True):
    """This version works when there is one row per curve or one per
    class.  The other cases are passed to special versions.
    """
    if not rows:
        rows = 'all'

    if table == 'ec_localdata':
        return make_localdata_upload_file(data, rows)

    if table == 'ec_galrep':
        return make_galrep_upload_file(data, rows)

    if table == 'ec_torsion_growth':
        return make_torsion_growth_upload_file(data, rows)

    include_id = include_id and (table == 'ec_curvedata')

    filename = os.path.join(UPLOAD_DIR, ".".join([table,rows]))
    allcurves = (table != 'ec_classdata')
    with open(filename, 'w') as outfile:
        print("Writing data for table {} to file {}".format(table, filename))
        if not allcurves:
            print(" (only outputting one curve per isogeny class)")

        cols = table_cols(table, include_id)
        schema = schemas[table]
        if 'id' in cols:
            schema['id'] = 'bigint'

        # Write header lines: (1) column names; (2) column types; (3) blank

        outfile.write("|".join(cols) + "\n")
        outfile.write("|".join([schema[col] for col in cols]) + "\n\n")

        n = 1
        for label, record in data.items():
            if table=='ec_iwasawa' and not 'iwdata' in record:
                continue
            if include_id:
                record['id'] = n
            if allcurves or record['number']==1:
                outfile.write(data_to_string(table, cols, record) +"\n")
            if n%10000==0:
                print("{} lines written so far...".format(n))
            n += 1
        n -= 1
        print("{} lines written to {}".format(n, filename))

def make_localdata_upload_file(data, rows=None):
    """
    This version is for ec_localdata only.  For each curve we output
    n lines where n is the number of bad primes.
    """
    if not rows:
        rows = 'all'
    table = 'ec_localdata'
    filename = os.path.join(UPLOAD_DIR, ".".join([table,rows]))
    with open(filename, 'w') as outfile:
        print("Writing data for table {} to file {}".format(table, filename))

        cols = table_cols(table, include_id=False)
        schema = schemas[table]

        # Write header lines: (1) column names; (2) column types; (3) blank

        outfile.write("|".join(cols) + "\n")
        outfile.write("|".join([schema[col] for col in cols]) + "\n\n")

        n = 1
        for label, record in data.items():
            for i in range(int(record['num_bad_primes'])):
                prime_record = {'label': record['label'], 'lmfdb_label': record['lmfdb_label'],
                                'prime': record['bad_primes'][i],
                                'tamagawa_number': record['tamagawa_numbers'][i],
                                'kodaira_symbol': record['kodaira_symbols'][i],
                                'reduction_type': record['reduction_types'][i],
                                'root_number': record['root_numbers'][i],
                                'conductor_valuation': record['conductor_valuations'][i],
                                'discriminant_valuation': record['discriminant_valuations'][i],
                                'j_denominator_valuation': record['j_denominator_valuations'][i],
                }
                outfile.write(data_to_string(table, cols, prime_record) +"\n")
                n += 1
                if n%10000==0:
                    print("{} lines written to {} so far...".format(n, filename))
        n -= 1
        print("{} lines written to {}".format(n, filename))

def make_galrep_upload_file(data, rows=None):
    """This version is for ec_galrep only.  For each curve we output n
    lines where n is the number of nonmaximal primes, so if there are
    no non-maximal primes for a curve then there is no line output for
    that curve.

    """
    if not rows:
        rows = 'all'
    table = 'ec_galrep'
    filename = os.path.join(UPLOAD_DIR, ".".join([table,rows]))
    with open(filename, 'w') as outfile:
        print("Writing data for table {} to file {}".format(table, filename))

        cols = table_cols(table, include_id=False)
        schema = schemas[table]

        # Write header lines: (1) column names; (2) column types; (3) blank

        outfile.write("|".join(cols) + "\n")
        outfile.write("|".join([schema[col] for col in cols]) + "\n\n")

        n = 1
        for label, record in data.items():
            for p, im in zip(record['nonmax_primes'], record['modp_images']):
                prime_record = {'label': record['label'], 'lmfdb_label': record['lmfdb_label'],
                                'prime': p,
                                'image': im,
                }
                outfile.write(data_to_string(table, cols, prime_record) +"\n")
                n += 1
                if n%10000==0:
                    print("{} lines written to {} so far...".format(n, filename))
        n -= 1
        print("{} lines written to {}".format(n, filename))

def make_torsion_growth_upload_file(data, rows=None):
    """This version is for ec_torsion_growth only.  For each curve we output one
    line for each field (of degree<24 currently) in which the torsion
    grows.

    """
    if not rows:
        rows = 'all'
    table = 'ec_torsion_growth'
    filename = os.path.join(UPLOAD_DIR, ".".join([table,rows]))
    with open(filename, 'w') as outfile:
        print("Writing data for table {} to file {}".format(table, filename))

        cols = table_cols(table, include_id=False)
        schema = schemas[table]

        # Write header lines: (1) column names; (2) column types; (3) blank

        outfile.write("|".join(cols) + "\n")
        outfile.write("|".join([schema[col] for col in cols]) + "\n\n")

        n = 1
        for label, record in data.items():
            if not 'torsion_growth' in record:
                continue
            for degree, dat in record['torsion_growth'].items():
                for field, torsion in dat:
                    field_record = {'label': record['label'], 'lmfdb_label': record['lmfdb_label'],
                                    'degree': degree,
                                    'field': field,
                                    'torsion': torsion,
                    }
                    outfile.write(data_to_string(table, cols, field_record) +"\n")
                    n += 1
                    if n%10000==0:
                        print("{} lines written to {} so far...".format(n, filename))
        n -= 1
        print("{} lines written to {}".format(n, filename))

def make_all_upload_files(data, tables=all_tables, rows=None, include_id=False):
    for table in tables:
        make_table_upload_file(data, table, rows=rows, include_id=include_id)

