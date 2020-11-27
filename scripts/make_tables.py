from files import schemas, HOME

def table_search_columns():
    for table in schemas:
        S = schemas[table]
        search_columns = dict([(v, [k for k in S if S[k]==v]) for v in Set(S.values())])
        print("name = '{}',\nsearch_columns = {}".format(table,search_columns))

sys.path.append(os.path.join(HOME, 'lmfdb'))
from lmfdb import db

db.create_table(name = 'ec_curvedata',
                search_columns = {
                    'text': ['label', 'lmfdb_label', 'iso', 'lmfdb_iso'],
                    'numeric': ['regulator'],
                    'numeric[]': ['ainvs', 'jinv', 'min_quad_twist_ainvs'],
                    'smallint': ['iso_nlabel', 'number', 'lmfdb_number', 'cm', 'num_bad_primes',
                                 'optimality', 'manin_constant', 'torsion', 'rank',
                                 'analytic_rank', 'signD', 'class_deg',
                                 'min_quad_twist_disc', 'faltings_index', 'faltings_ratio'],
                    'smallint[]':  ['isogeny_degrees', 'nonmax_primes', 'torsion_structure', 'torsion_primes', 'sha_primes'],
                    'integer': ['conductor', 'nonmax_rad', 'num_int_pts', 'sha'],
                    'integer[]': ['bad_primes'],
                    'bigint': ['degree'],
                    'boolean': ['semistable'],
                },
                label_col='label',
                sort=['conductor', 'iso_nlabel', 'number'],
                id_ordered=True
)

db.create_table(name = 'ec_localdata',
                search_columns = {
                    'text': ['label', 'lmfdb_label'],
                    'smallint': ['prime', 'tamagawa_number', 'kodaira_symbol', 'reduction_type', 'root_number', 'conductor_valuation', 'discriminant_valuation', 'j_denominator_valuation'],
                },
                label_col='label',
                sort=['label', 'prime'],
                id_ordered=False
)

db.create_table(name = 'ec_mwbsd',
                search_columns = {
                    'text': ['label', 'lmfdb_label'],
                    'numeric': ['special_value', 'real_period', 'area', 'sha_an'],
                    'integer': ['tamagawa_product'],
                    'smallint': ['ngens'],
                    'smallint[]': ['rank_bounds'],
                    'numeric[]': ['torsion_generators', 'xcoord_integral_points', 'gens', 'heights'],
                },
                label_col='label',
                sort=['label'],
                id_ordered=False
)

db.create_table(name = 'ec_classdata',
                search_columns = {
                    'text': ['iso', 'lmfdb_iso'],
                    'bigint': ['trace_hash'],
                    'smallint': ['class_size', 'class_deg'],
                    'smallint[]': ['isogeny_matrix', 'aplist', 'anlist'],
                },
                label_col='iso',
                sort=['iso'],
                id_ordered=False
)

db.create_table(name = 'ec_2adic',
                search_columns = {
                    'text': ['label', 'lmfdb_label', '2adic_label'],
                    'smallint': ['2adic_index', '2adic_log_level'],
                    'smallint[]': ['2adic_gens'],
                },
                label_col='label',
                sort=['label'],
                id_ordered=False
)

db.create_table(name = 'ec_galrep',
                search_columns = {
                    'text': ['label', 'lmfdb_label', 'image'],
                    'smallint': ['prime'],
                },
                label_col='label',
                sort=['label', 'prime'],
                id_ordered=False
)

db.create_table(name = 'ec_torsion_growth',
                search_columns = {
                    'text': ['label', 'lmfdb_label'],
                    'smallint': ['degree'],
                    'numeric[]': ['field'],
                    'smallint[]': ['torsion'],
                },
                label_col='label',
                sort=['label', 'degree'],
                id_ordered=False
)

db.create_table(name = 'ec_iwasawa',
                search_columns = {
                    'text': ['label', 'lmfdb_label'],
                    'smallint': ['iwp0'],
                    'jsonb': ['iwdata'],
                },
                label_col='label',
                sort=['label'],
                id_ordered=False
)
