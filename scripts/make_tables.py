from files import schemas, HOME

def table_search_columns():
    for table in schemas:
        S = schemas[table]
        search_columns = dict([(v, [k for k in S if S[k]==v]) for v in Set(S.values())])
        print("name = '{}',\nsearch_columns = {}".format(table,search_columns))

sys.path.append(os.path.join(HOME, 'lmfdb'))
from lmfdb import db

db.create_table(name = 'ec_q_curvedata',
                search_columns = {
                    'bigint[]': ['ainvs', 'jinv', 'min_quad_twist_ainvs'],
                    'text': ['label', 'lmfdb_label', 'iso', 'lmfdb_iso'],
                    'numeric': ['regulator'],
                    'smallint': ['iso_nlabel', 'number', 'lmfdb_number', 'cm', 'num_bad_primes',
                                 'optimality', 'manin_constant', 'torsion', 'rank',
                                 'analytic_rank', 'sha', 'signD', 'class_deg',
                                 'min_quad_twist_disc', 'faltings_index',
                                 'faltings_ratio'],
                    'integer[]': ['bad_primes'],
                    'integer': ['conductor', 'nonmax_rad', 'num_int_pts', 'degree'],
                    'boolean': ['semistable'],
                    'smallint[]':  ['isogeny_degrees', 'nonmax_primes', 'torsion_structure', 'torsion_primes', 'sha_primes'],
                },
                label_col='label',
                sort=['conductor', 'iso_nlabel', 'number'],
                id_ordered=True
)

db.create_table(name = 'ec_q_localdata',
                search_columns = {
                    'smallint': ['prime', 'tamagawa_number', 'kodaira_symbol', 'reduction_type', 'root_number', 'conductor_valuation', 'discriminant_valuation', 'j_denominator_valuation'],
                    'text': ['label', 'lmfdb_label'],
                },
                label_col='label',
                sort=['label', 'prime'],
                id_ordered=False
)

db.create_table(name = 'ec_q_mwbsd',
                search_columns = {
                    'integer[]': ['torsion_generators', 'xcoord_integral_points', 'gens'],
                    'text': ['label', 'lmfdb_label'],
                    'numeric': ['special_value', 'real_period', 'area', 'sha_an'],
                    'smallint': ['tamagawa_product', 'ngens'],
                    'smallint[]': ['rank_bounds'],
                    'numeric[]': ['heights'],
                },
                label_col='label',
                sort=['label'],
                id_ordered=False
)

db.create_table(name = 'ec_q_classdata',
                search_columns = {
                    'smallint': ['class_size', 'class_deg'],
                    'bigint': ['trace_hash'],
                    'smallint[]': ['isogeny_matrix', 'aplist', 'anlist'],
                    'text': ['iso', 'lmfdb_iso'],
                },
                label_col='iso',
                sort=['iso'],
                id_ordered=False
)

db.create_table(name = 'ec_q_2adic',
                search_columns = {
                    'smallint': ['2adic_index', '2adic_log_level'],
                    'smallint[]': ['2adic_gens'],
                    'text': ['label', 'lmfdb_label', '2adic_label'],
                },
                label_col='label',
                sort=['label'],
                id_ordered=False
)

db.create_table(name = 'ec_q_galrep',
                search_columns = {
                    'smallint': ['prime'],
                    'text': ['label', 'lmfdb_label', 'image'],
                },
                label_col='label',
                sort=['label', 'prime'],
                id_ordered=False
)

db.create_table(name = 'ec_q_torsion_growth',
                search_columns = {
                    'smallint': ['degree'],
                    'bigint[]': ['field'],
                    'smallint[]': ['torsion'],
                    'text': ['label', 'lmfdb_label'],
                },
                label_col='label',
                sort=['label', 'degree'],
                id_ordered=False
)

db.create_table(name = 'ec_q_iwasawa',
                search_columns = {
                    'smallint': ['iwp0'],
                    'jsonb': ['iwdata'],
                    'text': ['label', 'lmfdb_label'],
                },
                label_col='label',
                sort=['label'],
                id_ordered=False
)
