import sys
import os
from files import HOME

sys.path.append(os.path.join(HOME, 'lmfdb'))
from lmfdb import db

db.create_table(name = 'ec_curvedata',
                search_columns = {
                    'text': ['label', 'lmfdb_label', 'iso', 'lmfdb_iso'],
                    'numeric': ['regulator'],
                    'numeric[]': ['ainvs', 'jinv', 'min_quad_twist_ainvs'],
                    'smallint': ['iso_nlabel', 'number', 'lmfdb_number', 'cm', 'num_bad_primes',
                                 'optimality', 'manin_constant', 'torsion', 'rank',
                                 'analytic_rank', 'signD', 'class_deg', 'class_size',
                                 'min_quad_twist_disc', 'faltings_index', 'faltings_ratio'],
                    'smallint[]':  ['isogeny_degrees', 'nonmax_primes', 'torsion_structure', 'torsion_primes', 'sha_primes'],
                    'integer': ['conductor', 'nonmax_rad', 'num_int_pts', 'sha'],
                    'integer[]': ['bad_primes'],
                    'bigint': ['degree'],
                    'boolean': ['semistable', 'potential_good_reduction'],
                },
                label_col='label',
                sort=['conductor', 'iso_nlabel', 'number'],
                id_ordered=True
)

db.create_table(name = 'ec_localdata',
                search_columns = {
                    'text': ['label', 'lmfdb_label'],
                    'smallint': ['tamagawa_number', 'kodaira_symbol', 'reduction_type', 'root_number', 'conductor_valuation', 'discriminant_valuation', 'j_denominator_valuation'],
                    'integer': ['prime']
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
                    'text': ['label', 'lmfdb_label', 'twoadic_label'],
                    'smallint': ['twoadic_index', 'twoadic_log_level'],
                    'smallint[]': ['twoadic_gens'],
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

#######  Indexes ##############

db.ec_curvedata.create_index(['isogeny_degrees'], type='gin')
db.ec_curvedata.create_index(['nonmax_primes'], type='gin')
db.ec_curvedata.create_index(['ainvs'], type='btree')
db.ec_curvedata.create_index(['cm'], type='btree')
db.ec_curvedata.create_index(['conductor', 'iso_nlabel', 'lmfdb_number'], type='btree')
db.ec_curvedata.create_index(['iso'], type='btree')
db.ec_curvedata.create_index(['jinv', 'id'], type='btree')
db.ec_curvedata.create_index(['label'], type='btree')
db.ec_curvedata.create_index(['label', 'number'], type='btree')
db.ec_curvedata.create_index(['lmfdb_label'], type='btree')
db.ec_curvedata.create_index(['lmfdb_label', 'lmfdb_number'], type='btree')
db.ec_curvedata.create_index(['lmfdb_label', 'number'], type='btree')
db.ec_curvedata.create_index(['lmfdb_iso'], type='btree')
db.ec_curvedata.create_index(['lmfdb_number'], type='btree')
db.ec_curvedata.create_index(['number'], type='btree')
db.ec_curvedata.create_index(['rank'], type='btree')
db.ec_curvedata.create_index(['rank', 'number'], type='btree')
db.ec_curvedata.create_index(['sha', 'id'], type='btree')
db.ec_curvedata.create_index(['sha', 'rank', 'id'], type='btree')
db.ec_curvedata.create_index(['sha', 'rank', 'torsion', 'id'], type='btree')
db.ec_curvedata.create_index(['torsion'], type='btree')
db.ec_curvedata.create_index(['torsion_structure'], type='btree')
db.ec_curvedata.create_index(['nonmax_rad', 'id'], type='btree')
db.ec_curvedata.create_index(['id'], type='btree')
db.ec_curvedata.create_index(['semistable'], type='btree')
db.ec_curvedata.create_index(['semistable', 'conductor', 'iso_nlabel', 'lmfdb_number'], type='btree')
db.ec_curvedata.create_index(['potential_good_reduction'], type='btree')
db.ec_curvedata.create_index(['potential_good_reduction', 'conductor', 'iso_nlabel', 'lmfdb_number'], type='btree')
db.ec_curvedata.create_index(['class_size'], type='btree')
db.ec_curvedata.create_index(['class_deg'], type='btree')

db.ec_classdata.create_index(['iso'], type='btree')
db.ec_classdata.create_index(['lmfdb_iso'], type='btree')

db.ec_localdata.create_index(['label'], type='btree')
db.ec_localdata.create_index(['lmfdb_label'], type='btree')

db.ec_mwbsd.create_index(['label'], type='btree')
db.ec_mwbsd.create_index(['lmfdb_label'], type='btree')

db.ec_2adic.create_index(['label'], type='btree')
db.ec_2adic.create_index(['lmfdb_label'], type='btree')

db.ec_galrep.create_index(['label'], type='btree')
db.ec_galrep.create_index(['lmfdb_label'], type='btree')

db.ec_torsion_growth.create_index(['label'], type='btree')
db.ec_torsion_growth.create_index(['lmfdb_label'], type='btree')

db.ec_iwasawa.create_index(['label'], type='btree')
db.ec_iwasawa.create_index(['lmfdb_label'], type='btree')
