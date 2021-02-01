import os
HOME = os.getenv("HOME")
UPLOAD_DIR = os.path.join(HOME, "ecq-upload")
import sys
sys.path.append(os.path.join(HOME, 'lmfdb'))
from lmfdb import db

all_tables = [db.ec_curvedata, db.ec_localdata, db.ec_mwbsd,
          db.ec_classdata, db.ec_2adic, db.ec_galrep,
          db.ec_torsion_growth, db.ec_iwasawa]

main_tables = [db.ec_curvedata, db.ec_localdata, db.ec_mwbsd,
          db.ec_classdata, db.ec_2adic, db.ec_galrep]

all_ranges = ["{}0000-{}9999".format(n,n) for n in range(50)]
iwasawa_ranges = all_ranges[:15]
growth_ranges = all_ranges[:40]

# This one updates existing rows with revised data.  Any new rows are
# ignored!

# NB This can *only* be used for tables where the label uniquely identifies the row to be changed!  These are:

# db.ec_curvedata, db.ec_mwbsd, db.ec_classdata, db.ec_2adic, db.ec_iwasawa

# but NOT: db.ec_localdata, db.ec_galrep, db.ec_torsion_growth

# To update the latter it is necessary to delete the old rows and then use copy_from() instead of update_from_file().


def update_range(r, tables=all_tables, base_dir=UPLOAD_DIR):
    for t in tables:
        basefile = ".".join([t.search_table,r])
        file = os.path.join(UPLOAD_DIR, basefile)
        print("updating {} from {}".format(t.search_table,file))
        t.update_from_file(file)

# Use this one to add rows.  NB They must be new rows, else we end up
# with duplicate rows.  There should *not* be an 'id' column.  So a
# script to update everything from scratch, deleting all the content
# first, would go like this for all tables in all_tables or a subset.
#
# # Delete the old data in each table:
#
# for t in tables:
#     t.delete({})
#
# # Read the data files in ranges of 10000, write out the upload data
# # files, and upload the data in these to each table:
# 
# for r in all_ranges:
#     data = read_data(ranges=[r])
#     for t in tables:
#         make_table_upload_file(data, t.search_table, rows=r, include_id=False)
#     add_data(r, tables=tables)
#
# # OR: read the data files in all ranges, write out the upload data
# # files for the whole range, and upload the data in these to each
# # table:
# 
# data = read_data(ranges=all_ranges)
# for t in tables:
#     make_table_upload_file(data, t.search_table, rows='all', include_id=False)
# add_data('all', tables=tables)
#

def add_data(r, tables=all_tables, base_dir=UPLOAD_DIR):
    for t in tables:
        basefile = ".".join([t.search_table,r])
        file = os.path.join(UPLOAD_DIR, basefile)
        print("updating {} from {}".format(t.search_table,file))
        t.copy_from(file)
