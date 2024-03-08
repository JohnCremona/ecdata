# Testing that the sha values in curvedata files are the same as
# Sage's E.sha().an() computes for rank 0 curves

# log of runs done:

# (0) conductors 1-10k (to test): range 00000-09999 : done

# (1) 2357-smooth curves: all2357.over500k: 13153 curves, done
# test_range('all2357.over500k')

# (2)  swdb primes<3e8: 5e5-1e6  1182 curves done
#                       1e6-1e7 13262 curves done
#                       1e7-1e8 78343 curves done
#                       1e8-2e8 67459 curves done
#                       2e8-3e8 60235 curves done
# for r in ['curvedata.p.5e5-1e6', 'curvedata.p.1e6-1e7', 'curvedata.p.1e7-1e8', 'curvedata.p.1e8-2e8', 'curvedata.p.2e8-3e8']: test_range(r, SWDB_DIR)

def test_sha(curve):
     for col in ['rank', 'sha', 'ainvs', 'conductor']:
         curve[col] = decode(col, curve[col])
     if curve['rank']>0:
         return True
     if curve['conductor'] < 500000:
         E = EllipticCurve(curve['Clabel'])
     else:
         E = EllipticCurve(curve['ainvs'])
     return E.sha().an()==curve['sha']

def test_range(range, base_dir=ECDATA_DIR):
    data = read_data(base_dir, ranges=[range], file_types=['curvedata'], raw=True)
    print(f"Number of curves read: {len(data)}")
    ncurves = 0
    for label, curve in data.items():
        if decode('rank',curve['rank'])>0:
            continue
        t = test_sha(curve)
        s  = 'OK' if t else 'WRONG****************************'
        if not t:
            print(f"{label}: {s}")
        ncurves +=1
        if ncurves%1000==0:
            print(f"{ncurves} curves processed, last one {label}")
    print(f"{ncurves} curves processed, last one {label}")

