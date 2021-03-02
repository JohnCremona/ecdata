# Manage child Magma processes:
#
# User gets a Magma instance using get_magma(), which restarts after
# magma_count uses (default 100).  This avoid the possible problems
# with using Magma on hundreds of thousands of curves, without the
# overhead of starting a new one every time.
#
# Also, the scripts for 2adic and mod p galois images are read in to
# the Magam instance.
#

from sage.all import Magma
from twoadic import init_2adic
from galrep import init_galrep

magma = Magma()
init_2adic(magma)
init_galrep(magma)

magma_count = 0
magma_count_max = 100 # restart magma after this number of uses
MagmaEffort = 100000   # 1000 not enough for  282203479a1 (rank 2)

def get_magma(verbose=False):
    global magma, magma_count, magma_count_max
    if magma_count == magma_count_max:
        if verbose:
            print("Restarting Magma")
        magma.quit()
        magma = Magma()
        init_2adic(magma)
        init_galrep(magma)
        magma_count = 1
    else:
        if verbose:
            print("Reusing Magma (count={})".format(magma_count))
        magma_count += 1
    return magma

