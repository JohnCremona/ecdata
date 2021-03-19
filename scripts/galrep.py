# Sage interface to Sutherland's Magma script for Galois images

GALREP_SCRIPT_DIR = "/home/jec/galrep"

def init_galrep(mag, script_dir=GALREP_SCRIPT_DIR):
    """
    Load the 2adic magma script into this magma process
    """
    mag.eval('cwd:=GetCurrentDirectory();')
    mag.eval('ChangeDirectory("{}");'.format(script_dir))
    mag.eval('load "nfgalrep.m";')
    mag.eval('ChangeDirectory(cwd);')

def get_galrep_data(E, mag):
    """
    Use Magma script to compute mod-p Galois image data

    E is an elliptic curve over Q

    mag is a magma process.

    NB before calling this, the caller must have called init_galrep()

    Returns a (possibly empty) list of string representing image codes
    """
    return str(mag.ComputeQGaloisImage(E)).split()


