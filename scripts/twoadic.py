# Sage interface to 2adic Magma script

from codec import parse_twoadic_string
import os

TWOADIC_SCRIPT_DIR = "/home/jec/ecdata/scripts"

def init_2adic(mag, script_dir=TWOADIC_SCRIPT_DIR):
    """
    Load the 2adic magma script into this magma process
    """
    script = os.path.join(script_dir, "2adic.m")
    mag.eval('load "{}";'.format(script))

def get_2adic_data(E, mag):
    """
    Use 2adic.m Magma script to compute the 2-adic image data

    E is an elliptic curve over Q

    mag is a magma process.

    NB before calling this, the caller must have called init_2adic()
    """
    if E.has_cm():
        s = "inf inf [] CM"
    else:
        s = str(mag.make_2adic_string(E))
    return parse_twoadic_string(s)


