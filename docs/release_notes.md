# ecdata release notes (since 2001)

## Last major update: 2019-07-22

### 2021

- 1 March 2021: updated paricurves files after previous allgens
  update.  (These files are used to create the PARI/GP elldata
  optional package.)

### 2020

- 27 November 2020: updated allgens files with Minkowski-reduced
  generators mod torsion, of minimal naive height mod torsion. Thanks
  to R. Rathbun.

### 2019

- 26 November 2019: updated opt_man files after reaching conductor
  400000 in optimality checking.

- 28 October 2019: updated opt_man files after reaching conductor
  370000 in optimality checking.

  Added The Artistic License 2.0 to the repository.

- 23 October 2019: replaced all old torsion growth files with new
  ones, now split into conductor ranges of 10000 as with the other
  files (as well as by degree), and extended to degree 23.  The
  polynomials defining the extension fields are all polredabs-ed, and
  the new torsion structure over these fields has been independently
  checked (by me).  All the data was provided by Enrique Gonzalez
  Jimenez and Filip Najman.

- 11 September 2019: added one new isogeny class (with just 1 curve)
  406598c, relabelling the old 406598c to 406598d.  Checked that there
  are no more missing curve which either have |a4|, |a6|<=100 or which
  are in Stein-Watkins.  Thanks to Bill Allombert for spotting this
  omission.

- 2 September 2019: replaced the optimality/* files with opt_man/*
  files giving more information about optimality and the Manin
  constant, documented in a rewritten doc/manin.txt.

- 20 August 2019: switched labels of 3 pairs of curves in classes of
  two 2-isogenous curves for which the optimal curves was #2 (and now
  is #1).  Classes: 457226e, 461342l, 463834b

- 19 August 2019: switched labels of 10 pairs of curves in classes of
  two 2-isogenous curves for which the optimal curves was #2 (and now
  is #1).  Classes: 235470bb, 235746u, 258482a, 265706a, 333270bu,
  359282a, 369194a, 375410g, 377034t, 389774b

- 22 July 2019: Added all data for 410000-499999
                added missing data for 406000-407999
-  1 July 2019: Added all data for 400000-409999

### 2018
- 19 December 2017: Updated intpts for 207 curves

### 2017
- 23 November 2017: Added torsion growth data for 000000-399999

### 2016
- 17 October 2016: Added all data for 390000-399999
- 20 August 2016: Added all data for 380000-389999
-  7 February 2016: Added all data for 370000-379999

### 2015
- 31 October 2015: Added all data for 360000-369999
- 19 May 2015: Added all data for 350000-359999
- 11 February 2015: Added 2adic data (from Rouse)

### 2014
- 29 August 2014: Updated data for 240000-249999, 250000-259999 (both
   ranges 246400-99 and 252800-99 had been omitted entirely)
- 12 May 2014: Added all data for 340000-349999
- 25 March 2014: Added all data for 330000-339999
- 10 January 2014: Added all galdata for 320000-329999
- 10 February 2014: Added all data for 320000-329999
- 01 January 2014: Added all galdata for 300000-319999

### 2013
- 31 December 2013: Added all data for 310000-319999
- 30 December 2013: Added all data for 300000-309999
- 06 September 2013: Added galrep.* files (from Sutherland)
- 24 January 2013: Recomputed several allbigsha.* files and shas.html
- 14 January 2013:
  - Tidied up release notes
  - Recomputed allbsd.* files < 130000
  - added scripts to recreate shas.html, table.html

### 2012
- 22 October 2012: Added all data for 290000-299999
- 13 October 2012: Added all data for 280000-289999
- 27 September 2012: Added all data for 270000-279999
- 26 August 2012: Added all data for 260000-269999
- 20 August 2012: Added all data for 250000-259999
- July 2012 : Added all data for 240000-249999
- June 2012 : Added all data for 230000-239999
- April 2012: Added alllabels.* files and april2012 notes
  - Added alldegphi files, updated release notes
  - Switch 190096b1/2 since old b2 was optimal; updated release notes
  - For 48 isogeny classes of size 2 between 130k and 230k, swapped over curves #1 and #2
  - Added optimality.?? files, extra sections to INDEX and manin.txt
  - Remade allisog.* files for 0-12 to get correct order and filled matrices
- 27/03/12: Changed all allgens files for non-cyclic torsion groups so
 that the torsion generator orders match the group invariants.
- 02/03/12: Added all data for 220000-229999
- 30/01/12: Replaced corrupt paricurves.130000-139999
- 13/01/12: Added all data for 210000-219999

### 2011
- 01/12/11: Corrected paricurves.180000-189999 (186120bk2 was wrong)
- 21/11/11: Added all data for 200000-209999
- 26/10/11: Added all data for 190000-199999
- 15/09/11:
  - Added all data for 180000-189999
  - Corrected allbsd data for 130725bc and 135346c (thanks to R Andrew Ohana)
- 09/08/11: Added all data for 170000-179999
- 26/07/11: Added all data for 160000-169999
- 08/07/11: Added all data for 150000-159999
- 28/06/11:
  - Added all data for 140000-149999
  - corrected count.130000-139999 (N=135045 had 0 not 6)
- 21/04/11: Updated intpts.130000-139999 (data for 130416 was garbled),
  thanks to Geoff bailey.
- 15/04/11: Updated allgens.n0000-n9999 for n in 0..12 to include
   torsion generators and rank 0 curves.  Updated INDEX to
   reflect this.
- 12/04/11:
  - Corrected paricurves.130000-139999 (deleted extraneous Pythonic Ls)
  - Corrected L''(1) and |Sha| for 135014r1
- 08/04/11
  - Added all data for 130000-139999
  - Changed allgens format to include torsion (so there is one
    line in allgens* for each line in allcurves*).

### 2007
- 05/01/07: Corrected aplist file entries for bad primes over 100, for
   N divisible by more than one such prime (or by the square
   of such a prime). Files affected: aplist.* (except
   aplist.00000-09999)

### 2006
- 15/12/06: Corrected Sha(51522f3) from 4 to 9 (thanks to Tom Fisher).
   Files affected: allbigsha.[05]0000-*, allbsd.50000*, shas.html
- 24/09/06: Conductors divisible by 216 between 90000 and 99999 were
  missing, now included.  All files *.90000-99999 affected.
- 04/09/06: Saturation of given generators now proved for all curves!
- 04/09/06: Corrected generators for 26 curves (52266d1, 124132a1 and
  24 of conductor 122010) which were unsaturated.  Files
   affected are allgens.{5,12}*, paricurves.{5,12}* but NOT
   the allbsd or allbigsha files which are unchanged (they
   already had the correct regulators).
- 13/07/06: Corrected 20 entries in allisog.[0-7]* (curves incorrectly shown
   as self-isogenous), thanks to Tom Boothby
- 17/02/06: Corrected files degphi.1[012]0000-1[012]9999 (wrong codes)
- 16/02/06: Corrected files count.40000-49999 and count.60000-69999.

### 2005
- 20/11/05: Updated optimality claims to cover all N<50000
- 07/11/05: Minor corrections to {allbsd,allgens,paricurves}.120000-129999
- 05/11/05: Added all data for 120000-129999
- 01/11/05: Minor corrections to allbigsha* and shas.html
- 20/09/05: Minor corrections to files degphi.1[01]0000-1[01]9999,
   allbsd.[69]*, allgens.9*, allbigsha.9*, shas.html
- 18/09/05:
  - Added all data for 100000-119999
  - allgens files now use new letter codes
- 02/09/05: Changed all isogeny class codes to new simplified base-26 scheme:
 a,b,...,z,ba,bb,...,bz,ca,...
- 31/08/05: Added all data for 90000-99999
- 26/08/05: Added all data for 80000-89999
- 14/07/05:
  - Added all data for 70000-79999
  - Deleted trailing field " ***!!!***" from allbsd.* for
  - curves where field 11 (=analytic Sha) >1
  - Added allisog* files ("Table 6") giving isogeny matrix for each class
  - Removed curves.* and gens.*, as they are subsets of the allcurves.* and allgens.* files
- 20/06/05: Added all data for 60000-69999.
- 09/06/05:
  - Added all data for 50000-59999.
  - renamed files e.g. curves.00000-09999 etc;  no changes to
   content since no elliptic curves have conductor divisible by 1000.
- 27/05/05: Added all data for 40001-50000.
- 05/05/05: corrected generators for 14 curves in range 30001-40000
 which were not saturated, with consequent corrections to
 allgens.30001-40000, allbsd.30001-40000,
 allbigsha.30001-40000 and shas table.  (None of these were
  optimal curves so gens.* unchanged)
- 03/5/05: corrections to files allgens.30001-40000, gens.20001-30000,
  allgens.20001-30000, gens.30001-40000, allbsd.00001-10000,
  allgens.00001-10000 thanks to Geoff Bailey.  (Some duplications, and
  some wrong numbering of curves in an isogeny class).
- 26/4/05: added curve 25350III2 which had been missed
- 22/4/05:
  - added data for range 30001-40000
  - resorted files into uniform ranges of 10000
  - only include gzip-ed files
- 1/3/05:  added data for N=15299 (one curve) previously omitted
- 9/2/05:  added data for range 25001-30000 + certification of optimal
 curves extended from 8000 to 11000 (extended to 12000 on 13/02/05).
 Swapped the two curves in classes 15180,15624,15744 as the
  conditionally optimal one was second not first.

### 2004
- 21/6/04:  added data for range 20001-25000 + certification of generators
  for all curves so far
- 9/3/04:  corrections for 15810U3 -- previous generator was
  3*generator;  so |Sha|=9 not 1.

### 2003
- 4/4/03:  added data for range 17001-20000.
- 12/2/03:  added data for range 16001-17000.
- 17/1/03:  added data for range 15001-16000.

### 2002
- 25/10/02:  corrected torsion order for 4830N4 in allcurves.1-8000
  mention 2 new curves of rank 3 in INDEX file
- 8/10/02: Added data for the range 12001-15000, with the same status as
  for the range 8001-12000 detailed below (except that Mark Watkins has
  not -- yet -- covered this range). [He has now: 25/10/02.]
- 14/1/02: Curves in classes 8160R, 8585C, 11024B reordered.  For
  these levels the first curve in each class is guaranteed
  Gamma_0(N)-optimal.  For other levels up to 12000 where there is no
  guarantee, our first curve is the one which has been /conditionally/
  proved optimal by Mark Watkins.  It is now the case that for all
  levels in the range 8001-12000 the first curve in each class is
  either proved optimal, or at least agrees with Mark Watkins's
  conditional list of optimal curves.  The former holds whenever the
  degree of the modular parametrization is given in the curves file.

### 2001
- 12/12/01: for around 10% of classes in the range 8201-9000 the
  ordering has been changed.  In the new order the first curve in each
  class is more likely (though not guaranteed) to be optimal.  This
  change makes this range compatible with the range 9001-12000 in that
  the same strategy is used throughout for finding a curve from each
  newform.
- November 2001:  The tables are extended to level (conductor) 12000.
  - For levels up to 8000 fuller data exists, including determination
   of the optimal curve in each isogeny class and the degree of the
   modular parametrization.  The numbering of the curves in each class
   (with the optimal curve #1, except for class 990H) may be taken as
   canonical and standard in this range.
  - For levels 8001-12000 we give all curves in each isogeny class, but
   the first curve is not known to be optimal (unless the class only
   contains one curve, of course).  The numbering of the curves in
   each class containing more than one curve should NOT be taken as
   canonical and standard: it will change if the optimal curve turns
   out to be not the currently labelled first curve.  It is likely
   (based on past observation) that, where the first curve in a class
   is not optimal, it is related to the optimal curve by a 2-isogeny.
  - In the range 8001-12000 we also do not have the modular degree (in
   most cases).
  - For any individual level in the range 8001-12000, it is possible to
   carry out the extra computations needed to fill these gaps
   (determining the optimal curve and modular degree), but doing so
   for all curves in the range would probably take a few months with
   the current algorithms and programs.  An independent computation by
   Mark Watkins <watkins@math.psu.edu> has determined, conditionally
   on the Stevens conjecture, which curve in each class is optimal,
   and has determined its modular degree.  Our results agree for all
   levels up to 8000 and we expect them to agree in the higher range
   also.
- 17/9/01: Error corrected: added 5104B1 to files, omitted due to an
  error discovered while working on 2*5104=10208.  Further checking
  has shown that no further omissions were made at around the same
  time.

