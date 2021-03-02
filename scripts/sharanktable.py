# script used to create shas.html from allbsd.* files:

from sage.all import isqrt

# we do not call the output file "shas.html" so we can compare the new
# version with the old
HTML_FILENAME = "newshas.html"

MAX_RANK = 4
SHA_LIST = range(2, 35) + [37, 41, 43, 47, 50, 75]


def make_rankshatable(nmax=30, verbose=False):
    total_tab = {}
    rank_tab = [{} for r in range(MAX_RANK + 1)]
    range_tab = [{} for n in range(nmax)]
    total = 0
    rank_total = [0 for r in range(MAX_RANK + 1)]
    range_total = [0 for n in range(nmax)]

    for n in range(nmax):
        infilename = ''.join(["allbigsha/allbigsha.", str(n), "0000-", str(n), "9999"])
        if verbose:
            print("processing " + infilename)
        infile = open(infilename)
        for L in infile.readlines():
            _, _, _, _, r, _, S = L.split()
            r = int(r)
            S = int(S)
            s = int(isqrt(S))
            total_tab[s] = total_tab.get(s, 0) + 1
            rank_tab[r][s] = rank_tab[r].get(s, 0) + 1
            range_tab[n][s] = range_tab[n].get(s, 0) + 1
            total += 1
            rank_total[r] += 1
            range_total[n] += 1
        infile.close()

        if verbose:
            print("Totals for range {}0000-{}9999: {} (total {})".format(n, n, range_tab[n], range_total))

    if verbose:
        print("\nTotals for all ranks: {} (total {})\n".format(total_tab, total))

    outfilename = HTML_FILENAME
    outfile = open(outfilename, mode='w')

# header info for html file

    outfile.write('<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n')
    outfile.write("<HTML>\n")
    outfile.write("<HEAD>\n")
    outfile.write('<meta http-equiv="Content-Type" content="text/html;charset=ISO-8859-1">\n')
    outfile.write("<TITLE>Analytic Shas</TITLE>\n")
    outfile.write("</HEAD>\n")
    outfile.write("<BODY>\n")
#
# First table: by ranges, not subdivided by rank
#
    outfile.write("<H3 align=center>\n")
    outfile.write("Number of curves with nontrivial (analytic) Sha, by conductor range\n")
    outfile.write("</H3>\n")
    outfile.write("<P align=center>\n")
    outfile.write("Details in files allbigsha.00000-09999 etc.\n")
    outfile.write("</P>\n")
# Table head
    outfile.write("<table border=2 align=center cellpadding=3 rules=groups>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=33>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<thead>\n")
    outfile.write('<tr style="border: solid">\n')
    outfile.write("<th>N<th>all>1")
    for s in SHA_LIST:
        outfile.write("<th>%s<sup>2</sup>" % s)
    outfile.write("\n</tr>\n")
    outfile.write("</thead>\n")
# Table foot
    outfile.write("<tfoot >\n")
    outfile.write("<tr>\n")
    outfile.write("<th align=right>1-%s9999</th>\n" % str(nmax-1))
    outfile.write("<td align=right>%s</td>\n" % total)
    for s in SHA_LIST:
        outfile.write("<td align=right>%s</td>\n" % total_tab.get(s, 0))
    outfile.write("\n</tr>\n")

    outfile.write('<tr style="border: solid">\n')
    outfile.write("<th>N<th>all>1")
    for s in SHA_LIST:
        outfile.write("<th>%s<sup>2</sup>" % s)
    outfile.write("\n</tr>\n")
    outfile.write("</tfoot>\n")
# Table body
    outfile.write("<tbody>\n")
    for n in range(nmax):
        outfile.write("<tr>\n")
        if n == 0:
            outfile.write("<th align=right>1-9999</th>\n")
        else:
            outfile.write("<th align=right>%s0000-%s9999</th>\n" % (str(n), str(n)))
        outfile.write("<td align=right>%s</td>\n" % range_total[n])
        for s in SHA_LIST:
            outfile.write("<td align=right>%s</td>\n" % range_tab[n].get(s, '&nbsp;'))
        outfile.write("</tr>\n")

    outfile.write("</tbody>\n")
    outfile.write("</table>\n")
    outfile.write("<br>\n")

#
# Second table: by ranks, not subdivided by range
#
    outfile.write("<H3 align=center>\n")
    outfile.write("Number of curves with nontrivial (analytic) Sha, by rank\n")
    outfile.write("</H3>\n")
    outfile.write("<table border=2 align=center cellpadding=3 rules=groups>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=33>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<thead>\n")
    outfile.write('<tr style="border: solid">\n')
    outfile.write("<th>&nbsp;<th>all>1")
    for s in SHA_LIST:
        outfile.write("<th>%s<sup>2</sup>" % s)
    outfile.write("\n</tr>\n")
    outfile.write("</thead>\n")
# Table foot
    outfile.write("<tfoot >\n")
    outfile.write("<tr>\n")
    outfile.write("<th align=right>all ranks</th>\n")
    outfile.write("<td align=right>%s</td>\n" % total)
    for s in SHA_LIST:
        outfile.write("<td align=right>%s</td>\n" % total_tab.get(s, 0))
    outfile.write("\n</tr>\n")

    outfile.write('<tr style="border: solid">\n')
    outfile.write("<th>&nbsp;<th>all>1")
    for s in SHA_LIST:
        outfile.write("<th>%s<sup>2</sup>" % s)
    outfile.write("\n</tr>\n")
    outfile.write("</tfoot>\n")
# Table body
    outfile.write("<tbody>\n")
    for r in range(MAX_RANK+1):
        if rank_total[r] > 0:
            outfile.write("<tr>\n")
            outfile.write("<th align=right>r=%s</th>\n" % str(r))
            outfile.write("<td align=right>%s</td>\n" % rank_total[r])
            for s in SHA_LIST:
                outfile.write("<td align=right>%s</td>\n" % rank_tab[r].get(s, '&nbsp;'))
            outfile.write("</tr>\n")
    outfile.write("</tbody>\n")
    outfile.write("</table>\n")
    outfile.write("<br>\n")

# footer info for html file

    outfile.write("<p>\n")
    outfile.write('      <a href="http://validator.w3.org/check?uri=referer"><img border="0"\n')
    outfile.write('          src="http://www.w3.org/Icons/valid-html401"\n')
    outfile.write('          alt="Valid HTML 4.01!" height="31" width="88"></a>\n')
    outfile.write("</p>\n")
    outfile.write("</BODY>\n")
    outfile.write("</HTML>\n")


    if verbose:
        for r in range(MAX_RANK + 1):
            if rank_total[r] > 0:
                print("Totals for rank {}: {} (total {})".format(r, rank_tab[r], rank_total[r]))



