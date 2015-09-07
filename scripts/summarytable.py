# script used to create table.html from allbsd.* files:

#countrank.awk:

# cat curves.*0000-*9999 |
# gawk -v FIRST=$1 -v LAST=$2 'BEGIN{printf("Curve numbers by rank in the range %d...%d:\nrank:\t0\t1\t2\t3\t4\n",FIRST,LAST);}\
# ($1>=FIRST)&&($1<=LAST){r[$5]+=1;rt+=1;}\
#  END {printf("number:\t%d\t%d\t%d\t%d\t%d\nTotal  number: %d\n",
#              r[0],r[1],r[2],r[3],r[4],rt);printf("<tr>\n<th align=right>%d-%d</th>\n<td align=right>%d</td>\n",FIRST,LAST,rt);for(i=0;i<5;i++){printf("<td align=right>%d</td>\n",r[i]);}printf("</tr>\n");}'

#countcurves.awk:

# cat allcurves.*0000-*9999 |
# gawk -v FIRST=$1 -v LAST=$2 'BEGIN{printf("Numbers of isogeny and isomorphism classes in the range %d...%d:\n",FIRST,LAST);ncu=0;ncl=0;}\
# ($1>=FIRST)&&($1<=LAST){;ncu+=1;if($3==1){ncl+=1;}}\
#  END {printf("<tr>\n<th align=right>%d-%d</th>\n<td align=right>%d</td>\n<td align=right>%d</td>\n</tr>\n",FIRST,LAST,ncl,ncu);}'


# we do not call the output file "table.html" so we can compare the new
# version with the old
HTML_FILENAME = "newtable.html"

MAX_RANK = 4

def make_table(nmax=30, verbose=False):
    total_tab = {}
    rank_tab = [{} for r in range(MAX_RANK+1)]
    range_tab = [{} for n in range(nmax)]
    total = 0
    rank_total = [0 for r in range(MAX_RANK+1)]
    range_total = [0 for n in range(nmax)]
    range_total_all = [0 for n in range(nmax)]
    total_all = 0

    for n in range(nmax):
        infilename = "allcurves/allcurves."+str(n)+"0000-"+str(n)+"9999"
        range_total_all[n] = len(file(infilename).readlines())
        total_all += range_total_all[n]

    for n in range(nmax):
        infilename = "curves/curves."+str(n)+"0000-"+str(n)+"9999"
        if verbose:
            print "processing %s"%infilename
        infile = file(infilename)
        for L in infile.readlines():
            N, cl, num, ainvs, r, t, d = L.split()
            r = int(r)
            total_tab[r] = total_tab.get(r,0)+1
            range_tab[n][r] = range_tab[n].get(r,0)+1
            total +=1
            rank_total[r] +=1
            range_total[n] +=1
        infile.close()

        if verbose:
            print "Totals for range %s0000-%s9999: %s (total %s)"%(n,n,range_tab[n],range_total)

    if verbose:
        print
        print "Totals for all: %s (total %s)"%(total_tab,total)
        print

    outfilename = HTML_FILENAME
    outfile = file(outfilename, mode='w')

# header info for html file

    outfile.write('<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n')
    outfile.write("<HTML>\n")
    outfile.write("<HEAD>\n")
    outfile.write("<TITLE>Elliptic Curves</TITLE>\n")
    outfile.write("</HEAD>\n")
    outfile.write("<BODY>\n")

# Table 1 header

    outfile.write("<H3 align=center>\n")
    outfile.write("Number of isogeny classes of curves of conductor N < %s0000, sorted by rank\n"%nmax)
    outfile.write("</H3>\n")
    outfile.write("<table border=2 align=center  cellpadding=3 rules=groups>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=5>\n")
    outfile.write("<thead>\n")
    outfile.write("<tr>\n")
    outfile.write("<th>N<th>all<th>r=0<th>r=1<th>r=2<th>r=3<th>r=4\n")
    outfile.write("</tr>\n")
    outfile.write("</thead>\n")

#Table 1 footer

    outfile.write("<tfoot>\n")
    outfile.write("<tr>\n")
    outfile.write("<th align=right>1-%s9999</th>\n"%(nmax-1))
    outfile.write("<td align=right>%s</td>\n"%total)
    for r in range(MAX_RANK+1):
        outfile.write("<td align=right>%s</td>\n"%total_tab[r])
    outfile.write("</tr>\n")
    outfile.write("</tfoot>\n")

#Table 1 body

    outfile.write("<tbody>\n")
    for n in range(nmax):
        outfile.write("<tr>\n")
        if n==0:
            outfile.write("<th align=right>1-9999</th>\n")
        else:
            outfile.write("<th align=right>%s0000-%s9999</th>\n"%(n,n))
        outfile.write("<td align=right>%s</td>\n"%range_total[n])
        for r in range(MAX_RANK+1):
            outfile.write("<td align=right>%s</td>\n"%range_tab[n].get(r,0))
        outfile.write("</tr>\n")

    outfile.write("</tbody>\n")
    outfile.write("</table>\n")

    outfile.write("<P></P>\n")
    outfile.write("<HR>\n")
    outfile.write("<P></P>\n")



# Table 2 header

    outfile.write("<H3 align=center>\n")
    outfile.write("Total number of curves of conductor N < %s0000\n"%nmax)
    outfile.write("</H3>\n")
    outfile.write("<table border=2 align=center cellpadding=3 rules=groups>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<colgroup span=1>\n")
    outfile.write("<thead>\n")
    outfile.write("<tr>\n")
    outfile.write("<th>N</th>\n")
    outfile.write("<th># isogeny classes</th>\n")
    outfile.write("<th># isomorphism classes</th>\n")
    outfile.write("</tr>\n")
    outfile.write("</thead>\n")

# Table 2 footer

    outfile.write("<tfoot>\n")
    outfile.write("<tr>\n")
    outfile.write("<th align=right>1-%s9999</th>\n"%(nmax-1))
    outfile.write("<td align=right>%s</td>\n"%total)
    outfile.write("<td align=right>%s</td>\n"%total_all)
    outfile.write("</tr>\n")
    outfile.write("</tfoot>\n")

# Table 2 body

    outfile.write("<tbody>\n")
    for n in range(nmax):
        outfile.write("<tr>\n")
        if n==0:
            outfile.write("<th align=right>1-9999</th>\n")
        else:
            outfile.write("<th align=right>%s0000-%s9999</th>\n"%(n,n))
        outfile.write("<td align=right>%s</td>\n"%range_total[n])
        outfile.write("<td align=right>%s</td>\n"%range_total_all[n])
        outfile.write("</tr>\n")

    outfile.write("</tbody>\n")
    outfile.write("</table>\n")


# html footer

    outfile.write("<P></P>\n")
    outfile.write("<HR>\n")
    outfile.write("<P></P>\n")
    outfile.write("    <p>\n")
    outfile.write('      <a href="http://validator.w3.org/check?uri=referer"><img border="0"\n')
    outfile.write('          src="http://www.w3.org/Icons/valid-html401"\n')
    outfile.write('          alt="Valid HTML 4.01!" height="31" width="88"></a>\n')
    outfile.write("    </p>\n")
    outfile.write("</BODY>\n")
    outfile.write("</HTML>\n")

    outfile.close()
