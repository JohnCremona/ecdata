# Script to fix Lr1 and Sha values when we had forgotten to divide by r! (361725<=N<=370000)

squares = [n*n for n in srange(1,100)]

def fix_bsd(infilename):
    """
    If rank is >1 check to see if Lr1 and Sha need to be divided by r!

    Input line fields:

    conductor iso number ainvs rank torsion tamagawa_product real_period special_value regulator sha_an

    Sample input line:

    11 a 1 [0,-1,1,-10,-20] 0 5 5 1.2692093042795534217 0.25384186085591068434 1 1.00000000000000000000
    """
    infile = file(infilename,mode='r')
    outfile = file(infilename+".fixed", mode='w')
    for line in file(infilename).readlines():
        data = line.split()
        r=int(data[4])
        s = round(float(data[10]))
        if not s in squares:
            print("fixing %s%s%s: rank = %s, Sha = %s" % (data[0],data[1],data[2],data[4],data[10]))
            rfact = factorial(r)
            data[10]=str(RR(data[10])/rfact)
            data[8]=str(RR(data[8])/rfact)
            print("fixed, new Sha = %s" % data[10])
            line = ' '.join(data)+'\n'
        outfile.write(line)

    infile.close()
    outfile.close()

