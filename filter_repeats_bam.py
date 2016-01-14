#!usr/bin/python

# read bam file from stream

import re
import sys, getopt

def distance(infile, outfile=None, maxDist=1e5, multi=False, nosingle=False, debug=False):
    if infile != 'stdin':
        raise ValueError('Function can only be used with standard input')

    if outfile != None:
        sam_match = re.search('\.sam',outfile)
    
    if outfile != None:
        o = open(outfile, 'w')
        
    read = ''
    read2 = ''
    readbin = []
    readcount = 0
    validcount = 0
    while True:
        line = sys.stdin.readline()        
        
        if not line:
            break

        # skip header lines and write to file
        header_match = re.match('@',line)
        if header_match != None:
            if outfile != None:
                o.write(line)   # write header to file if output is SAM format
            else:
                print line.replace('\n','')
            continue
        
        readcount += 1
        line2 = line.replace('\n','')
        cols = line2.split('\t')
        cols[1] = int(cols[1])

        # skip unmapped reads
        unmapped = re.search('\*',cols[2])
        if unmapped != None:
            continue

        if read == '':
            read = cols[0]
            read2 = read + '/2'
            readbin += [line2]
        elif read != cols[0] and read2 != cols[0]:
            if debug:
                print '\n# unparsed alignments:', len(readbin)
                
            final, single = get_valid_align(readbin, maxDist, debug)

            if final != None:
                if len(final) > 1 and not multi:
                    if debug:
                        print 'No unique alignments!'
                else:
                    #validcount += len(final)

                    if outfile != None:
                        if not single:      # if not singleton
                            validcount += output_sam(final, o)
                        else:
                            if not nosingle:
                                validcount += len(final)
                                for i in xrange(len(final)):
                                    o.write('\t'.join(final[i]) + '\n')
                    else:
                        for i in xrange(len(final)):
                            if not single:
                                validcount += output_sam(final)
                            else:
                                if not nosingle:
                                    validcount += len(final)
                                    for i in xrange(len(final)):
                                        print '\t'.join(final[i])
            else:
                if debug:
                    print '\nNo valid alignments!'
            
            read = cols[0]
            read2 = read + '/2'
            readbin = [line2]
        else:
            readbin += [line2]

        if readcount % 5e5 == 0:
            if outfile != None:
                print str(readcount) + ' alignments processed ...'
            else:
                sys.stderr.write(str(readcount) + ' alignments processed ...\n')

        if debug:
            if validcount > 100:
                if outfile != None:
                    o.close()
                break

    if debug:
        if outfile != None:
            o.close()
        return
    
    if len(readbin) > 0:
        final, single = get_valid_align(readbin, maxDist, debug)
        if final != None:
            if len(final) > 1 and not multi:
                if debug:
                    print '\nNo valid alignments!'
            else:
                #validcount += len(final)

                if outfile != None:
                    if not single:      # if not singleton
                        validcount += output_sam(final, o)
                    else:
                        if not nosingle:
                            validcount += len(final)
                            for i in xrange(len(final)):
                                o.write('\t'.join(final[i]) + '\n')
                else:                    
                    if not single:
                        validcount += output_sam(final)
                    else:
                        if not nosingle:
                            validcount += len(final)
                            for i in xrange(len(final)):
                                print '\t'.join(final[i])
    if outfile != None:
        o.close()
    else:
        sys.stderr.write('\nTotal alignments processed:' + str(readcount) + '\n')
        sys.stderr.write('# valid alignments:' + str(validcount) + '\n')

# function to parse alignments
# output sam file with reads in valid alignments
def get_valid_align(readbin, maxDist=1e5, debug=False):
    # lists to separate by left/right read, forward/backward alignment
    lfreads = []
    lbreads = []
    rfreads = []
    rbreads = []
    single = False  # flag to denote singleton alignment(s)
    
    for i in xrange(len(readbin)):        
        cols = readbin[i].split('\t')
        id_match = re.search('/2',cols[0])

        # criteria for left and right reads:
        # right read id has pattern '/2'
        if id_match == None:    # left read
            if int(cols[1]) & 0x10 == 0:    # forward alignment
                lfreads += [cols]
            else:
                lbreads += [cols]
        else:                   # right read
            if int(cols[1]) & 0x10 == 0:    # forward alignment
                rfreads += [cols]
            else:
                rbreads += [cols]

    if debug:
        print '# forward alignments for read 1:',len(lfreads),
        print '# backward alignments for read 1:', len(lbreads)
        print '# forward alignments for read 2:',len(rfreads),
        print '# backward alignments for read 2:', len(rbreads)

    # singleton alignments
    # pick alignment(s) with minimum mismatches
    sreads = []
    if len(rfreads) + len(rbreads) == 0:
        sreads = lfreads + lbreads
    elif len(lfreads) + len(lbreads) == 0:
        sreads = rfreads + rbreads

    if len(sreads) > 0:
        single = True
        s_mm = [ int(re.search('\d', sreads[a][11]).group(0)) \
                 for a in xrange(len(sreads)) ]
        smin_mm = min(s_mm)
        filt_indices = [ a for a in xrange(len(sreads)) if \
                         s_mm[a] == smin_mm]
        sfilt = [sreads[a] for a in filt_indices]

        if debug:
            print '# of valid singleton alignments:',len(sfilt)
        return sfilt, single
    
    valid = []
    lvalid = []
    rvalid = []
    lmin_mm = 10
    rmin_mm = 10
    
    # find suitable partners for left forward reads    
    for i in xrange(len(lfreads)):
        [chrom, start, tags] = [lfreads[i][2],lfreads[i][3], lfreads[i][11]]
        l_mm = int(re.search('\d',tags).group(0))   # mismatches in left forward read
        
        # get right backward reads on same chromosome
        mates = [ a for a in xrange(len(rbreads)) if rbreads[a][2] == chrom ]
        if len(mates) == 0:
            continue

        rlen = [ re.search('\d+',rbreads[a][5]).group(0) \
                 for a in xrange(len(rbreads))] # length of reads
        dist = [ (float(rbreads[a][3]) + float(rlen[a]) - float(start)) for a in mates ]    # distances between pairs
        r_mm = [ int(re.search('\d', rbreads[a][11]).group(0)) for a in mates] # mismatches in right mate
        sum_mm = [ (a + l_mm) for a in r_mm ]  # total mismatches in pair
        lmin_mm = min(sum_mm)
        
        # get suitable mates: proper distance and minimum mismatches
        filt_indices = [ a for a in xrange(len(mates)) if  \
                         dist[a] < maxDist and dist[a] > 0 and \
                         sum_mm[a] == lmin_mm ]
        filtmates = [ mates[a] for a in filt_indices ]
        
        if len(filtmates) > 0:
            for j in xrange(len(filtmates)):
                temp = [ '\t'.join(lfreads[i]), '\t'.join(rbreads[ filtmates[j] ]), \
                         str(dist[ filt_indices[j] ]), str(sum_mm[ filt_indices[j] ])]
                lvalid += [temp]

    # find suitable partners for right forward reads
    for i in xrange(len(rfreads)):
        [chrom, start, tags] = [rfreads[i][2], rfreads[i][3], rfreads[i][11]]
        r_mm = int(re.search('\d',tags).group(0))   # mismatches in right forward read
        
        # get right backward reads on same chromosome
        mates = [ a for a in xrange(len(lbreads)) if lbreads[a][2] == chrom ]
        if len(mates) == 0:
            continue

        llen = [ re.search('\d+',lbreads[a][5]).group(0) \
                 for a in xrange(len(lbreads))] # length of reads
        dist = [ (float(lbreads[a][3]) + float(llen[a]) - float(start)) for a in mates ]    # distances between pairs
        l_mm = [ int(re.search('\d', lbreads[a][11]).group(0)) for a in mates] # mismatches in right mate
        sum_mm = [ (a + r_mm) for a in l_mm ]  # total mismatches in pair
        rmin_mm = min(sum_mm)
        
        # get suitable mates: proper distance and minimum mismatches
        filt_indices = [ a for a in xrange(len(mates)) if  \
                         dist[a] < maxDist and dist[a] > 0 and \
                         sum_mm[a] == rmin_mm ]
        filtmates = [ mates[a] for a in filt_indices ]
        
        if len(filtmates) > 0:
            for j in xrange(len(filtmates)):
                temp = [ '\t'.join(rfreads[i]), '\t'.join(lbreads[ filtmates[j] ]), \
                         str(dist[ filt_indices[j] ]), str(sum_mm[ filt_indices[j] ])]
                rvalid += [temp]
                
    if lmin_mm < rmin_mm:
        valid = lvalid
    elif rmin_mm < lmin_mm:
        valid = rvalid
    elif lmin_mm == rmin_mm:
        valid = lvalid + rvalid

    if debug:
        print '# of valid alignments:',len(valid)
    
    if len(valid) == 0:
        return None, single
    else:
        return valid, single

def filter_mismatch(alignlist, debug=False):
    min_mismatch = 10
    for i in xrange(len(alignlist)):
        tags = [a for a in alignlist[i] if re.search('XA',a) != None]
        if len(tags) == 0:
            alignlist[i].append('0')
        else:
            mismatch = [int(re.search('\d',a).group(0)) for a in tags]
            temp_sum = sum(mismatch)
            if temp_sum < min_mismatch:
                min_mismatch = temp_sum
            alignlist[i].append( str( temp_sum ) )

    filtlist = []
    for i in xrange(len(alignlist)):
        if int(alignlist[i][-1]) == min_mismatch:
            filtlist += [alignlist[i]]
    if debug:
        print '# of filtered alignments:',len(filtlist)
    return filtlist

def output_sam(align, handle=None):
    sam_align = []
    for i in xrange(len(align)):
        sam_align += [align[i][0]] + [align[i][1]]
    sam_align = sorted(list(set(sam_align)))    # remove duplicates if any
    
    if handle != None:
        handle.write('\n'.join(sam_align) + '\n')
    else:
        print '\n'.join(sam_align)
    return len(sam_align)
    
def usage():
    print 'python parse_repeats.py -i <infile> -o <outfile>',
    print '[-d/--distance maxDist] [-n/--debug] [-m/--multi] [-s/--nosingle]'
    
def main(argv):
    infile = ''
    outfile = None
    mode = ''
    debug = False
    maxDist = 1e5
    multi = False
    nosingle = False
    try:
        opts, args = getopt.getopt(argv,'hi:o:d:nm',\
                                   ['help','ifile=','ofile=',\
                                    'distance=','debug','multi','nosingle'])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)
        
    for opt, arg in opts:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ('-i','--ifile'):
            infile = arg
        elif opt in ('-o','--ofile'):
            outfile = arg
        elif opt in ('-d','--distance'):
            mode = 'distance'
            maxDist = float(arg)
        elif opt in ('-n','--debug'):
            debug = True
        elif opt in ('-m','--multi'):
            multi = True
        elif opt in ('-s','--nosingle'):
            nosingle = True
   
    if mode == 'distance':
        distance('stdin', outfile, maxDist, multi, nosingle, debug)   
    
if __name__ == '__main__':
    main(sys.argv[1:])
