import argparse
import sys

import gzip
def readFile(infile):
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin
        
def writeFile(outfile):
    if outfile.endswith((".gz","gzip")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

import pysam
def readSam(insamfile):
    if insamfile.endswith(".bam"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam.gz"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam"):
        insam = pysam.AlignmentFile(insamfile,'r')
    else:
        raise ValueError("the input sam/bam file is not end with sam or bam!")
    return insam
        
def writeSam(outsamfile,header):
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pair", default="", type=str, help="singleton file")
    parser.add_argument('-s', "--single", default="", type=str, help="single-end file")
    parser.add_argument('-o', "--outfile", required=True, help="The output file. Automatically recognize sam or bam")
    args = parser.parse_args()

    if args.pair == "" and args.single == "":
        sys.exit("please give a pair-end bam or single-end bam\n")

    fout = writeFile(args.outfile)
    fdis = writeFile("{0}.len".format(args.outfile))
    read_file_counter = 0
    frags = {}

    if args.single != "":
        for f in args.single.split(","):
            fin = readSam(f)
            for read in fin:
                read_file_counter += 1
                chrom = read.reference_name
                if read.is_reverse:
                    start = read.reference_start
                    end = read.reference_end - 5
                else:
                    start = read.reference_start + 4
                    end = read.reference_end
                rlen = end - start
                fdis.write("{0}\n".format(rlen))
                barcode = read.get_tag("CB")
                frag_id = "{0}\t{1}\t{2}\t{3}".format(chrom,start,end,barcode)
                if frag_id not in frags:
                    frags[frag_id] = 0
                frags[frag_id] += 1
            fin.close()
        
    if args.pair != "":
        for f in args.pair.split(","):
            fin = readSam(f)
            for read1 in fin:
                read2 = next(fin)
                read_file_counter += 1
                chrom = read1.reference_name
                if read1.is_reverse:
                    start = read2.reference_start + 4
                    end = read1.reference_end - 5
                else:
                    start = read1.reference_start + 4
                    end = read2.reference_end - 5
                rlen = end - start
                fdis.write("{0}\n".format(rlen))
                barcode = read1.get_tag("CB")
                frag_id = "{0}\t{1}\t{2}\t{3}".format(chrom,start,end,barcode)
                if frag_id not in frags:
                    frags[frag_id] = 0
                frags[frag_id] += 1
            fin.close()
        fdis.close()

for k in frags.keys():
    fout.write("{0}\t{1}\n".format(k,frags[k]))
fout.close()
