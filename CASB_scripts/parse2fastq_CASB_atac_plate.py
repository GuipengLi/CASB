import pysam
import argparse
import os
import gzip
import regex
'''
read 1 contain barcodes
read 2 contain gRNA
read 1 contain replicate index at the begining
'''

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return (bases)
    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-r1", "--read1",default="test.R1.fq.gz", help="reads 1")
    parser.add_argument("-r2", "--read2",default="test.R2.fq.gz", help="reads 2")
    parser.add_argument("-o", "--outfile",default="cellbarcode_cabs.txt", help="output file")
    parser.add_argument("-p", "--prefix",default="output prefix", help="output fastq prefix")
    args = parser.parse_args()
    read1 = args.read1
    read2 = args.read2
    outfiles = [ args.prefix+"_"+x for x in [read1,read2] ]
    outf = args.outfile
    grnadict = {}
    c,n,i,j,k = 0,0,0,0,0
    kk = 0
    with pysam.FastxFile(read1) as fin_1, pysam.FastxFile(read2) as fin_2, open(outf, 'w') as fout, gzip.open(outfiles[0],'wb') as fout1, gzip.open(outfiles[1],'wb') as fout2:
        for read1 in fin_1:
            read2 = fin_2.next()
            c += 1
            if c % 100000 == 0:
                print('processed %d' %c)
            if read1.name != read2.name:
                print("read names not equal.")
                os.exit()
            else:
                reada = read1.sequence
                p2=regex.findall(r'([ATCG]{20}GTTTTAGAGCT){s<=1}',reada)
                if p2:
                    i += 1
                    cabs = p2[0][-30:-11]
                    line = cabs+'\n'
                    fout.write(line)
                else:
                    fout1.write(str(read1)+'\n')
                    fout2.write(str(read2)+'\n')
                    j += 1

    print('processed %d' %c)
    print('extract_plate, %d, %d, %d' %(c,i,j))

if __name__=="__main__":
    main() 
