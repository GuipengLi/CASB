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
indexseq = ['ACAGTG','GCCAAT','CAGATC','ACTTGA']

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return (bases)
    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-i1", "--read0",default="test.I1.fq.gz", help="reads 0")
    parser.add_argument("-r1", "--read1",default="test.R1.fq.gz", help="reads 1")
    parser.add_argument("-r2", "--read2",default="test.R2.fq.gz", help="reads 2")
    parser.add_argument("-r3", "--read3",default="test.R3.fq.gz", help="reads 3")
    parser.add_argument("-o", "--outfile",default="cellbarcode_cabs.txt", help="output file")
    parser.add_argument("-p", "--prefix",default="output prefix", help="output fastq prefix")
    args = parser.parse_args()
    read0 = args.read0
    read1 = args.read1
    read2 = args.read2
    read3 = args.read3
    outfiles = [ args.prefix+"_"+x for x in [read0,read1,read2,read3] ]
    outf = args.outfile
    grnadict = {}
    c,n,i,j,k = 0,0,0,0,0
    kk = 0
    with pysam.FastxFile(read0) as fin_0, pysam.FastxFile(read1) as fin_1, pysam.FastxFile(read2) as fin_2, pysam.FastxFile(read3) as fin_3, open(outf, 'w') as fout, gzip.open(outfiles[0],'wb') as fout0, gzip.open(outfiles[1],'wb') as fout1, gzip.open(outfiles[2],'wb') as fout2, gzip.open(outfiles[3],'wb') as fout3:
        for read0 in fin_0:
            read1 = fin_1.next()
            read2 = fin_2.next()
            read3 = fin_3.next()
            c += 1
            if c % 5000000 == 0:
                print('processed %d' %c)
            if read1.name != read3.name:
                print("read names not equal.")
                os.exit()
            else:
                reada = read1.sequence
                cb = read2.sequence
                p2=regex.findall(r'([ATCG]{20}GTTTTAGAGCT){s<=1}',reada)
                if p2:
                    i += 1
                    cabs = p2[0][-30:-11]
                    line = cb+'\t'+cabs+'\n'
                    fout.write(line)
                else:
                    fout0.write(str(read0)+'\n')
                    fout1.write(str(read1)+'\n')
                    fout2.write(str(read2)+'\n')
                    fout3.write(str(read3)+'\n')
                    j += 1

    print('processed %d' %c)
    print('%d, %d, %d' %(c,i,j))

if __name__=="__main__":
    main() 
