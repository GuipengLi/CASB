import pysam
import argparse
import os
import gzip
import regex
'''
extract barcodes according to read1(first 26bp: 16bp 10X bc, 10bp UMI) and read2(8bp sample bc) prefix.
perfect match with cell barcodes.
allow 1 mismatch to sample barcode.
output UMI count.
'''

def in_sample(bc, samples):
    for sp in samples:
        p=regex.findall(r'(%s){s<=1}'%bc,sp)
        if p:
            return(sp)
    return(bc)
    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-r1", "--read1",default="test.R1.fq.gz", help="reads 1")
    parser.add_argument("-r2", "--read2",default="test.R2.fq.gz", help="reads 2")
    parser.add_argument("-b", "--barcode", default="LF_CABS_scRNA_CR_merge/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", help="scRNA barcodes")
    parser.add_argument("-s", "--sample", default="samples_barcode.txt", help="sample barcodes")
    parser.add_argument("-o", "--out", default="output", help="output file prefix name")
    args = parser.parse_args()
    read1 = args.read1
    read2 = args.read2
    outfile = args.out
    barcodes = list()
    with gzip.open(args.barcode, 'r') as fin:
        for line in fin:
            bc=line.strip().split('-')[0]
            if bc not in barcodes:
                barcodes.append(bc)
    samples = list()
    with open(args.sample, 'r') as fin:
        for line in fin:
            bc=line.strip().split('\t')[2]
            if bc not in samples:
                samples.append(bc)

    mat = {}
    c,i,j,k = 0,0,0,0
    with pysam.FastxFile(read1) as fin_a, pysam.FastxFile(read2) as fin_b:
        for reada in fin_a:
            readb = fin_b.next()
            c += 1
            if c % 5000000 == 0:
                print('processed %d' %c)
            if reada.name != readb.name:
                print("read names not equal.")
                os.exit()
            else:
                read1 = reada.sequence[0:26]
                sbc = readb.sequence[20:28]
                cbc = read1[0:16]
                umi = read1[16:26]
                line = cbc + '\t'+ umi + '\t' + sbc
                if line not in mat:
                    i += 1
                    mat[line] = 1
                else:
                    mat[line] += 1

    with open(outfile, 'w') as fout, open(outfile+'.scRNA_flt','w') as fout2:
        for line in mat:
            fout.write(line+'\t'+str(mat[line])+'\n')
            bc = line.split('\t')[0]
            sbc = line.split('\t')[2]
            umi = line.split('\t')[1]
            if bc in barcodes:
                sbc = in_sample(sbc, samples)
                newline = bc + '\t'+ umi + '\t' + sbc
                j += 1
                fout2.write(newline+'\t'+str(mat[line])+'\n')
                
    print('Total: %d, combination: %d, scRNA: %d' %(c,i,j))

if __name__=="__main__":
    main() 
