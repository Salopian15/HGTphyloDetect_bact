import os
import sys
import time
from Bio import SeqIO


def main() :
    """
    Performs a remote Blastp search against NCBI online database for each gene in the input file.

    Args:
        Infile.fasta: Must be a protein fasta file.

    Returns:
        Blastp search file in tabular format for each individual protein sequence.
    """
    start_time = time.time()
    name = sys.argv[1:][0]
    genes = list()
    geneSeq = dict()
    with open(name, 'r') as handleGene :
        for record in SeqIO.parse(handleGene, "fasta") :
            gene = str(record.id)
            sequence = str(record.seq)
            geneSeq[gene] = sequence
            genes.append(gene)

    for gene in genes :
        print('This is gene:', gene)
        blastp_file = './blastp_files/%s.txt' % gene
        if not os.path.exists(blastp_file) or os.stat(blastp_file).st_size == 0:
            with open('./%s.fasta' % gene, 'w') as outfile :
                outfile.write('>'+gene+'\n')
                outfile.write(geneSeq[gene]+'\n')
            if not os.path.exists("./blastp_files/") :
                os.system('mkdir ./blastp_files')
            # Execute shell command line in Python using the os.system function
            myCmd = "blastp -db nr -remote -query ./%s.fasta -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue bitscore length pident' -out %s" % (gene, blastp_file)
            os.system(myCmd)
            os.remove('./%s.fasta' % gene)
        else:
            print('Blastp file already exists for gene:', gene)

    elapsed_time = (time.time() - start_time)
    print("Finished------------------")
    print("Running the Blast process: %.4fs" %(elapsed_time))


if __name__ == "__main__" :
    main()
