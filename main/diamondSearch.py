import os
import sys
import time
from Bio import SeqIO


def main() :
    """
    This function performs a diamond search for genes in a given file.

    The function reads a file containing gene sequences in FASTA format,
    performs a diamond search for each gene, and saves the results in separate
    files. If a diamond file already exists for a gene, it skips the search
    for that gene.

    Args:
        Infile.fasta: A file containing protein sequences in FASTA format.

    Returns:
        Diamond search file in tabular format for each individual protein sequence.

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
        diamond_file = './diamond_files/%s.txt' % gene
        if not os.path.exists(diamond_file) or os.stat(diamond_file).st_size == 0:
            with open('./%s.fasta' % gene, 'w') as outfile :
                outfile.write('>'+gene+'\n')
                outfile.write(geneSeq[gene]+'\n')
            if not os.path.exists("./diamond_files/") :
                os.system('mkdir ./diamond_files')
            print("Running Diamond")
            # Dependent on diamond database being installed
            myCmd = "diamond blastp -d /ibers/ernie/scratch/jac180/diamond_db/uniref50db.dmnd -q ./%s.fasta --max-target-seqs 250 --outfmt 6 qseqid sseqid evalue bitscore length pident staxids sskingdoms sphylums -o %s" % (gene, diamond_file)
            os.system(myCmd)
            os.remove('./%s.fasta' % gene)
        else:
            print('Diamond file already exists for gene:', gene)

    elapsed_time = (time.time() - start_time)
    print("Finished------------------")
    print("Running the Diamond process: %.4fs" %(elapsed_time))


if __name__ == "__main__" :
    main()








