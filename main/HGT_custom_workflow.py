#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys, os, warnings, math, csv, argparse
import pandas as pd
from Bio import SeqIO, BiopythonWarning, Entrez
from ete3 import NCBITaxa



# Path to databases need to be specified
# Not terribly sure how I go about doing this for the time being
# Perhaps just leave it to the user to specify the database paths??



def parse_arguments():
    parser = argparse.ArgumentParser(description="Modified version of HGTPhyloDetect close workflow for HGT events, takes protein fasta file and iterates through each sequence outputting a likelihood of HGT origin for each", epilog="Author: Jack A. Crosby, Aberystwyth University/Queens University Belfast")
    parser.add_argument("input_file", help="Input file path, should be a fasta file of protein sequences")
    parser.add_argument("--bitscore_parameter", type=float, default=100, help="Bitscore parameter, default is 100")
    parser.add_argument("--HGTIndex", type=float, default=0.5, help="HGT Index, default is 0.5")
    parser.add_argument("--out_pct", type=float, default=0.8, help="Out Pct, default is 0.8")
    parser.add_argument("--tax_level", type=str, default="family", choices=["superkingdom", "kingdom", "phylum", "subphylum", "class", "order", "family", "genus", "species"], help="Taxonomic level, organisms outisde of this level will be classified as 'outgroup', default is family.")
    parser.add_argument("--search", type=str, default="blastp", choices=["blastp", "diamond", "mmseqs"], help="Search methods, blastp uses remote nr search, diamond & mmseqs use local database for search, default is blastp.")
    return parser.parse_args()




def parse_NCBI(gene):
    try :
        with open("./blastp_files/%s.txt" % gene, "r") as filename :
            index = None
            data = filename.readlines()
            for line in data :
                if line.strip("\n").endswith("found") :
                    index = data.index(line)

            blast_results = data[index+1:-1]
    except Exception as e:
        print(Exception, e)
        return None
    
    accession_number = list()
    accession_bitscore = dict()
    for blast in blast_results :
        accession = blast.strip("\n").split("\t")[1]
        accession_number.append(accession)
        accession_bitscore[accession] = float(blast.strip('\n').split("\t")[3])
    return accession_number, accession_bitscore



def getTaxid(accession):
    # Retrieving data in the GenBank using only the GenBank code accession in biopython
    # https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    # https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    # https://biopython.org/DIST/docs/api/Bio.Entrez-module.html
    Entrez.email = "abcd@ncbi.org" # jac180@aber.ac.uk

    # https://www.biostars.org/p/304175/
    # get tax id using only the GenBank code accession in biopython
    # handle = Entrez.efetch(db='protein', id="NP_012706.1", rettype='gb')
    handle = Entrez.efetch(db='protein', id=accession, rettype='gb')
    record = SeqIO.read(handle,'genbank')
    # print(record.features[0].qualifiers)
    if record.features[0].qualifiers['db_xref'][0].split(":")[0] == 'taxon':
        taxid = record.features[0].qualifiers['db_xref'][0].split(":")[1] # the type is a string
        organism = record.features[0].qualifiers['organism'][0]
    # seq = record.seq
    # print(taxid,organism)

    return taxid



def protein_search(gene, search, geneSeq):
    match search:
        case "blastp":
        
            if os.path.exists("./blastp_files/%s.txt" % gene) :
                accession_number,accession_bitscore = parse_NCBI(gene)
                print('Blast file exists')
                return accession_number,accession_bitscore
            else :
                # Need to install blast!
                with open('./%s.fasta' % gene, 'w') as outfile :
                    outfile.write('>'+gene+'\n')
                    outfile.write(geneSeq[gene]+'\n')
                if not os.path.exists("./blastp_files/") :
                    os.system('mkdir ./blastp_files')
                myCmd = "blastp -db nr -remote -query ./%s.fasta -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue bitscore length pident' -out ./blastp_files/%s.txt" %(gene, gene)
                os.system(myCmd)
                os.remove('./%s.fasta' % gene)
                accession_number,accession_bitscore = parse_NCBI(gene)
                return accession_number,accession_bitscore
        
        case "diamond":
            
            if os.path.exists("./diamond_files/%s.txt" % gene) :
                #accession_number,accession_bitscore = parse_NCBI(gene)
                print('Yes, diamond file already exists, nice!')
            else :
                # Need to install diamond!
                with open('./%s.fasta' % gene, 'w') as outfile :
                    outfile.write('>'+gene+'\n')
                    outfile.write(geneSeq[gene]+'\n')
            if not os.path.exists("./diamond_files/") :
                os.system('mkdir ./diamond_files')
            myCmd = "diamond blastp -d nr -q ./%s.fasta --max-target-seqs 250 -outfmt 6 qseqid sseqid evalue bitscore length pident -o ./diamond_files/%s.txt" %(gene, gene)
            os.system(myCmd)
            os.remove('./%s.fasta' % gene)
            accession_number,accession_bitscore = parse_NCBI(gene)
            return accession_number,accession_bitscore
        
        case "mmseqs":
            
            if os.path.exists("./mmseqs_files/%s.txt" % gene) :
                #accession_number,accession_bitscore = parse_NCBI(gene)
                print('Yes, mmseqs file already exists, nice!')
            else :
                
                with open('./%s.fasta' % gene, 'w') as outfile :
                    outfile.write('>'+gene+'\n')
                    outfile.write(geneSeq[gene]+'\n')
            if not os.path.exists("./mmseqs_files/") :
                os.system('mkdir ./mmseqs_files')
            myCmd = "mmseqs search ./%s.fasta /path/to/nr_database ./mmseqs_files/%s.tsv --format-output 'query, target, evalue, bitscore, length, pident, taxlineage'" %(gene, gene)
            os.system(myCmd)
            os.remove('./%s.fasta' % gene)
            accession_number,accession_bitscore = parse_NCBI(gene)
            return accession_number,accession_bitscore
        
        case _:
            
            print(f"Unknown search method: {search}")
            sys.exit(1)
    
    
    
def check_search(gene, search, geneSeq):
        match search:
            case "blastp":
                if os.path.exists("./blastp_files/%s.txt" % gene) :
                    accession_number,accession_bitscore = parse_NCBI(gene)
                    print('Blast file exists')
                else :
                    accession_number,accession_bitscore = protein_search(gene, search, geneSeq)
            case "diamond":
                if os.path.exists("./diamond_files/%s.txt" % gene) :
                    accession_number,accession_bitscore = parse_NCBI(gene)
                    print('Diamond file exists')
                else :
                    accession_number,accession_bitscore = protein_search(gene, search, geneSeq)
            case "mmseqs":
                if os.path.exists("./mmseqs_files/%s.txt" % gene) :
                    accession_number,accession_bitscore = parse_NCBI(gene)
                    print('MMseqs file exists')
                else :
                    accession_number,accession_bitscore = protein_search(gene, search, geneSeq)
        return accession_number,accession_bitscore


        
def alt_taxid(accession, search, gene):
    
    match search:
        
        case "blastp":
            # Continue as normal
            print("Blast")
            exit()
        case "diamond":
            # Grab Taxonomy from the diamond files
            print("Diamond")
            
        case "mmseqs":
            # Grab Taxonomy from the mmseqs files
            print("MMseqs")
            def parse_mmseqs_file(gene):
                file_path = f"./mmseqs_files/{gene}.tsv"
                taxonomic_info = []
                with open(file_path, 'r') as file:
                    for line in file:
                        line = line.strip().split('\t')
                        taxonomic_info.append(line[-1])
                return taxonomic_info

            accession_number, accession_bitscore = parse_NCBI(gene)
            taxonomic_info = parse_mmseqs_file(gene)
        case _:
            print(f"Unknown search method: {search}")
            sys.exit(1)



def search_protein_comb(gene, search, geneSeq):
    
    match search:
    
        case "blastp":
    
            if os.path.exists("./blastp_files/%s.txt" % gene):
                accession_number, accession_bitscore = parse_NCBI(gene)
                print('Blast file exists')
            else:
                with open('./%s.fasta' % gene, 'w') as outfile:
                    outfile.write('>' + gene + '\n')
                    outfile.write(geneSeq[gene] + '\n')
                if not os.path.exists("./blastp_files/"):
                    os.system('mkdir ./blastp_files')
                myCmd = "blastp -db nr -remote -query ./%s.fasta -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue bitscore length pident' -out ./blastp_files/%s.txt" % (gene, gene)
                os.system(myCmd)
                os.remove('./%s.fasta' % gene)
                accession_number, accession_bitscore = parse_NCBI(gene)
        
        case "diamond":
       
            if os.path.exists("./diamond_files/%s.txt" % gene):
                accession_number, accession_bitscore = parse_NCBI(gene)
                print('Diamond file exists')
            else:
                with open('./%s.fasta' % gene, 'w') as outfile:
                    outfile.write('>' + gene + '\n')
                    outfile.write(geneSeq[gene] + '\n')
                if not os.path.exists("./diamond_files/"):
                    os.system('mkdir ./diamond_files')
                myCmd = "diamond blastp -d /databasepath -q ./%s.fasta --max-target-seqs 250 -outfmt 6 qseqid sseqid evalue bitscore length pident staxids -o ./diamond_files/%s.txt" % (gene, gene)
                os.system(myCmd)
                os.remove('./%s.fasta' % gene)
                accession_number, accession_bitscore = parse_NCBI(gene)
        
        case "mmseqs":
          
            if os.path.exists("./mmseqs_files/%s.txt" % gene):
                accession_number, accession_bitscore = parse_NCBI(gene)
                print('MMseqs file exists')
            else:
                with open('./%s.fasta' % gene, 'w') as outfile:
                    outfile.write('>' + gene + '\n')
                    outfile.write(geneSeq[gene] + '\n')
                if not os.path.exists("./mmseqs_files/"):
                    os.system('mkdir ./mmseqs_files')
                myCmd = "mmseqs search ./%s.fasta /path/to/nr_database ./mmseqs_files/%s.tsv --format-output 'query, target, evalue, bitscore, length, pident, taxlineage'" % (gene, gene)
                os.system(myCmd)
                os.remove('./%s.fasta' % gene)
                accession_number, accession_bitscore = parse_NCBI(gene)
        
        case _:
            
            print(f"Unknown search method: {search}")
            print("Please choose from blastp, diamond or mmseqs")
            print("Exiting...")
            sys.exit(1)
    
    return accession_number, accession_bitscore

    

def main(bitscore_parameter=100, HGTIndex=0.5, out_pct=0.8) :
    """
    Runs HGTPhyloDetect but slightly modified, takes arguments regarding the
    taxonomic level to determine the outgroup and the search method to use.
        
    """
    
    #name = sys.argv[1:][0]
    parse_arguments()
    #check_args(sys.argv)
    name = parse_arguments().input_file
    bitscore_parameter = parse_arguments().bitscore_parameter
    HGTIndex = parse_arguments().HGTIndex
    out_pct = parse_arguments().out_pct
    tax_level = parse_arguments().tax_level.lower()
    search = parse_arguments().search.lower()
    genes = list()
    geneSeq = dict()
    HGT = list()
    warnings.simplefilter('ignore', BiopythonWarning)

    with open(name, 'r') as handleGene :
        for record in SeqIO.parse(handleGene, "fasta") :
            gene = str(record.id)
            sequence = str(record.seq)
            geneSeq[gene] = sequence
            genes.append(gene)

    n=0
    for gene in genes :
        n += 1
        print("This is gene %d------------------" %(n))
        print(gene)
    
        if os.path.exists("./blastp_files/%s.txt" % gene) :
            accession_number,accession_bitscore = parse_NCBI(gene)
            print('Blast file exists')
            
        elif os.path.exists("./diamond_files/%s.txt" % gene) :
            accession_number,accession_bitscore = parse_NCBI(gene)
            print('Diamond file exists')
            
        elif os.path.exists("./mmseqs_files/%s.txt" % gene) :
            accession_number,accession_bitscore = parse_NCBI(gene)
            print('MMseqs file exists')
            
        else :
            accession_number, accession_bitscore= protein_search(gene, search, geneSeq)
            # Need to install blast!
            #with open('./%s.fasta' % gene, 'w') as outfile :
                #outfile.write('>'+gene+'\n')
                #outfile.write(geneSeq[gene]+'\n')
            #if not os.path.exists("./blastp_files/") :
                #os.system('mkdir ./blastp_files')
            #myCmd = "blastp -db nr -remote -query ./%s.fasta -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue bitscore length pident' -out ./blastp_files/%s.txt" %(gene, gene)
            #os.system(myCmd)
            #os.remove('./%s.fasta' % gene)
            #accession_number,accession_bitscore = parse_NCBI(gene)

        ncbi = NCBITaxa()
        recipient_accession = list()
        outgroup_accession = list()
        recipient_accession_bitscore = dict()
        outgroup_accession_bitscore = dict()
        recipient_species = list()
        outgroup_species = list()

        try :
            
            # Pull taxonomic information from the gene accessions 
            
            gene_taxid = getTaxid(gene)
            gene_lineage = ncbi.get_lineage(gene_taxid)
            gene_lineage2ranks = ncbi.get_rank(gene_lineage)
            gene_ranks2lineage = dict((rank, taxid) for (taxid, rank) in gene_lineage2ranks.items())
            gene_taxonomy_alignment = gene_ranks2lineage
            #gene_superkingdom = gene_taxonomy_alignment['superkingdom']
            #gene_kingdom = gene_taxonomy_alignment['kingdom']
            #gene_phylum = gene_taxonomy_alignment['phylum']
            #gene_subphylum = gene_taxonomy_alignment['subphylum']
            #gene_class = gene_taxonomy_alignment['class']
            #gene_order = gene_taxonomy_alignment['order']
            #gene_family = gene_taxonomy_alignment['family']
            #gene_genus = gene_taxonomy_alignment['genus']
            #gene_species = gene_taxonomy_alignment['species']
            gene_taxlevel = gene_taxonomy_alignment[tax_level]
            # print(gene_kingdom)
            # print(type(gene_kingdom)) # <class 'int'>
            # print(gene_subphylum)
        except :
            print('Attention: please check the gene accession id!')
            break

            
        # Get the taxonomic information for each accession in the blast file
        
        for accession in accession_number[:200] :
            try :
                taxid = getTaxid(accession)
                print('Taxid by BLAST:', taxid)
            except Exception as e:
                print(Exception, e)
                continue

            # Get the lineage of the taxid
            
            try :
                lineage = ncbi.get_lineage(taxid)
                lineage2ranks = ncbi.get_rank(lineage)
                ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
                #print(ranks2lineage)

                taxid2name = ncbi.get_taxid_translator(lineage)
                taxonomy_alignment = ranks2lineage
            except :
                print('Warning: %s taxid not found!' % str(taxid))
                continue


            # Check if the taxonomic level of the gene is the same as the taxonomic level of the accession
            
            try :
                #if taxonomy_alignment['order'] == gene_order : #replace first with subphylum
                #    recipient_accession.append(accession)
                #    recipient_species.append(taxonomy_alignment['species'])
                
                
                #if taxonomy_alignment['phylum'] == gene_phylum and taxonomy_alignment['order'] != gene_order : #kingdom and subphylum
                #    outgroup_accession.append(accession)
                #    outgroup_species.append(taxonomy_alignment['species'])
                
                # This section should work with the taxonomic level selected by the user, need to incorporate at a later date
                
                if taxonomy_alignment[tax_level] != gene_taxonomy_alignment[tax_level] :
                    outgroup_accession.append(accession)
                    outgroup_species.append(taxonomy_alignment['species'])
                if taxonomy_alignment[tax_level] == gene_taxonomy_alignment[tax_level] :
                    recipient_accession.append(accession)
                    recipient_species.append(taxonomy_alignment['species'])
                
            #except if taxonomy_alignment[tax_level] == null :
                #print('Warning: %s taxid not found!' % str(taxid))
                #continue        
            except Exception as e:
                print(Exception, e)
                continue

        # Get the bitscore of the recipient and outgroup accessions, and append to the dictionary

        for accession_id in recipient_accession :
            recipient_accession_bitscore[accession_id] = accession_bitscore[accession_id]
            #print('Recipient accession:', recipient_accession_bitscore)
        for accession_id in outgroup_accession :
            outgroup_accession_bitscore[accession_id] = accession_bitscore[accession_id]
            #print('Outgroup accession:', outgroup_accession_bitscore)
        if recipient_accession_bitscore :
            max_recipient_organism_accession_key = max(recipient_accession_bitscore,key=recipient_accession_bitscore.get)
            max_recipient_organism_bitscore = recipient_accession_bitscore[max_recipient_organism_accession_key]
            #print('Recipient accession:', max_recipient_organism_accession_key, max_recipient_organism_bitscore)
        print("Outgroups", outgroup_accession_bitscore)
        if outgroup_accession_bitscore :
            max_outgroup_accession_key = max(outgroup_accession_bitscore,key=outgroup_accession_bitscore.get)
            print(max_outgroup_accession_key)
            max_outgroup_bitscore = outgroup_accession_bitscore[max_outgroup_accession_key]
            #print('Outgroup accession:', max_outgroup_accession_key, max_outgroup_bitscore)
            if max_outgroup_accession_key :
                max_taxid = getTaxid(max_outgroup_accession_key)
                max_lineage = ncbi.get_lineage(max_taxid)
                max_lineage2ranks = ncbi.get_rank(max_lineage)
                max_ranks2lineage = dict((rank, taxid) for (taxid, rank) in max_lineage2ranks.items())
                try :
                    #max_taxid2name = ncbi.get_taxid_translator([max_ranks2lineage['superkingdom'], max_ranks2lineage['phylum'], max_ranks2lineage['order'], max_ranks2lineage['species']])
                    max_taxid2name = ncbi.get_taxid_translator([max_ranks2lineage[tax_level]])
                    print(max_taxid2name)
                except Exception as e:
                    print(Exception, e)
                    continue

                print(gene)
                print(max_recipient_organism_bitscore)
                print(max_outgroup_bitscore)
                print(max_taxid2name)

                if recipient_species :
                    recipient_species_number = len(set(recipient_species))
                if outgroup_species :
                    outgroup_species_number = len(set(outgroup_species))

                HGT_index = format(max_outgroup_bitscore/max_recipient_organism_bitscore, '.4f')
                Outg_pct = format(outgroup_species_number/(outgroup_species_number+recipient_species_number), '.4f')

                print('HGT index: %s' % str(HGT_index))
                print('Out_pct: %s' % str(Outg_pct))
                if max_outgroup_bitscore>=bitscore_parameter and float(HGT_index)>=HGTIndex and float(Outg_pct)>=out_pct :
                    print('This is a HGT event')
                    #taxonomy = max_taxid2name[max_ranks2lineage['phylum']] + '/' + max_taxid2name[max_ranks2lineage['class']]
                    taxonomy = max_taxid2name[max_ranks2lineage[tax_level]]
                    item = [gene, max_outgroup_bitscore, Outg_pct, HGT_index, taxonomy]
                    HGT.append(item)
                else :
                    print('This is not a HGT event')
                    item = [gene, max_outgroup_bitscore, Outg_pct, HGT_index, 'No']
                    HGT.append(item)

    outfile = open("./output_close_HGT.tsv", "wt")
    tsv_writer = csv.writer(outfile, delimiter="\t")
    column = ['Gene/Protein', 'Bitscore', 'Out_pct', 'HGT index', 'Donor taxonomy']
    tsv_writer.writerow(column)
    for HGT_info in HGT :
        tsv_writer.writerow(HGT_info)
    outfile.close()



if __name__== "__main__":
    
    main()
