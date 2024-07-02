#!//ibers/ernie/home/jac180/.conda/envs/HGTPhylo/bin/python
import sys, os, warnings, math, csv, argparse, subprocess
import pandas as pd
from Bio import SeqIO, BiopythonWarning, Entrez
from ete3 import NCBITaxa
#from taxoniq import Taxonomy


# Path to databases need to be specified
# Not terribly sure how I go about doing this for the time being
# Perhaps just leave it to the user to specify the database paths??


def parse_arguments():
    parser = argparse.ArgumentParser(description="Modified version of HGTPhyloDetect close workflow for HGT events, takes protein fasta file and iterates through each sequence outputting a likelihood of HGT origin for each", epilog="Author: Jack A. Crosby, Aberystwyth University/Queens University Belfast")
    parser.add_argument("input_file", help="Input file path, should be a fasta file of protein sequences")
    parser.add_argument("--bitscore_parameter", type=float, default=100, help="Bitscore parameter, default is 100")
    parser.add_argument("--HGTIndex", type=float, default=0.5, help="HGT Index, default is 0.5")
    parser.add_argument("--out_pct", type=float, default=0.8, help="Out Pct, default is 0.8")
    parser.add_argument("-t", "--tax_level", type=str, default="family", choices=["superkingdom", "kingdom", "phylum", "subphylum", "class", "order", "family", "genus", "species"], help="Taxonomic level, organisms outisde of this level will be classified as 'outgroup', default is family.")
    parser.add_argument("-s", "--search", type=str, default="blastp", choices=["blastp", "diamond", "mmseqs"], help="Search methods, blastp uses remote nr search, diamond & mmseqs use local database for search, default is blastp.")
    parser.add_argument("-u", "--update", action="store_true", help="Update the NCBI taxonomy database")
    parser.add_argument("-q", "--query_tax", type=int, help="Taxid associated with the query sequence")
    return parser.parse_args()



def search_protein_comb(search, name): #gene geneseq
    
    match search:
    
        case "blastp":
    
            #if os.path.exists("./blastp_files/%s.txt" % gene):
            #    accession_number, accession_bitscore = parse_NCBI(gene)
            #    print('Blast file exists')
            #else:
            #    with open('./%s.fasta' % gene, 'w') as outfile:
            #        outfile.write('>' + gene + '\n')
            #        outfile.write(geneSeq[gene] + '\n')
            #    if not os.path.exists("./blastp_files/"):
            #        os.system('mkdir ./blastp_files')
            #    myCmd = "blastp -db nr -remote -query ./%s.fasta -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue bitscore length pident' -out ./blastp_files/%s.txt" % (gene, gene)
            #    os.system(myCmd)
            #    os.remove('./%s.fasta' % gene)
            #   accession_number, accession_bitscore = parse_NCBI(gene)
            print("Invalid")
            
        case "diamond":
            
            if os.path.exists(f"{os.path.splitext(name)[0]}.tsv"):
                #accession_number, accession_bitscore = parse_NCBI(gene)
                print(f'Diamond file found for {os.path.splitext(name)[0]}')
            else:
                
                outf = str(name.split(".")[0] + ".tsv")
                myCmd =f'diamond blastp -d /fast-scratch/jac180/nr/nr_db.dmnd -q {name} --max-target-seqs 250 --outfmt 6 qseqid sseqid evalue bitscore length pident staxids -o {outf}'
                myCmd = str(myCmd)
                os.system(myCmd)
                my1Cmd = [
                    "diamond", "blastp",
                    "-d", "/fast-scratch/jac180/nr/nr_db.dmnd",
                    "-q", name,
                    "--max-target-seqs", "250",
                    "-outfmt", "6 qseqid sseqid evalue bitscore length pident staxids",
                    "--out", outf
                ]

                #result = subprocess.run(myCmd, capture_output=True, text=True)

                # Check if the command was successful
                #if result.returncode == 0:
                #    print("Command executed successfully!")
                #else:
                #    print("Error:", result.stderr)
                #    sys.exit()
                #os.remove('./%s.fasta' % gene)
                #accession_number, accession_bitscore = parse_NCBI(gene)
        

            

        case "mmseqs":
          
            if os.path.exists(f'MMSeqs file found for {os.path.splitext(name)[0]}.tsv'):
                print(f'MMSeqs file found for {os.path.splitext(name)[0]}.tsv')
            else:
                outf = str(name.split(".")[0] + ".tsv")

                myCmd = f"mmseqs search {name} /path/to/nr_database {name}.tsv --format-output 'query, target, evalue, bitscore, length, pident, taxid'"
                myCmd = str(myCmd)
                os.system(myCmd)
            
        case _:
            
            print(f"Unknown search method: {search}")
            print("Please choose from blastp, diamond or mmseqs")
            print("Exiting...")
            sys.exit(1)
    
    #return accession_number, accession_bitscore



def load_diamond_results(combined_file, gene):
    # Load diamond results file into dataframe
    results=pd.read_csv(combined_file, sep='\t', header=None)
    # Filter the results for the gene of interest
    gene_results = results[results[0] == gene]
    # Return filtered results
    return gene_results



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



def get_taxid(gene_results, accession_number):
    df = pd.read_csv(gene_results, sep='\t', header=None)
    filtered_results = df[df[1] == accession_number]
    taxid = filtered_results[6].str.split(';').str[-1].values[0]
    
    # Filter the gene results for the specific accession number
    #filtered_results = gene_results[gene_results[0] == accession_number]
    # Extract the taxid from the filtered results
    #taxid = filtered_results[6].str.split(';').str[-1].values[0]
    return taxid



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

def process_gene(gene_data):
    gene, sequence, ncbi, tax_level, bitscore_parameter, HGTIndex, out_pct, result_file, qtaxid = gene_data
    
    print(f'Processing gene: {gene}')
    
    # Check for the presence of results files
    if os.path.exists(f"./blastp_files/{gene}.txt"):
        accession_number, accession_bitscore = parse_NCBI(gene)
        print('Blast file exists')
    elif os.path.exists(result_file):
        gene_results = load_diamond_results(result_file, gene)
        accession_number = gene_results[1]
        accession_bitscore = gene_results[3]
    else:
        print('No results file found')
        return None

    # Get taxids for the gene
    try:
        gene_taxid = qtaxid
        gene_lineage = ncbi.get_lineage(gene_taxid)
        gene_lineage2ranks = ncbi.get_rank(gene_lineage)
        gene_ranks2lineage = dict((rank, taxid) for (taxid, rank) in gene_lineage2ranks.items())
        gene_taxonomy_alignment = gene_ranks2lineage
        gene_taxlevel = gene_taxonomy_alignment[tax_level]
    except Exception as e:
        print(f"Error processing gene taxonomy: {e}")
        return None

    # Process accessions
    recipient_accession = []
    outgroup_accession = []
    recipient_species = []
    outgroup_species = []

    tax_alignments = taxonomy_alignment
    for accession in accession_number[:200]:
        try:
            taxid = get_taxid(result_file, accession)
            lineage = ncbi.get_lineage(taxid)
            lineage2ranks = ncbi.get_rank(lineage)
            ranks2lineage = dict((rank,taxid) for (taxid,rank) in lineage2ranks.items())
            taxonomy_alignment = ranks2lineage

            if taxonomy_alignment[tax_level] == gene_taxonomy_alignment[tax_level]:
                recipient_accession.append(accession)
                recipient_species.append(taxonomy_alignment['species'])
            else:
                outgroup_accession.append(accession)
                outgroup_species.append(taxonomy_alignment['species'])
        except Exception as e:
            continue

    # Process bitscores
    recipient_accession_bitscore = {acc: score for acc, score in zip(accession_number, accession_bitscore) if acc in recipient_accession}
    outgroup_accession_bitscore = {acc: score for acc, score in zip(accession_number, accession_bitscore) if acc in outgroup_accession}

    if not recipient_accession_bitscore or not outgroup_accession_bitscore:
        return None

    max_recipient_organism_accession_key = max(recipient_accession_bitscore, key=recipient_accession_bitscore.get)
    max_recipient_organism_bitscore = recipient_accession_bitscore[max_recipient_organism_accession_key]
    max_outgroup_accession_key = max(outgroup_accession_bitscore, key=outgroup_accession_bitscore.get)
    max_outgroup_bitscore = outgroup_accession_bitscore[max_outgroup_accession_key]

    # Calculate HGT index and Outgroup percentage
    recipient_species_number = len(set(recipient_species))
    outgroup_species_number = len(set(outgroup_species))
    HGT_index = format(max_outgroup_bitscore/max_recipient_organism_bitscore, '.4f')
    Outg_pct = format(outgroup_species_number/(outgroup_species_number+recipient_species_number), '.4f')

    # Determine if it's a HGT event
    if max_outgroup_bitscore >= bitscore_parameter and float(HGT_index) >= HGTIndex and float(Outg_pct) >= out_pct:
        max_taxid = get_taxid(result_file, max_outgroup_accession_key)
        max_lineage = ncbi.get_lineage(max_taxid)
        max_lineage2ranks = ncbi.get_rank(max_lineage)
        max_ranks2lineage = dict((rank, taxid) for (taxid, rank) in max_lineage2ranks.items())
        max_taxid2name = ncbi.get_taxid_translator([max_ranks2lineage[tax_level]])
        taxonomy = max_taxid2name[max_ranks2lineage[tax_level]]
        return [gene, max_outgroup_bitscore, Outg_pct, HGT_index, taxonomy]
    else:
        return [gene, max_outgroup_bitscore, Outg_pct, HGT_index, 'No']

def batch_fetch_taxonomy(taxids):
    # Fetch lineages for all taxids
    taxid_to_lineage = {taxid: ncbi.get_lineage(taxid) for taxid in taxids}
    
    # Fetch ranks for all unique taxids in all lineages
    unique_taxids = set(taxid for lineage in taxid_to_lineage.values() for taxid in lineage)
    taxid_to_rank = ncbi.get_rank(unique_taxids)
    
    # Convert lineages to rank-to-taxid mappings
    taxid_to_taxonomy_alignment = {}
    for taxid, lineage in taxid_to_lineage.items():
        ranks2lineage = {taxid_to_rank[t]: t for t in lineage if t in taxid_to_rank}
        taxid_to_taxonomy_alignment[taxid] = ranks2lineage
    
    return taxid_to_taxonomy_alignment

def main():
    """
    Runs HGTPhyloDetect but slightly modified, takes arguments regarding the
    taxonomic level to determine the outgroup and the search method to use.
        
    """
    
    # Initialise args, NCBI database and print used parameters
    
    parse_arguments()
    name = parse_arguments().input_file
    bitscore_parameter = parse_arguments().bitscore_parameter
    HGTIndex = parse_arguments().HGTIndex
    out_pct = parse_arguments().out_pct
    tax_level = parse_arguments().tax_level.lower()
    search = parse_arguments().search.lower()
    update = parse_arguments().update
    qtaxid = parse_arguments().query_tax
    ncbi = NCBITaxa()
    if update:
        ncbi.update_taxonomy_database()
    genes = list()
    geneSeq = dict()
    HGT = list()
    warnings.simplefilter('ignore', BiopythonWarning)
    # Print the table header
    print("Input Parameters:")
    print("-----------------")
    print(f"{'Input File':<20} | {name}")
    print(f"{'Bitscore Parameter':<20} | {bitscore_parameter}")
    print(f"{'HGT Index':<20} | {HGTIndex}")
    print(f"{'Outgroup Percentage':<20} | {out_pct}")
    print(f"{'Taxonomic Level':<20} | {tax_level}")
    print(f"{'Search Method':<20} | {search}")
    print("-----------------")
    
    # Read in the gene sequences from the input file
    # Store the gene sequences in a dictionary
    # Store the gene names in a list
    
    with open(name, 'r') as handleGene:
        for record in SeqIO.parse(handleGene, "fasta"):
            gene = str(record.id)
            sequence = str(record.seq)
            geneSeq[gene] = sequence
            genes.append(gene)
    
    n=0
    search_protein_comb(search,name)
    result_file=str(f"{os.path.splitext(name)[0]}.tsv")
    for gene in genes:
        n+=1
        print(f'This is gene {n}')
        print(f'Gene: {gene}')
        
        # Check for the presnece of results files
        # If they exist, parse the results
        # If they don't exist, run the search
        if os.path.exists("./blastp_files/%s.txt" % gene):
            accession_number, accession_bitscore = parse_NCBI(gene)
            print('Blast file exists')
        elif os.path.exists(result_file):
            gene_results = load_diamond_results(result_file, gene)
            accession_number = gene_results[1]
            accession_bitscore = gene_results[3]
        else:
            print('No results file found')
            break
    
        # Initialise ncbi database
        # Initialise lists and dictionaries relating to ingroup and outgroups
        #ncbi = NCBITaxa(dbfile='~/.etetoolkit/taxa.sqlite')
        ncbi = NCBITaxa()
        recipient_accession = list()
        recipient_accession = list()
        outgroup_accession = list()
        recipient_accession_bitscore = dict()
        outgroup_accession_bitscore = dict()
        recipient_species = list()
        outgroup_species = list()

        
        # Get taxids for each gene in the results file
        # Get the taxonomic lineage for each gene
        # Get the taxonomic level of each gene
        
        try :
            # Pull taxonomic information from the gene accessions 
            
            #gene_taxid = getTaxid(gene)
            gene_taxid = qtaxid
            #gene_taxid =get_taxid(result_file,gene)
            gene_lineage = ncbi.get_lineage(gene_taxid)
            gene_lineage2ranks = ncbi.get_rank(gene_lineage)
            gene_ranks2lineage = dict((rank, taxid) for (taxid, rank) in gene_lineage2ranks.items())
            gene_taxonomy_alignment = gene_ranks2lineage
            gene_taxlevel = gene_taxonomy_alignment[tax_level]
            print("Gene Taxonomy Information:")
            print("--------------------------")
            for rank, taxid in gene_taxonomy_alignment.items():
                print(f"{rank.capitalize():<20} | {taxid}")
            print("--------------------------")
        except Exception as e:
            print(f"Error type: {e.__class__.__name__}, Message: {e}")
            print(gene_taxonomy_alignment)
            break

       
        # Get the taxonomic lineage for each gene in the results file
        #print(accession_number)
        #print(accession_bitscore)
        taxids = [9606, 10090]  # Human and Mouse taxids (Had these to hand)
         
        taxonomy_alignments = batch_fetch_taxonomy(taxids)
        print(taxonomy_alignments)

        for accession in accession_number[:200] :
            try :
                taxid = get_taxid(result_file, accession)
                #print('NCBI taxid:', taxid)
            except Exception as e:
                print(f"Error type: {e.__class__.__name__}, Message: {e}")
                print('line 286 fail')
                break
            try :
                lineage = ncbi.get_lineage(taxid)
                lineage2ranks = ncbi.get_rank(lineage)
                ranks2lineage = dict((rank,taxid) for (taxid,rank) in lineage2ranks.items())
                taxonomy_alignment = ranks2lineage
                #print(taxonomy_alignment, 'query tax alignment')
                #taxlevel = taxonomy_alignment[tax_level]
            except :
                #print(f'Error: {accession} not found in NCBI taxonomy database')
                continue
            #print('Line 297')
            try :
                if taxonomy_alignment[tax_level] == gene_taxonomy_alignment[tax_level]:
                    recipient_accession.append(accession)
                    recipient_species.append(taxonomy_alignment['species'])
                if taxonomy_alignment[tax_level] != gene_taxonomy_alignment[tax_level]:
                    outgroup_accession.append(accession)
                    outgroup_species.append(taxonomy_alignment['species'])
            except Exception as e:
                #print(Exception, e, 'line 306')
                continue
        
        
        recipient_accession_bitscore = {}
        for accession_id, bitscore in zip(accession_number, accession_bitscore):
            if accession_id in recipient_accession:
                recipient_accession_bitscore[accession_id] = bitscore

        outgroup_accession_bitscore = {}
        for accession_id, bitscore in zip(accession_number, accession_bitscore):
            if accession_id in outgroup_accession:
                outgroup_accession_bitscore[accession_id] = bitscore
        #for accession_id in recipient_accession :
            #recipient_accession_bitscore[accession_id] = accession_bitscore[accession_id]
       #     recipient_accession_bitscore = dict(zip(accession_number, accession_bitscore))
       # for accession_id in outgroup_accession :
            #outgroup_accession_bitscore[accession_id] = accession_bitscore[accession_id]
        #    outgroup_accession_bitscore = dict(zip(accession_number, accession_bitscore))
        if recipient_accession_bitscore :
            max_recipient_organism_accession_key = max(recipient_accession_bitscore,key=recipient_accession_bitscore.get)
            max_recipient_organism_bitscore = recipient_accession_bitscore[max_recipient_organism_accession_key]
        if outgroup_accession_bitscore :
            max_outgroup_accession_key = max(outgroup_accession_bitscore,key=outgroup_accession_bitscore.get)
            max_outgroup_bitscore = outgroup_accession_bitscore[max_outgroup_accession_key]
            if max_outgroup_accession_key :
                max_taxid = get_taxid(result_file, max_outgroup_accession_key)
                max_lineage = ncbi.get_lineage(max_taxid)
                max_lineage2ranks = ncbi.get_rank(max_lineage)
                max_ranks2lineage = dict((rank, taxid) for (taxid, rank) in max_lineage2ranks.items())
                try :
                    max_taxid2name = ncbi.get_taxid_translator([max_ranks2lineage[tax_level]])
                except Exception as e:
                    #print(Exception, e)
                    continue
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

    outfile = open(f"./output_{tax_level}_HGT.tsv", "wt")
    tsv_writer = csv.writer(outfile, delimiter="\t")
    column = ['Gene/Protein', 'Bitscore', 'Out_pct', 'HGT index', 'Donor taxonomy']
    tsv_writer.writerow(column)
    for HGT_info in HGT :
        tsv_writer.writerow(HGT_info)
    outfile.close()



if __name__== "__main__":
    
    main()                