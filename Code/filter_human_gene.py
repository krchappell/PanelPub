import argparse
import pandas as pd
import os

def filter_human_gene():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gene_info", type=str)
    parser.add_argument("-c", "--chunk_size", type=int, default=10000)
    parser.add_argument("-o", "--output", type=str, default = "gene_info_human.txt")

    # #tax_id, GeneID, Symbol, LocusTag, Synonyms, dbXrefs, chromosome, map_location, description, type_of_gene, Symbol_from_nomenclature_authority, 
    # Full_name_from_nomenclature_authority, Nomenclature_status, Other_designations, Modification_date, Feature_type
    
    args = parser.parse_args()
        
    gene_info_file = pd.read_csv(args.gene_info, compression="gzip", sep="\t", chunksize=args.chunk_size)

    output = os.path.join(os.getcwd(), args.output)
    first = True
    for chunk in gene_info_file:
        chunk_hs = chunk[chunk['#tax_id']==9606]
        
        if first:
            chunk_hs.to_csv(output, sep="\t", index=False, mode="w")
            first = False
        else:
            chunk_hs.to_csv(output, sep="\t", index=False, header=False, mode="a")

filter_human_gene()
        
