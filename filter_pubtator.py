import requests
import json
import argparse
import pandas as pd
import os
from collections import Counter

def filter_pubtator():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--pubtator_file", type=str)
    parser.add_argument("-p", "--pmid", type=str)
    parser.add_argument("-c", "--chunk_size", type=int, default=10000)
    parser.add_argument("-o", "--output", type=str, default = "PTC_filter.txt")

    args = parser.parse_args()
    
    found_counts = Counter()
    
    pubtator_file = pd.read_csv(
        args.pubtator_file,
        compression="gzip",
        sep="\t",
        chunksize=args.chunk_size,
        header=None,
        names=['PMID', 'Type', 'ID', 'Name', 'Norm'],
        dtype=str,
        keep_default_na=False, 
        na_filter=False
    )#, on_bad_lines='skip')
    # pmid_file = pd.read_csv(args.pmid, sep="\t", header=None)
    # pmids = list(set(pmid_file.iloc[:,0].tolist()))
    with open(args.pmid, "r") as f:
        pmids = {line.strip() for line in f if line.strip()}
        # pmids = sorted({int(x) for x in pmids}, reverse=True)
    
    output = os.path.join(os.getcwd(), args.output)
    first = True
    for chunk in pubtator_file:
        chunk['PMID'] = chunk['PMID'].astype(str).str.strip()
        chunk_pmid = chunk[chunk['PMID'].isin(pmids)]
        
        if first:
            chunk_pmid.to_csv(output, sep="\t", index=False, mode="w")
            first = False
        else:
            chunk_pmid.to_csv(output, sep="\t", index=False, header=False, mode="a")
        
        found_counts.update(chunk_pmid['PMID'].tolist())

    missing_ids = sorted(pmids - set(found_counts.keys()))
    
    missing_output = output.replace(".txt", "_missing.txt")
    with open(missing_output, "w", encoding="utf-8") as f:
        f.write("\n".join(missing_ids))

filter_pubtator()
        