import argparse
import pandas as pd
import os
from collections import Counter

def filter_litvar():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--litvar", type=str)
    parser.add_argument("-r", "--rs", type=str)
    parser.add_argument("-c", "--chunk_size", type=int, default=10000)
    parser.add_argument("-o", "--output", type=str, default = "LitVar2_filter.txt")

    args = parser.parse_args()
    
    found_counts = Counter()
    
    litvar_file = pd.read_csv(args.litvar, compression="gzip", sep="\t", chunksize=args.chunk_size)
    with open(args.rs, "r") as f:
        rs_nums = {line.strip() for line in f if line.strip()}
    
    output = os.path.join(os.getcwd(), args.output)
    first = True
    for chunk in litvar_file:
        chunk['rsid'] = chunk['rsid'].astype(str).str.strip()
        chunk_rs = chunk[chunk['rsid'].isin(rs_nums)]
        
        if first:
            chunk_rs.to_csv(output, sep="\t", index=False, mode="w")
            first = False
        else:
            chunk_rs.to_csv(output, sep="\t", index=False, header=False, mode="a")
            
    missing_ids = sorted(rs_nums - set(found_counts.keys()))
    
    missing_output = output.replace(".txt", "_missing.txt")
    with open(missing_output, "w", encoding="utf-8") as f:
        f.write("\n".join(missing_ids))

filter_litvar()
        
