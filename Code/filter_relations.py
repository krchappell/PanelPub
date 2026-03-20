import requests
import json
import argparse
import pandas as pd
import os

def filter_relations():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--relations", type=str)
    parser.add_argument("-m", "--mesh", type=str)
    parser.add_argument("-c", "--chunk_size", type=int, default=10000)
    parser.add_argument("-o", "--output", type=str, default = "Relations_filter.txt")

    args = parser.parse_args()
    relations = pd.read_csv(args.relations, sep="\t", chunksize=args.chunk_size, header=None, names=['PMID', 'Type', 'Concept1', 'Concept2'], on_bad_lines='skip')

    mesh_filter = pd.read_csv(args.mesh, sep="\t", header=None)
    mesh = list(set(mesh_filter.iloc[:,0].tolist()))

    output = os.path.join(os.getcwd(), args.output)
    first = True
    for chunk in relations:
        chunk[['Type1', 'ID1']] = chunk['Concept1'].str.split('|', expand=True)
        chunk['ID1'] = chunk['ID1'].replace("MESH:", "", regex=True) 
        chunk[['Type2', 'ID2']] = chunk['Concept2'].str.split('|', expand=True)
        chunk['ID2'] = chunk['ID2'].replace("MESH:", "", regex=True) 
        chunk_mesh1 = chunk[chunk['ID1'].isin(mesh)]
        chunk_mesh2 = chunk[chunk['ID2'].isin(mesh)]
        chunk_mesh = pd.concat([chunk_mesh1, chunk_mesh2], ignore_index=True)
        if first:
            chunk_mesh.to_csv(output, sep="\t", index=False, mode="w")
            first = False
        else:
            chunk_mesh.to_csv(output, sep="\t", index=False, header=False, mode="a")
        
filter_relations()
