# PanelPub
Leveraging PubTator 3.0 (https://www.ncbi.nlm.nih.gov/research/pubtator3/) and LitVar 2.0 (https://www.ncbi.nlm.nih.gov/research/litvar2/) to retrieve genetic annotations from PubMed articles and compare them to virtual gene panels of PanelApp (https://panelapp-aus.org/ and https://panelapp.genomicsengland.co.uk/).

<p align="center">
  <img width="750" height="563" alt="gene-panel-PubTator_v2" src="https://github.com/user-attachments/assets/84ad1544-1415-4c88-86a0-76c98aff4558" />
</p>

## Data sources

### PubTator 3.0
Download PT3 genes, variants, and relations from:

<b>https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/gene2pubtator3.gz</b>

<b>https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/mutation2pubtator3.gz</b>

<b>https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/relation2pubtator3.gz</b>

### LitVar 2.0
Download LitVar 2.0 data from: 

<b>https://ftp.ncbi.nlm.nih.gov/pub/lu/LitVar/litvar2_variants.json.gz</b>

### GO-BP gene sets
Download Gene Ontology - Biological Process data from: 

<b>https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Hs/c5.go.bp.v2025.1.Hs.entrez.gmt</b>

### HGNC-approved symbols
Download HUGO Gene Nomenclature Committee (HGNC) Approved Symbols from: 

<b>https://www.genenames.org https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt</b>

### MeSH data
Download MeSH descriptors and Supplementary Concepts from: 

<b>https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2025.gz</b>

<b>https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/supp2025.gz</b>

### PMIDs
1. Make sure eutils is downloade

```bash
sudo apt install ncbi-entrez-direct
```

2. On the command line, run the following, where "query" is your PubMed query 

```bash
esearch -db pubmed -query "query" | efetch -format uilist > file.txt
```

<b>Example</b>
```
esearch -db pubmed -query "Hypogonadism/genetics[mesh] OR 'hypogonad* hypogonad*'[tiab]" | efetch -format uilist > CHH.txt
```

### Non-human orthologs
1. Download ortholog data from https://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz
2. In R, run the following to have human genes in one column (hs_gene) and non-human ortholog info in the others (ortholog_gene & ortholog_tax)

```R
library(data.table)
library(dplyr)
ortholog_df <- fread("path/to/gene_orthologs")
ortho <- ortholog_df %>%
  rename(tax_id = `#tax_id`) %>%
  filter(tax_id == 9606 | Other_tax_id == 9606) %>%
  mutate(hs_tax = if_else(tax_id == 9606, tax_id, Other_tax_id),
         hs_gene = if_else(tax_id == 9606, GeneID, Other_GeneID),
         ortholog_gene = if_else(tax_id == 9606, Other_GeneID, GeneID),
         ortholog_tax = if_else(tax_id == 9606, Other_tax_id, tax_id)) %>%
  select(hs_gene, ortholog_tax, ortholog_gene)
fwrite(ortho, "path/to/gene_orthologs_MOD.txt", sep = "\t")
```

### Gene synonyms
1. Download gene synonyms from https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
2. In the directory where gene_info.gz is located, run the following python script in the command line to retrieve human genes; gene synonyms will be '|'-separated in the column 'Synonyms'

```bash
python3 filter_human_gene.py -g gene_info.gz
```
Alternatively, download on the main page: <b>gene_info_human.txt</b>
