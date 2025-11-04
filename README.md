# PanelPub
Leveraging PubTator 3.0 (https://www.ncbi.nlm.nih.gov/research/pubtator3/) and LitVar 2.0 (https://www.ncbi.nlm.nih.gov/research/litvar2/) to retrieve genetic annotations from PubMed articles and compare them to virtual gene panels of PanelApp (https://panelapp-aus.org/ and https://panelapp.genomicsengland.co.uk/).

<p align="center">
  <img width="750" height="563" alt="gene-panel-PubTator_v2" src="https://github.com/user-attachments/assets/1b841249-67e9-454f-ab96-40c3ed798432" />
</p>

## Data sources
### Gene synonyms
1. Download gene synonyms from https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
2. On the command line and in the directory where the gene_info.gz file is located, run the following to retrieve human genes; gene synonyms will be '|'-separated in the column 'Synonyms'

```bash
zcat gene_info.gz | grep -e "^#tax_id" -e "^9606" | awk '{print $2, $3, $5}' > gene_synonyms_hsapiens.txt
```
3. Alternatively, download on the main page: gene_synonyms_hsapiens.txt

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

### GO-BP gene sets
Download Gene Ontology - Biological Process data from https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Hs/c5.all.v2025.1.Hs.symbols.gmt
