# PubTator
Leveraging PubTator 3.0 (https://www.ncbi.nlm.nih.gov/research/pubtator3/) and LitVar 2.0 (https://www.ncbi.nlm.nih.gov/research/litvar2/) to retrieve genetic annotations from PubMed articles and comparing them to virtual gene panels of PanelApp (https://panelapp-aus.org/ and https://panelapp.genomicsengland.co.uk/).

## Data sources
### Gene synonyms
1. Go to https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/9606/
2. Select columns 'Gene ID', 'Symbol', 'Synonyms'
3. Check the box next to 'Gene ID'; all rows should be checked
4. Click Download, Download Table, then Download

### Non-human orthologs
1. Download ortholog data from https://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz
2. In R, run the following to have human genes in one column (hs_gene) and non-human ortholog info in the others (ortholog_gene & ortholog_tax)

```R
ortholog_df <- fread("C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/gene_orthologs")
ortho <- ortholog_df %>%
  rename(tax_id = `#tax_id`) %>%
  filter(tax_id == 9606 | Other_tax_id == 9606) %>%
  mutate(hs_tax = if_else(tax_id == 9606, tax_id, Other_tax_id),
         hs_gene = if_else(tax_id == 9606, GeneID, Other_GeneID),
         ortholog_gene = if_else(tax_id == 9606, Other_GeneID, GeneID),
         ortholog_tax = if_else(tax_id == 9606, Other_tax_id, tax_id)) %>%
  select(hs_gene, ortholog_tax, ortholog_gene)
fwrite(ortho, "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/gene_orthologs_MOD.txt", sep = "\t")
```

### GO-BP gene sets
Download Gene Ontology - Biological Process data from https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Hs/c5.all.v2025.1.Hs.symbols.gmt
