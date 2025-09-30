# PubTator
Leveraging PubTator 3.0 (https://www.ncbi.nlm.nih.gov/research/pubtator3/) and LitVar 2.0 (https://www.ncbi.nlm.nih.gov/research/litvar2/) to retrieve genetic annotations from PubMed articles and comparing them to virtual gene panels of PanelApp (https://panelapp-aus.org/ and https://panelapp.genomicsengland.co.uk/).

## Data sources
### Gene synonyms
1. Go to https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/9606/
2. Select columns: 'Gene ID', 'Symbol', 'Synonyms'
3. Check the box next to 'Gene ID'; all rows should be checked
4. Click Download, Download Table, then Download

### Non-human orthologs
1. Download ortholog data from https://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz
