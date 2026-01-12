# PanelPub Data Extraction ####

## load dependencies ####
library(data.table)
library(dplyr)
library(fgsea)
library(ggplot2)
library(ggpubr)
library(ggvenn)
library(httr)
library(jsonlite)
library(patchwork)
library(pubmed.mineR)
library(stringr)
library(tidytext)
library(tidyverse)
library(xlsx)

## NCBI orthologs ####

# ortholog data from https://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz
ortholog_df <- fread("path/to/gene_orthologs.gz")

ortho <- ortholog_df %>%
  rename(tax_id = `#tax_id`) %>%
  filter(tax_id == 9606 | Other_tax_id == 9606) %>%
  mutate(hs_tax = if_else(tax_id == 9606, tax_id, Other_tax_id),
         hs_gene = if_else(tax_id == 9606, GeneID, Other_GeneID),
         ortholog_gene = if_else(tax_id == 9606, Other_GeneID, GeneID),
         ortholog_tax = if_else(tax_id == 9606, Other_tax_id, tax_id)) %>% 
  select(hs_gene, ortholog_tax, ortholog_gene)

fwrite(ortho, "path/to/gene_orthologs_MOD.txt", sep = "\t")

## PanelApp England & Australia ####

# England
panelapp_ids_eng <- PanelIDs_PanelApp("ENG")
panelapp_genes_eng <- lapply(panelapp_ids_eng, PanelGenes_PanelApp, panelapp = "ENG")
panelapp_genes_eng_df <- data.frame(rbindlist(panelapp_genes_eng, fill = TRUE))

panelapp_genes_eng_df_tabreplace <- panelapp_genes_eng_df %>%
  mutate_if(is.character, str_replace_all, "\t", " ")

panelapp_genes_eng_df_tabreplace <- panelapp_genes_eng_df_tabreplace %>%
  mutate(color = gsub("Expert Review ", "", 
                             str_extract(panelapp_genes_eng_df$evidence, pattern = "Expert Review [A-Za-z]+"))) %>%
  mutate(color = if_else(is.na(color), "Red", color),
         color = case_when(color == "green" ~ "Green",
                           color == "amber" ~ "Amber",
                           color == "red" ~ "Red"))

date_PA <- "04-11-2025"
fwrite(panelapp_genes_eng_df_tabreplace, sep = "\t",
       file = paste0("path/to/PanelAppEngland_panels_",
                     date_PA,
                     ".txt"))

# Australia
panelapp_ids_aus <- PanelIDs_PanelApp("AUS")
panelapp_genes_aus <- lapply(panelapp_ids_aus, PanelGenes_PanelApp, panelapp = "AUS")
panelapp_genes_aus_df <- data.frame(rbindlist(panelapp_genes_aus, fill = TRUE))

panelapp_genes_aus_df_tabreplace <- panelapp_genes_aus_df %>%
  mutate_if(is.character, str_replace_all, "\t", " ")

panelapp_genes_aus_df_tabreplace <- panelapp_genes_aus_df_tabreplace %>%
  mutate(color = gsub("Expert Review ", "",
                      str_extract(panelapp_genes_aus_df$evidence, pattern = "Expert Review [A-Za-z]+"))) %>%
  mutate(color = if_else(is.na(color), "Red", color),
         color = case_when(color == "green" ~ "Green",
                           color == "amber" ~ "Amber",
                           color == "red" ~ "Red"))

date_PA <- "04-11-2025"
fwrite(panelapp_genes_aus_df_tabreplace, sep = "\t",
       file = paste0("path/to/PanelAppAustralia_panels_", 
                     date_PA,
                     ".txt"))

# Combined
date_PA <- "04-11-2025"
panel_app_all <- panelapp_genes_eng_df_tabreplace %>%
  mutate(panelapp = "ENG",
         hgnc_release = as.character(hgnc_release)) %>%
  bind_rows(panelapp_genes_aus_df_tabreplace %>%
              mutate(panelapp = "AUS",
                     hgnc_release = as.character(hgnc_release)))
dim(panel_app_all)

fwrite(panel_app_all, sep = "\t",
       file = paste0("path/to/PanelAppENG-AUS_panels_", 
                     date_PA, 
                     ".txt"))

panel_app_clean <- panel_app_all %>%
  mutate(entry_id = 1:nrow(panel_app_all)) %>%
  separate_longer_delim(cols = publications, delim = "~") %>%
  separate_longer_delim(cols = publications, delim = ";")

nrow(panel_app_clean)
a <- nrow(panel_app_clean[grepl("[^0-9]", panel_app_clean$publications) & panel_app_clean$publications != "", ])
b <- nrow(panel_app_clean[grepl("^[0-9]+$", panel_app_clean$publications) & panel_app_clean$publications != "", ])
c <- nrow(panel_app_clean[panel_app_clean$publications == "", ])
sum(a, b, c)

# write non-PMID publications for manual PMID curation
panel_app_key <- panel_app_clean %>%
  mutate(publications2 = gsub("PMID: |PMID:|PMID|PMIDS: |PMIS: |PubMed: |PubMed:|PubMed ", "", publications),
                publications2 = gsub("^([0-9]+).*", "\\1", publications2),
                publications2 = gsub("[\\[(](?=\\d)", "", publications2, perl = TRUE),
                publications2 = gsub("(?<=\\d)[\\])]", "", publications2, perl = TRUE)) %>%
  mutate(publications2 = str_trim(publications2)) %>%
  filter(grepl("[^0-9]", publications2) & publications2 != "") %>%
  group_by(publications2) %>%
  summarize(entries = paste_unique(entry_id)) %>%
  arrange(publications2)

fwrite(panel_app_key, sep = "\t",
       file = "path/to/PanelApp_Key.txt")

# clean up PMIDs
panel_app_all_replace <- panel_app_clean %>%
  mutate(publications2 = gsub("PMID: |PMID:|PMID|PMIDS: |PMIS: |PubMed: |PubMed:|PubMed ", "", publications),
         publications2 = gsub("^([0-9]+).*", "\\1", publications2),
         publications2 = gsub("[\\[(](?=\\d)", "", publications2, perl = TRUE),
         publications2 = gsub("(?<=\\d)[\\])]", "", publications2, perl = TRUE)) %>%
  mutate(publications2 = str_trim(publications2))

# read in manually curated key and join PMIDs
key <- fread("path/to/PanelApp_Key_05-11-2025.txt", sep = "\t")
key_long <- key %>%
  separate_longer_delim(cols = "entries", delim = ";") %>%
  group_by(entries) %>%
  summarize(PMID = paste_unique(PMID)) %>%
  mutate(PMID = gsub("^NA;", "", PMID),
                PMID = gsub(";NA$", "", PMID),
                PMID = gsub(";NA;", "", PMID))

# fix PMIDs with errors (links from PanelApp pages)
panel_app_all_replace_join <- panel_app_all_replace %>%
  left_join(key_long %>%
                     separate_longer_delim(entries, delim = ";") %>%
                     mutate(entries = as.numeric(entries)),
                   by = join_by(entry_id == entries)) %>%
  mutate(publications2 = if_else(is.na(PMID), publications2, PMID)) %>%
  select(-PMID) %>%
  mutate(publications2 = if_else(publications2 == "None", NA, publications2),
                publications2 = str_replace(publications2, "\n", ";"),
                # 29721915 redirected to 30740741 (duplicate)
                # 31270415 is a abstract collection from the 51st European Society of Human Genetics Conference
                # 35234647 redirected to 35770779 (duplicate)
                publications2 = if_else(publications2 == "29721915", "30740741", publications2),
                publications2 = if_else(publications2 == "35234647", "35770779", publications2),
                
                publications2 = if_else(publications2 == "0587156", "30587156", publications2),
                publications2 = if_else(publications2 == "0714330", "30714330", publications2),
                publications2 = if_else(publications2 == "0932188", "10932188", publications2),
                publications2 = if_else(publications2 == "0961548", "30961548", publications2),
                publications2 = if_else(publications2 == "0980531", "10980531", publications2),
                
                publications2 = if_else(publications2 == "104693112", "10469312", publications2),
                publications2 = if_else(publications2 == "188806880", "18806880", publications2),
                publications2 = if_else(publications2 == "201788319", "20178831", publications2),
                publications2 = if_else(publications2 == "2194867118000911", "21948671;18000911", publications2),
                publications2 = if_else(publications2 == "239975631", "23997563", publications2),
                publications2 = if_else(publications2 == "240396609", "24039609", publications2),
                publications2 = if_else(publications2 == "251920460", "25192046", publications2),
                publications2 = if_else(publications2 == "2552958229174527", "25529582;29174527", publications2),
                publications2 = if_else(publications2 == "300068544", "30068544", publications2),
                publications2 = if_else(publications2 == "32232222962", "32222962", publications2),
                publications2 = if_else(publications2 == "328840387", "32884387", publications2),
                publications2 = if_else(publications2 == "329228291", "32928291", publications2),
                publications2 = if_else(publications2 == "0", "NA", publications2), 
                publications2 = if_else(publications2 == "NA" | publications2 == "", NA, publications2),
                publications2 = gsub("^NA;", "", publications2),
                publications2 = gsub(";NA$", "", publications2),
                publications2 = gsub(";NA;", "", publications2),
                color = gsub("Expert Review ", "", str_extract(panel_app_all_replace$evidence, pattern = "Expert Review [A-Za-z]+")),
                color = if_else(is.na(color), "Red", color)) %>%
  group_by(entry_id) %>%
  summarize(across(.cols = everything(), paste_unique)) %>%
  select(-entry_id)

panel_app_all_replace_join <- panel_app_all_replace_join %>%
  mutate(publications2 = gsub(";NA$", "", publications2))

dim(panel_app_all)
dim(panel_app_all_replace_join)

panel_app_all_replace_join <- panel_app_all_replace_join %>%
  mutate(color = if_else(color == "green", "Green", color))

fwrite(panel_app_all_replace_join, sep = "\t",
       file = paste0("path/to/PanelAppENG-AUS_panels_PUB-EDIT_", 
                     date_PA, 
                     ".txt"))

# remove 'Removed' genes and duplicates
panel_app_all_dupe <- panel_app_all_replace_join %>%
  filter(color != "Removed") %>%
  group_by(panelapp, panel, id, entity_name) %>%
  summarize(across(everything(), paste_unique), .groups = "drop") 
  
# join NCBI/entrez gene ID
hugo <- fread("path/to/hgnc_total_approved_symbols.txt", sep = "\t")

panel_app_all_dupe_hugo <- panel_app_all_dupe %>%
  left_join(hugo %>%
              select(hgnc_id, entrez_id), by = "hgnc_id")

fwrite(panel_app_all_dupe_hugo, sep = "\t",
       file = paste0("path/to/PanelAppENG-AUS_panels_PUB-EDIT_NO-DUPLICATES_", 
                     date_PA, 
                     ".txt"))

## Genes ####
### read in data ####
#### PanelApp ####
panel_app_all <- fread("path/to/PanelAppENG-AUS_panels_PUB-EDIT_NO-DUPLICATES_04-11-2025.txt", sep = "\t")
panel_app_eng <- panel_app_all[panel_app_all$panelapp == "ENG", ]
panel_app_aus <- panel_app_all[panel_app_all$panelapp == "AUS", ]
panel_app_combo <- panel_app_all %>%
  mutate(disease = case_when(panelapp == "AUS" & id == 69 ~ "CDH", 
                             (panelapp == "ENG" & id == 92 | panelapp == "ENG" & id == 650) ~ "CHH",
                             panelapp == "AUS" & id == 98 ~ "DBA",
                             (panelapp == "ENG" & id == 544 | panelapp == "AUS" & id == 78) ~ "CHO",
                             (panelapp == "ENG" & id == 283 | panelapp == "AUS" & id == 263) ~ "CKD",
                             (panelapp == "ENG" & id == 9 | panelapp == "AUS" & id == 99) ~ "DSD",
                             (panelapp == "ENG" & id == 478 | panelapp == "AUS" & id == 3763) ~ "FET",
                             (panelapp == "ENG" & id == 312 | panelapp == "AUS" & id == 3894) ~ "HOT",
                             panelapp == "ENG" & id == 480 ~ "HPT",
                             (panelapp == "ENG" & id == 285 | panelapp == "AUS" & id == 250) ~ "INT",
                             (panelapp == "ENG" & id == 112 | panelapp == "AUS" & id == 203) ~ "MIT",
                             (panelapp == "ENG" & id == 85 | panelapp == "AUS" & id == 3120) ~ "NEU",
                             panelapp == "ENG" & id == 482 ~ "RKT",
                             (panelapp == "ENG" & id == 548 | panelapp == "AUS" & id == 199) ~ "TKD")) %>%
  mutate(id = disease) %>%
  filter(!is.na(disease)) %>%
  distinct(disease, entity_name, .keep_all = TRUE)

#### Disease context ####
pmids_1_all <- fread("path/to/All_1_PMIDs.txt", sep = "\t", header = FALSE)
length(pmids_1_all$V1)
length(unique(pmids_1_all$V1))

pmids_1 <- list.files(path = "path/to/",
                      pattern = "[A-Z]+_1_PMIDs.txt", full.names = TRUE)
names(pmids_1) <- substr(list.files(path = "path/to/",
                                    pattern = "[A-Z]+_1_PMIDs.txt"), 1, 3)

pmids_1_df <- rbindlist(lapply(pmids_1, fread, col.names = "PMID"), idcol = "disease")
length(pmids_1_df$PMID)
length(unique(pmids_1_all$V1))

df_disease <- pmids_1_df %>%
  group_by(disease) %>%
  group_split(.keep = FALSE)
names(df_disease) <- pmids_1_df %>%
  group_by(disease) %>% 
  group_keys() %>%
  pull(disease)
lapply(df_disease, nrow)

#### PubTator 3.0 - FTP ####
pt3ftp_gene <- fread("path/to/PT3_Genes_1_PMIDs.txt", sep = "\t")
pt3ftp_var <- fread("path/to/PT3_Variants_1_PMIDs.txt", sep = "\t")
sum(pt3ftp_gene$PMID %in% pmids_1_all$V1)/nrow(pt3ftp_gene)
sum(pt3ftp_var$PMID %in% pmids_1_all$V1)/nrow(pt3ftp_var)
length(unique(pt3ftp_gene$PMID))
length(unique(pmids_1_all$V1))
length(unique(pmids_1_all$V1)) - length(unique(pt3ftp_gene$PMID))
(length(unique(pmids_1_all$V1)) - length(unique(pt3ftp_gene$PMID)))/length(unique(pmids_1_all$V1)) # closer to 1 = more missing

pt3ftp_df_all <- pt3ftp_gene %>%
  select(PMID, ID) %>%
  separate_longer_delim(ID, delim = ";") %>%
  rename(Genes = ID) %>%
  bind_rows(pt3ftp_var %>%
              mutate(Variants = paste0(ID, "$", Name)) %>%
              select(PMID, Variants)
  ) %>%
  group_by(PMID) %>%
  reframe(Genes = paste_unique(Genes, "~"), Variants = paste_unique(Variants, sep = "~"))

#### PubTator 3.0 - API ####
annotatoR(ids = unique(pmids_1_df$PMID))
pt3api_df_all <- fread("path/to/PT_Master.txt", sep = "\t") %>%
  filter(PMID %in% as.character(unique(pmids_1_df$PMID)))
pt3api_df_all$PMID <- as.numeric(pt3api_df_all$PMID)
length(unique(pt3api_df_all$PMID))

pt3api_df_all <- pt3api_df_all %>%
  separate_longer_delim(Genes, delim = "|") %>%
  mutate(Genes = ifelse((Genes == "" | grepl("NA~", Genes)), NA, Genes),
         Variants = ifelse((Variants == "" | grepl("NA~", Variants)), NA, Variants)) %>%
  separate_wider_delim(Genes, delim = "~", names = c("Genes", "Gene_names")) %>%
  select(PMID, Genes, Variants) %>%
  group_by(PMID) %>%
  reframe(Genes = paste_unique(Genes, "~"), Variants = paste_unique(Variants, sep = "~"))
dim(pt3api_df_all)

ftp_not_in_api <- unique(pt3ftp_df_all$PMID)[!unique(pt3ftp_df_all$PMID) %in% unique(pt3api_df_all$PMID)]
api_not_in_ftp <- unique(pt3api_df_all$PMID)[!unique(pt3api_df_all$PMID) %in% unique(pt3ftp_df_all$PMID)]
length(ftp_not_in_api) # all ftp PMIDS in api
length(api_not_in_ftp)# api PMIDs not in ftp = 189859
length(unique(pt3ftp_df_all$PMID)) # number of unique ftp PMIDs
length(unique(pt3api_df_all$PMID)) - length(api_not_in_ftp) # number api PMIDs - api PMIDs not in ftp = number unique ftp PMIDs

api_gene <- pt3api_df_all %>% # api Genes
  filter(!is.na(Genes)) %>%
  pull(PMID)

ftp_gene <- pt3ftp_df_all %>% # ftp Genes
  filter(!is.na(Genes)) %>%
  pull(PMID)

api_gene_not_in_ftp <- api_gene[!api_gene %in% ftp_gene] # PMIDs of api Genes not in ftp Genes
fpp_gene_not_in_api <- ftp_gene[!ftp_gene %in% api_gene] # PMIDs of ftp Genes not in api Genes

pt3_df_Genes <- pt3api_df_all %>%     # PT3 FTP and API datasets combined 
  mutate(source = "API") %>%
  bind_rows(pt3ftp_df_all %>%
              mutate(source = "FTP")) %>%
  separate_longer_delim(Genes, delim = "~") %>%
  mutate(Genes = as.integer(Genes)) %>%
  group_by(PMID, source) %>%
  reframe(Genes = paste_unique(sort(Genes), sep = "~"))

pt3_genes_not_equal <- pt3_df_Genes %>%  # where FTP and API Genes are not equal 
  pivot_wider(names_from = source, values_from = Genes) %>%
  filter(FTP != API | (is.na(FTP) & !is.na(API)) | (!is.na(FTP) & is.na(API)))
nrow(pt3_genes_not_equal)

pt3_genes_not_equal_na <- pt3_df_Genes %>%  # where FTP is NA and API is not, or vice versa 
  pivot_wider(names_from = source, values_from = Genes) %>%
  filter((is.na(FTP) & !is.na(API)) | (!is.na(FTP) & is.na(API)))
nrow(pt3_genes_not_equal_na)
nrow(subset(pt3_genes_not_equal_na, is.na(pt3_genes_not_equal_na$API) & # where API is NA, not FTP
              !is.na(pt3_genes_not_equal_na$FTP)))
nrow(subset(pt3_genes_not_equal_na, is.na(pt3_genes_not_equal_na$FTP) & # where FTP is NA, not API
              !is.na(pt3_genes_not_equal_na$API)))

#### PubTator Central ####
ptc_gene <- fread("path/to/PTC_Genes_1_PMIDs.txt", sep = "\t")
ptc_var <- fread("path/to/PTC_Variants_1_PMIDs.txt", sep = "\t")
sum(ptc_gene$PMID %in% pmids_1_all$V1)/nrow(ptc_gene)
sum(ptc_var$PMID %in% pmids_1_all$V1)/nrow(ptc_var)
length(unique(ptc_gene$PMID))
length(unique(pmids_1_all$V1))
length(unique(pmids_1_all$V1)) - length(unique(ptc_gene$PMID))
(length(unique(pmids_1_all$V1)) - length(unique(ptc_gene$PMID))) / length(unique(pmids_1_all$V1))

ptc_df_all <- ptc_gene %>%
  select(PMID, ID) %>%
  
  separate_longer_delim(ID, delim = ";") %>%
  rename(Genes = ID) %>%
  bind_rows(ptc_var %>%
              mutate(Variants = paste0(ID, "$", Name)) %>%
              select(PMID, Variants)
  ) %>%
  group_by(PMID) %>%
  reframe(Genes = paste_unique(Genes, "~"), Variants = paste_unique(Variants, sep = "~"))  

#### PubTator 3.0 Relations ####
## MeSH Descriptors for each Disease Context
mesh_pub <- fread("path/to/MeSH_PubTator_List_10-11-2025.txt", sep = "\t")
mesh_pub_ready <- unique(str_trim(mesh_pub$MeSH))

## Data containing MeSH Descriptors and their MeSH unique IDs
mesh_ids <- fread("path/to/mesh_terms_2025.txt", sep = "\t", header = FALSE)

mesh_names_ids <- mesh_ids %>% 
  filter(V1 %in% mesh_pub_ready)

mesh_names_ids %>%
  pull(V2) %>%
  list() %>%
  fwrite("path/to/mesh_term_ids_Q1_10-11-2025.txt")

# run filter_relations.py -r relation2pubtator3_17-08-2025.gz -m mesh_term_ids_Q1_10-11-2025.txt -o PT3_Relations_1_PMIDs.txt

relations_query <- fread("path/to/PT3_Relations_1_PMIDs.txt", sep = "\t")
rel_gene <- relations_query %>%
  filter((Type1 == "Disease" & Type2 == "Gene") | (Type2 == "Disease" & Type1 == "Gene")) %>%
  mutate(MeSH = ifelse(Type1 == "Disease", ID1, ID2),
         Genes = ifelse(Type1 == "Gene", ID1, ID2)) %>%
  left_join(mesh_ids, by = join_by(MeSH == V2)) %>%
  left_join(mesh_pub %>%
              distinct(), by = join_by(V1 == MeSH)) %>%
  rename(MeSH_term = V1) %>%
  select(PMID, Type, MeSH, Genes, MeSH_term, disease)
rel_var <- relations_query %>%
  filter(PMID %in% rel_gene$PMID) %>%
  filter(Type1 == "DNAMutation" | Type2 == "DNAMutation") %>%
  mutate(Variants = ifelse(Type1 == "DNAMutation",  ID1, ID2)) %>%
  select(PMID, Variants)

rel_df_all <- rel_gene %>%
  bind_rows(rel_var) %>%
  group_by(PMID) %>%
  reframe(Genes = paste_unique(Genes, "~"), Variants = paste_unique(Variants, sep = "~"),
          disease = paste_unique(disease))  



df_all <- pt3ftp_df_all %>% 
  mutate(source = "FTP") %>% 
  bind_rows(pt3api_df_all %>%
              mutate(source = "API")) %>% 
  bind_rows(ptc_df_all %>% 
              mutate(source = "PTC")) %>% 
  bind_rows(rel_df_all %>% 
              mutate(source = "REL"))

fwrite(df_all, sep = "\t", "path/to/PT_All_1_PMIDs.txt")

### variants ####
df_all_var <- pt3ftp_df_all %>% 
  mutate(source = "FTP") %>% 
  bind_rows(pt3api_df_all %>%
              mutate(source = "API")) %>% 
  bind_rows(ptc_df_all %>% 
              mutate(source = "PTC")) %>% 
  bind_rows(rel_df_all %>% 
              mutate(source = "REL")) %>%
  filter(!is.na(Variants)) %>%
  mutate(Variants = ifelse(source == "API", gsub("\\|", "~", gsub("~", "$", Variants)), Variants))

rs_nums <- df_all_var %>%
  separate_longer_delim(Variants, delim = "~") %>%
  separate_wider_delim(Variants, delim = "$", names = c("varID", "varName"), too_few = "align_start") %>%
  mutate(varID = ifelse(grepl("^RS#:", varID) & !grepl("CorrespondingGene:", varID), 
                        paste0("rs", str_extract(varID, "(?<=RS#:)[0-9]+")), 
                        varID)) %>%
  filter(grepl("^rs[0-9]", varID)) %>%
  mutate(varID = gsub("##", "", varID)) %>%
  select(varID) %>%
  distinct()
dim(rs_nums)

fwrite(rs_nums, "path/to/LitVar_1_RS.txt", sep = "\t")

### gene counts ####
arg_list <- list(
  # CDH
  list(pa = panel_app_aus, pa_id = 69),
  # CHH
  list(pa = panel_app_eng, pa_id = 92), list(pa = panel_app_eng, pa_id = 650), list(pa = panel_app_combo, pa_id = "CHH"),
  # CHO
  list(pa = panel_app_eng, pa_id = 544), list(pa = panel_app_aus, pa_id = 78), list(pa = panel_app_combo, pa_id = "CHO"),
  # CKD
  list(pa = panel_app_eng, pa_id = 283), list(pa = panel_app_aus, pa_id = 263), list(pa = panel_app_combo, pa_id = "CKD"),
  # DBA
  list(pa = panel_app_aus, pa_id = 98),
  # DSD
  list(pa = panel_app_eng, pa_id = 9), list(pa = panel_app_aus, pa_id = 99), list(pa = panel_app_combo, pa_id = "DSD"),
  # FET
  list(pa = panel_app_eng, pa_id = 478), list(pa = panel_app_aus, pa_id = 3763), list(pa = panel_app_combo, pa_id = "FET"),
  # HOT
  list(pa = panel_app_eng, pa_id = 312), list(pa = panel_app_aus, pa_id = 3894), list(pa = panel_app_combo, pa_id = "HOT"),
  # HPT 
  list(pa = panel_app_eng, pa_id = 480),
  # INT
  list(pa = panel_app_eng, pa_id = 285), list(pa = panel_app_aus, pa_id = 250), list(pa = panel_app_combo, pa_id = "INT"),
  # MIT
  list(pa = panel_app_eng, pa_id = 112), list(pa = panel_app_aus, pa_id = 203), list(pa = panel_app_combo, pa_id = "MIT"),
  # NEU
  list(pa = panel_app_eng, pa_id = 85), list(pa = panel_app_aus, pa_id = 3120), list(pa = panel_app_combo, pa_id = "NEU"),
  # RKT
  list(pa = panel_app_eng, pa_id = 482),
  # TKD
  list(pa = panel_app_eng, pa_id = 548), list(pa = panel_app_aus, pa_id = 199), list(pa = panel_app_combo, pa_id = "TKD"))

#### PubTator 3.0 - FTP ####
pt3ftp_df_disease <- lapply(df_disease, function(a) left_join(x = a, y = pt3ftp_df_all))
pt3ftp_df_disease <- unlist(lapply(pt3ftp_df_disease, function(i) list(i, i, i)), recursive = FALSE)

names(pt3ftp_df_disease) <- case_when(grepl("1$", names(pt3ftp_df_disease)) ~ gsub("1$", "_ENG", names(pt3ftp_df_disease)), 
                                      grepl("2$", names(pt3ftp_df_disease)) ~ gsub("2$", "_AUS", names(pt3ftp_df_disease)),
                                      grepl("3$", names(pt3ftp_df_disease)) ~ gsub("3$", "_COMBO", names(pt3ftp_df_disease)))

pt3ftp_df_disease$CDH_ENG <- pt3ftp_df_disease$CDH_COMBO <- pt3ftp_df_disease$DBA_ENG <- pt3ftp_df_disease$DBA_COMBO <- pt3ftp_df_disease$HPT_AUS <- pt3ftp_df_disease$HPT_COMBO <- pt3ftp_df_disease$RKT_AUS <- pt3ftp_df_disease$RKT_COMBO <- NULL

unlist(lapply(pt3ftp_df_disease, nrow))
unlist(lapply(pt3ftp_df_disease, function(x) nrow(x[!is.na(x$Genes), ])))
length(unique(pt3ftp_df_all$PMID)) # unique article count

pt3ftp_gene_counts <- Map(function(df, args) {do.call(gene_countR, c(list(df), args))}, pt3ftp_df_disease, arg_list)

#### PubTator 3.0 - API ####
pt3api_df_disease <- lapply(df_disease, function(a) left_join(x = a, y = pt3api_df_all))
pt3api_df_disease <- unlist(lapply(pt3api_df_disease, function(i) list(i, i, i)), recursive = FALSE)

names(pt3api_df_disease) <- case_when(grepl("1$", names(pt3api_df_disease)) ~ gsub("1$", "_ENG", names(pt3api_df_disease)), 
                                      grepl("2$", names(pt3api_df_disease)) ~ gsub("2$", "_AUS", names(pt3api_df_disease)),
                                      grepl("3$", names(pt3api_df_disease)) ~ gsub("3$", "_COMBO", names(pt3api_df_disease)))

pt3api_df_disease$CDH_ENG <- pt3api_df_disease$CDH_COMBO <- pt3api_df_disease$DBA_ENG <- pt3api_df_disease$DBA_COMBO <- pt3api_df_disease$HPT_AUS <- pt3api_df_disease$HPT_COMBO <- pt3api_df_disease$RKT_AUS <- pt3api_df_disease$RKT_COMBO <- NULL

unlist(lapply(pt3api_df_disease, nrow))
unlist(lapply(pt3api_df_disease, function(x) nrow(x[!is.na(x$Genes), ])))

pt3api_gene_counts <- Map(function(df, args) {do.call(gene_countR, c(list(df), args))}, pt3api_df_disease, arg_list)

#### PubTator Central ####
ptc_df_disease <- lapply(df_disease, function(a) left_join(x = a, y = ptc_df_all))
ptc_df_disease <- unlist(lapply(ptc_df_disease, function(i) list(i, i, i)), recursive = FALSE)

names(ptc_df_disease) <- case_when(grepl("1$", names(ptc_df_disease)) ~ gsub("1$", "_ENG", names(ptc_df_disease)), 
                                   grepl("2$", names(ptc_df_disease)) ~ gsub("2$", "_AUS", names(ptc_df_disease)),
                                   grepl("3$", names(ptc_df_disease)) ~ gsub("3$", "_COMBO", names(ptc_df_disease)))

ptc_df_disease$CDH_ENG <- ptc_df_disease$CDH_COMBO <- ptc_df_disease$DBA_ENG <- ptc_df_disease$DBA_COMBO <- ptc_df_disease$HPT_AUS <- ptc_df_disease$HPT_COMBO <- ptc_df_disease$RKT_AUS <- ptc_df_disease$RKT_COMBO <- NULL

ptc_gene_counts <- Map(function(df, args) {do.call(gene_countR, c(list(df), args))}, ptc_df_disease, arg_list)

#### PubTator 3.0 Relations ####
rel_df_disease <- rel_df_all %>% 
  separate_longer_delim(disease, delim = ";") %>%
  group_by(disease) %>%
  group_split(.keep = FALSE)
names(rel_df_disease) <- rel_df_all %>%
  separate_longer_delim(disease, delim = ";") %>%
  group_by(disease) %>% 
  group_keys() %>%
  pull(disease)
lapply(rel_df_disease, nrow)
lapply(rel_df_disease, function(x) nrow(x[!is.na(x$Genes), ])) # all PMIDS contain Gene annotation as expected

rel_df_disease <- lapply(rel_df_disease, function(a) left_join(x = a, y = rel_df_all))
rel_df_disease <- unlist(lapply(rel_df_disease, function(i) list(i, i, i)), recursive = FALSE)

names(rel_df_disease) <- case_when(grepl("1$", names(rel_df_disease)) ~ gsub("1$", "_ENG", names(rel_df_disease)), 
                                   grepl("2$", names(rel_df_disease)) ~ gsub("2$", "_AUS", names(rel_df_disease)),
                                   grepl("3$", names(rel_df_disease)) ~ gsub("3$", "_COMBO", names(rel_df_disease)))

rel_df_disease$CDH_ENG <- rel_df_disease$CDH_COMBO <- rel_df_disease$DBA_ENG <- rel_df_disease$DBA_COMBO <- rel_df_disease$HPT_AUS <- rel_df_disease$HPT_COMBO <- rel_df_disease$RKT_AUS <- rel_df_disease$RKT_COMBO <- NULL

rel_gene_counts <- Map(function(df, args) {do.call(gene_countR, c(list(df), args))}, rel_df_disease, arg_list)

### stats ####
#### PubTator 3.0 - FTP ####
pt3ftp_stats_list <- rbindlist(lapply(pt3ftp_gene_counts, function(x) x$stats), idcol = "comp")
pt3ftp_stats_list <- pt3ftp_stats_list %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "PT3FTP") %>%
  select(-comp)

#### PubTator 3.0 - API ####
pt3api_stats_list <- rbindlist(lapply(pt3api_gene_counts, function(x) x$stats), idcol = "comp")
pt3api_stats_list <- pt3api_stats_list %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "PT3API") %>%
  select(-comp)

#### PubTator Central ####
ptc_stats_list <- rbindlist(lapply(ptc_gene_counts, function(x) x$stats), idcol = "comp")
ptc_stats_list <- ptc_stats_list %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "PTC") %>%
  select(-comp)

#### PubTator 3.0 Relations ####
rel_stats_list <- rbindlist(lapply(rel_gene_counts, function(x) x$stats), idcol = "comp")
rel_stats_list <- rel_stats_list %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "REL") %>%
  select(-comp)

all_stats_list <- pt3ftp_stats_list %>% 
  bind_rows(pt3api_stats_list) %>%
  bind_rows(ptc_stats_list) %>% 
  bind_rows(rel_stats_list)

fwrite(all_stats_list, sep = "\t", "path/to/stats_df_20-11-2025.txt")

all_stats_list_table <- all_stats_list %>%
  filter(grepl("tpr", stat)) %>%
  mutate(value2 = ifelse(stat == "tpr", value, NA)) %>% 
  group_by(disease, panelapp) %>%
  fill(value2, .direction = "downup") %>%
  filter(row_number() != 1) %>% 
  ungroup() %>%
  pivot_wider(names_from = source, values_from = c(value, value2)) %>% 
  select(disease, panelapp, value2_PT3FTP, value_PT3FTP, value2_PT3API, value_PT3API,
         value2_PTC, value_PTC, value2_REL, value_REL) %>%
  filter(!is.na(value2_PT3FTP)) %>% 
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         panelapp = factor(panelapp, levels = c("COMBO", "ENG", "AUS"))) %>%
  filter(panelapp == "COMBO" | disease %in% c("CDH", "DBA", "HPT", "RKT"))

fwrite(all_stats_list_table, sep = "\t", 
       "path/to/tpr_table_20-11-2025.txt")

### genes ####
#### PubTator 3.0 - FTP ####
pt3ftp_genes_list <- rbindlist(lapply(pt3ftp_gene_counts, function(x) x$gene_count), idcol = "comp")
pt3ftp_genes_list <- pt3ftp_genes_list %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp_instance"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "PT3FTP") %>%
  select(-comp)

#### PubTator 3.0 - API ####
pt3api_genes_list <- rbindlist(lapply(pt3api_gene_counts, function(x) x$gene_count), idcol = "comp")
pt3api_genes_list <- pt3api_genes_list %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp_instance"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "PT3API") %>%
  select(-comp)

#### PubTator Central ####
ptc_genes_list <- rbindlist(lapply(ptc_gene_counts, function(x) x$gene_count), idcol = "comp")
ptc_genes_list <- ptc_genes_list %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp_instance"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "PTC") %>%
  select(-comp)

#### PubTator 3.0 Relations ####
rel_genes_list <- rbindlist(lapply(rel_gene_counts, function(x) x$gene_count), idcol = "comp")
rel_genes_list <- rel_genes_list %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp_instance"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "REL") %>%
  select(-comp)

#### HGNC ####
hugo <- fread("path/to/hgnc_complete_set.txt", sep = "\t")
hugo_filt <- hugo %>%
  select(entrez_id, hgnc_id, locus_group)

all_genes_list <- pt3ftp_genes_list %>% 
  bind_rows(pt3api_genes_list) %>% 
  bind_rows(ptc_genes_list) %>% 
  bind_rows(rel_genes_list) %>%
  left_join(hugo_filt, by = join_by(Genes == entrez_id))

fwrite(all_genes_list, sep = "\t", "path/to/genes_df_20-11-2025.txt")

####missing ####
pt3ftp_missing_list <- rbindlist(lapply(pt3ftp_gene_counts, function(x) list(missing = x$missing)),
                                 idcol = "comp") %>%
  mutate(source = "PT3FTP")
pt3api_missing_list <- rbindlist(lapply(pt3api_gene_counts, function(x) list(missing = x$missing)),
                                 idcol = "comp") %>%
  mutate(source = "PT3API")
ptc_missing_list <- rbindlist(lapply(ptc_gene_counts, function(x) list(missing = x$missing)),
                              idcol = "comp") %>%
  mutate(source = "PTC")
rel_missing_list <- rbindlist(lapply(rel_gene_counts, function(x) list(missing = x$missing)),
                              idcol = "comp") %>%
  mutate(source = "REL")

all_missing_list <- pt3ftp_missing_list %>%
  bind_rows(pt3api_missing_list) %>%
  bind_rows(ptc_missing_list) %>%
  bind_rows(rel_missing_list) %>%
  separate_longer_delim(cols = "missing", delim = ";") %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp")) %>%
  mutate(entity_id = paste0(disease, "_", missing)) %>%
  group_by(entity_id) %>%
  reframe(disease = paste_unique(disease), missing = paste_unique(missing), source = paste_unique(source))

fwrite(all_missing_list, sep = "\t", "path/to/missing_df_20-11-2025.txt")

### missing genes ####
## PanelApp PMIDs to file
panelapp_eng_id = c(544, 92, 650, 283, 9, 480, 312, 478, 85, 285, 112, 482, 548)
panelapp_aus_id = c(78, 69, 263, 99, 3894, 3763, 3120, 250, 203, 199)

pa_genes_all <- panel_app_all %>% 
  filter((panelapp == "ENG" & id %in% panelapp_eng_id) | (panelapp == "AUS" & id %in% panelapp_aus_id)) %>% 
  select(panelapp, panel, id, entity_name, phenotypes, publications2, color, entrez_id)

pa_genes_pmids <- pa_genes_all %>%
  separate_longer_delim(cols = publications2, delim = ";") %>%
  filter(!is.na(publications2)) %>%
  pull(publications2) %>%
  unique()

# write PMIDs of PanelApp publications from missing genes for MeSH annotating 
fwrite(list(pa_genes_pmids), "path/to/PanelApp_Publications_PMIDs_05-11-2025.txt")

## join MeSH
lines <- readLines("path/to/PanelApp_Publications_PMIDs_MeSH_05-11-2025.txt")

pmid <- NULL
mh_list <- list()
current_pmid <- NULL

for (line in lines) {
  if (grepl("^PMID- ", line)) {
    current_pmid <- gsub("PMID- ", "", line)
    pmid <- c(pmid, current_pmid)
    mh_list[[current_pmid]] <- character(0)
  } else if (grepl("^MH  - ", line) && !is.null(current_pmid)) {
    term <- gsub("MH  - ", "", line)
    mh_list[[current_pmid]] <- c(mh_list[[current_pmid]], term)
  }
}

panelapp_sources_mesh <- data.frame(
  PMID = names(mh_list),
  MH = sapply(mh_list, function(x) paste(x, collapse = "~")),
  row.names = NULL,
  stringsAsFactors = FALSE
)

sum(panelapp_sources_mesh$PMID %in% pa_genes_pmids)/length(panelapp_sources_mesh$PMID)
sum(pa_genes_pmids %in% panelapp_sources_mesh$PMID)/length(pa_genes_pmids)
pa_genes_pmids[!pa_genes_pmids %in% panelapp_sources_mesh$PMID] # these PMIDs are not retrieved on PubMed

pa_genes_all_mesh <- pa_genes_all %>%
  separate_longer_delim(cols = publications2, delim = ";") %>%
  rename(publications = publications2) %>%
  left_join(panelapp_sources_mesh, by = join_by(publications == PMID)) %>%
  mutate(disease = case_when(id %in% c(544, 78) ~ "CHO", id %in% c(69) ~ "CDH", id %in% c(92, 650) ~ "CHH", 
                             id %in% c(283, 263) ~ "CKD", id %in% c(9, 99) ~ "DSD", id %in% c(480) ~ "HPT", 
                             id %in% c(312, 3894) ~ "HOT", id %in% c(478, 3763) ~ "FET", id %in% c(3120, 85) ~ "NEU", 
                             id %in% c(285, 250) ~ "INT", id %in% c(203, 112) ~ "MIT", id %in% c(482) ~ "RKT", 
                             id %in% c(199, 548) ~ "TKD")) %>%
  mutate(entity_id = paste0(disease, "_", entity_name))

fwrite(pa_genes_all_mesh, sep = "\t",
       "path/to/PanelApp_Publications_PMID_MeSH_05-11-2025.txt")

## missing gene MeSH
pa_genes_all_mesh_missing <- pa_genes_all_mesh %>% 
  left_join(all_missing_list, by = join_by(entity_id == entity_id, entity_name == missing, 
                                           disease == disease)) %>%
  select(-entity_id) %>%
  filter(!is.na(source))

pa_genes_mesh_phenos <- pa_genes_all_mesh_missing %>%
  select(source, disease, entity_name, phenotypes, MH, publications, color, entrez_id) %>%
  separate_longer_delim(cols = "phenotypes", delim = "~") %>%
  separate_longer_delim(cols = "MH", delim = "~") %>%
  mutate(MH = gsub("\\*", "", MH),
         phenotypes = str_extract(phenotypes, "^[^,]+"),
         MH = str_extract(MH, "^[^/]+")) %>%
  group_by(disease, entity_name) %>%
  reframe(missing_from = paste_unique(source), PMIDs = paste_unique(publications), phenos = paste_unique(phenotypes), 
          MHs = paste_unique(MH), color = paste_unique(color), entrez_id = paste_unique(entrez_id)) %>%
  arrange(disease, entity_name)

pa_genes_mesh_phenos[pa_genes_mesh_phenos == "NA"] <- ""
pa_genes_mesh_phenos$phenos <- gsub("NA", "", pa_genes_mesh_phenos$phenos)
pa_genes_mesh_phenos$phenos <- gsub("^;|;$", "", pa_genes_mesh_phenos$phenos)
pa_genes_mesh_phenos$PMIDs <- gsub("NA", "", pa_genes_mesh_phenos$PMIDs)
pa_genes_mesh_phenos$PMIDs <- gsub("^;|;$", "", pa_genes_mesh_phenos$PMIDs)
pa_genes_mesh_phenos$MHs <- gsub("NA", "", pa_genes_mesh_phenos$MHs)
pa_genes_mesh_phenos$MHs <- gsub("^;|;$", "", pa_genes_mesh_phenos$MHs)

gene_list <- fread("path/to/genes_df_20-11-2025.txt", sep = "\t")
pa_genes_mesh_phenos <- pa_genes_mesh_phenos %>%
  separate_longer_delim(PMIDs, delim = ";") %>%
  mutate(PMIDs = as.integer(PMIDs), entrez_id = as.integer(entrez_id)) %>%
  left_join(pmids_1_df %>%
              mutate(present = PMID), 
            by = join_by(PMIDs == PMID, disease == disease)) %>%
  left_join(gene_list %>%
              select(disease, Genes, PMID_human, PMID_ortho, PMID_var),
            by = join_by(disease == disease, entrez_id == Genes)) %>%
  rowwise() %>%
  mutate(query_pmid = paste_unique(c(strsplit(PMID_human, ";")[[1]], 
                                     strsplit(PMID_ortho, ";")[[1]],
                                     strsplit(PMID_var, ";")[[1]]))) %>%
  select(-PMID_human, -PMID_ortho, -PMID_var) %>% 
  group_by(disease, entity_name) %>%
  reframe(missing_from = missing_from, phenos = paste_unique(phenos), MHs = paste_unique(MHs), 
          color = paste_unique(color), PMID = paste_unique(PMIDs), present = paste_unique(present),
          query_pmid = paste_unique(query_pmid)) %>% 
  distinct(disease, entity_name, .keep_all = TRUE) %>%
  mutate(explanation = case_when(is.na(PMID) & is.na(phenos) ~ "Unknown",
                                 is.na(PMID) & !is.na(phenos) ~ "Phenotype not in query",
                                 !grepl("PT3FTP|PTC", missing_from) & !is.na(present) ~ "Text",
                                 !grepl("PT3FTP|PTC", missing_from) & is.na(present) ~ "Text;Phenotype not in query",
                                 (missing_from == "PT3FTP;PT3API;PTC;REL" & (!is.na(PMID) | !is.na(phenos))) |
                                   (!grepl("PT3FTP", missing_from) & grepl("PTC", missing_from)) | 
                                   (grepl("PT3FTP", missing_from) & !grepl("PTC", missing_from)) |
                                   (grepl("PT3FTP;PTC", missing_from)) ~ "Investigate")) %>%
  select(disease, entity_name, color, missing_from, phenos, MHs, PMID, present, explanation, query_pmid)

query_genes_pmids <- pa_genes_mesh_phenos %>%
  separate_longer_delim(cols = query_pmid, delim = ";") %>%
  filter(!is.na(query_pmid)) %>%
  select(query_pmid) %>%
  distinct()

# write PMIDs of PanelApp publications from missing genes for MeSH annotating 
fwrite(query_genes_pmids, "path/to/MissingFromQuery_PMIDs_20-11-2025.txt")

## join MeSH
lines2 <- readLines("path/to/MissingFromQuery_PMIDs_20-11-2025_MeSH.txt")

pmid2 <- NULL
mh_list2 <- list()
current_pmid2 <- NULL

for (line in lines2) {
  if (grepl("^PMID- ", line)) {
    current_pmid2 <- gsub("PMID- ", "", line)
    pmid2 <- c(pmid2, current_pmid2)
    mh_list2[[current_pmid2]] <- character(0)
  } else if (grepl("^MH  - ", line) && !is.null(current_pmid2)) {
    term <- gsub("MH  - ", "", line)
    mh_list2[[current_pmid2]] <- c(mh_list2[[current_pmid2]], term)
  }
}

queries_sources_mesh <- data.frame(
  PMID = names(mh_list2),
  MH = sapply(mh_list2, function(x) paste(x, collapse = "~")),
  row.names = NULL,
  stringsAsFactors = FALSE
)
queries_sources_mesh <- queries_sources_mesh %>% 
  mutate(MH = gsub("\\*", "", MH),
         MH = str_extract(MH, "^[^/]+"))

mesh_pheno_recurrence <- pa_genes_mesh_phenos %>%
  select(disease, MHs) %>%
  separate_longer_delim(MHs, delim = ";") %>%
  mutate(term = tolower(gsub("[0-9]", "", MHs))) %>%
  count(disease, term) %>%
  mutate(type = "mesh") %>%
  bind_rows(pa_genes_mesh_phenos %>%
              separate_longer_delim(phenos, delim = ";") %>%
              mutate(term = tolower(gsub("\\(.+\\)", "", gsub("[0-9]", "", phenos)))) %>%
              count(disease, term) %>%
              mutate(type = "panelapp_pheno")) %>%
  arrange(disease, desc(n))

fwrite(mesh_pheno_recurrence, sep = "\t", "path/to/mesh_pheno_count.txt")

### relations update ####
## MeSH Descriptors for each Disease Context
mesh_pub2 <- fread("path/to/MeSH_PubTator_List_v2_20-11-2025.txt", sep = "\t")
mesh_pub2_ready <- unique(str_trim(mesh_pub2$MeSH))

## Data containing MeSH Descriptors and their MeSH unique IDs
mesh_ids <- fread("path/to/mesh_terms_2025.txt", sep = "\t", header = FALSE)
supp_ids <- fread("path/to/supplementary_2025.txt", sep = "\t", header = FALSE)
mesh_names_ids2 <- mesh_ids %>% 
  filter(V1 %in% mesh_pub2_ready)
supp_names_ids2 <- supp_ids %>%
  filter(V1 %in% mesh_pub2_ready)
mesh_names_ids2 %>%
  bind_rows(supp_names_ids2) %>%
  pull(V2) %>%
  list() %>%
  fwrite("path/to/mesh_term_ids_Q2_20-11-2025.txt")

# run filter_relations.py -r relation2pubtator3_17-08-2025.gz -m mesh_term_ids_Q1_10-11-2025.txt -o PT3_Relations_1_PMIDs.txt

relations_query2 <- fread("path/to/PT3_Relations_2_PMIDs.txt", sep = "\t")
rel_gene2 <- relations_query2 %>%
  filter((Type1 == "Disease" & Type2 == "Gene") | (Type2 == "Disease" & Type1 == "Gene")) %>%
  mutate(MeSH = ifelse(Type1 == "Disease", ID1, ID2),
         Genes = ifelse(Type1 == "Gene", ID1, ID2)) %>%
  left_join(mesh_ids %>% 
              bind_rows(supp_ids), 
            by = join_by(MeSH == V2)) %>%
  left_join(mesh_pub2 %>%
              distinct(), by = join_by(V1 == MeSH)) %>%
  rename(MeSH_term = V1) %>%
  select(PMID, Type, MeSH, Genes, MeSH_term, disease)
rel_var2 <- relations_query2 %>%
  filter(PMID %in% rel_gene2$PMID) %>%
  filter(Type1 == "DNAMutation" | Type2 == "DNAMutation") %>%
  mutate(Variants = ifelse(Type1 == "DNAMutation", ID1, ID2)) %>%
  select(PMID, Variants)

rel_df_all2 <- rel_gene2 %>%
  bind_rows(rel_var2) %>%
  group_by(PMID) %>%
  reframe(Genes = paste_unique(Genes, "~"), Variants = paste_unique(Variants, sep = "~"),
          disease = paste_unique(disease))  

# PubTator 3.0 Relations
rel_df_disease2 <- rel_df_all2 %>% 
  separate_longer_delim(disease, delim = ";") %>%
  group_by(disease) %>%
  group_split(.keep = FALSE)
names(rel_df_disease2) <- rel_df_all2 %>%
  separate_longer_delim(disease, delim = ";") %>%
  group_by(disease) %>% 
  group_keys() %>%
  pull(disease)
lapply(rel_df_disease2, nrow)
lapply(rel_df_disease2, function(x) nrow(x[!is.na(x$Genes), ])) # all PMIDS contain Gene annotation as expected

rel_df_disease2 <- lapply(rel_df_disease2, function(a) left_join(x = a, y = rel_df_all2))
rel_df_disease2 <- unlist(lapply(rel_df_disease2, function(i) list(i, i, i)), recursive = FALSE)

names(rel_df_disease2) <- case_when(grepl("1$", names(rel_df_disease2)) ~ gsub("1$", "_ENG", names(rel_df_disease2)), 
                                    grepl("2$", names(rel_df_disease2)) ~ gsub("2$", "_AUS", names(rel_df_disease2)),
                                    grepl("3$", names(rel_df_disease2)) ~ gsub("3$", "_COMBO", names(rel_df_disease2)))

rel_df_disease2$CDH_ENG <- rel_df_disease2$CDH_COMBO <- rel_df_disease2$DBA_ENG <- rel_df_disease2$DBA_COMBO <- rel_df_disease2$HPT_AUS <- rel_df_disease2$HPT_COMBO <- rel_df_disease2$RKT_AUS <- rel_df_disease2$RKT_COMBO <- NULL

rel_gene_counts2 <- Map(function(df, args) {do.call(gene_countR, c(list(df), args))}, rel_df_disease2, arg_list)

rel_stats_list2 <- rbindlist(lapply(rel_gene_counts2, function(x) x$stats), idcol = "comp")
rel_stats_list2 <- rel_stats_list2 %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "REL2") %>%
  select(-comp)

all_stats_table_update <- all_stats_list %>%
  bind_rows(rel_stats_list2) %>%
  filter(grepl("tpr", stat)) %>%
  mutate(value2 = ifelse(stat == "tpr", value, NA)) %>% 
  group_by(disease, panelapp) %>%
  fill(value2, .direction = "downup") %>%
  filter(row_number() != 1) %>% 
  ungroup() %>%
  pivot_wider(names_from = source, values_from = c(value, value2)) %>% 
  select(disease, panelapp, value2_PT3FTP, value_PT3FTP, value2_PT3API, value_PT3API,
         value2_PTC, value_PTC, value2_REL, value_REL, value2_REL2, value_REL2) %>%
  filter(!is.na(value2_PT3FTP)) %>% 
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         panelapp = factor(panelapp, levels = c("COMBO", "ENG", "AUS"))) %>%
  mutate(across(.cols = value2_PT3FTP:value_REL2, \(x) round(x, 3))) %>%
  filter(panelapp == "COMBO" | disease %in% c("CDH", "DBA", "HPT", "RKT"))

rel_genes_list2 <- rbindlist(lapply(rel_gene_counts2, function(x) x$gene_count), idcol = "comp")
rel_genes_list2 <- rel_genes_list2 %>%
  separate_wider_delim(cols = comp, delim = "_", names = c("disease", "panelapp_instance"), cols_remove = FALSE) %>%
  mutate(disease = case_when(grepl("CHH_ENG", comp) ~ "CHH", 
                             grepl("CHH_AUS", comp) ~ "CHH_GMS",
                             grepl("CHH_COMBO", comp) ~ "CHH_COMBO",
                             !grepl("CHH", comp) ~ disease),
         source = "REL2") %>%
  select(-comp)

hugo <- fread("path/to/hgnc_complete_set.txt", sep = "\t")
hugo_filt <- hugo %>%
  select(entrez_id, hgnc_id, locus_group)

rel_genes_list2 <- rel_genes_list2 %>% 
  left_join(hugo_filt, by = join_by(Genes == entrez_id))

fwrite(rel_stats_list2, sep = "\t", 
       "path/to/stats_REL2_df_20-11-2025.txt")

fwrite(all_stats_table_update, sep = "\t", 
       "path/to/tpr_table_REL2_20-11-2025.txt")

fwrite(rel_genes_list2, sep = "\t", 
       "path/to/genes_REL2_df_20-11-2025.txt")
