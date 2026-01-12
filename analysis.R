# PanelPub Analysis ####
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

## read in Data ####

# annotations 
df_all <- fread("path/to/PT_All_1_PMIDs.txt", sep = "\t")
pmids_1_all <- fread("path/to/All_1_PMIDs.txt", sep = "\t", header = FALSE)

# extracted data
stats_list <- fread("path/to/stats_df_20-11-2025.txt", sep = "\t")
stats_list_rel2 <- fread("path/to/stats_REL2_df_20-11-2025.txt", sep = "\t")
stats_list <- stats_list %>%
  bind_rows(stats_list_rel2)

gene_list <- fread("path/to/genes_df_20-11-2025.txt", sep = "\t")
gene_list_rel2 <- fread("path/to/genes_REL2_df_20-11-2025.txt", sep = "\t")
gene_list <- gene_list %>%
  bind_rows(gene_list_rel2)

missing_list <- fread("path/to/missing_df_20-11-2025.txt", sep = "\t")

# PanelApp
panel_app_all <- fread("path/to/PanelAppENG-AUS_panels_PUB-EDIT_NO-DUPLICATES_04-11-2025.txt", sep = "\t")
panelapp_eng_id = c(544, 92, 650, 283, 9, 480, 312, 478, 85, 285, 112, 482, 548)
panelapp_aus_id = c(78, 69, 263, 98, 99, 3894, 3763, 3120, 250, 203, 199)
panel_app_all <- panel_app_all %>%
  filter((panelapp == "ENG" & id %in% panelapp_eng_id) | (panelapp == "AUS" & id %in% panelapp_aus_id)) %>%
  rename(Gene = entity_name, panelapp_instance = panelapp) %>%
  mutate(disease = case_when(panelapp_instance == "AUS" & id == 69 ~ "CDH", panelapp_instance == "ENG" & id == 92 ~ "CHH",
                             panelapp_instance == "ENG" & id == 650 ~ "CHH_GMS",
                             panelapp_instance == "AUS" & id == 98 ~ "DBA",
                             (panelapp_instance == "ENG" & id == 544 | panelapp_instance == "AUS" & id == 78) ~ "CHO",
                             (panelapp_instance == "ENG" & id == 283 | panelapp_instance == "AUS" & id == 263) ~ "CKD",
                             (panelapp_instance == "ENG" & id == 9 | panelapp_instance == "AUS" & id == 99) ~ "DSD",
                             (panelapp_instance == "ENG" & id == 478 | panelapp_instance == "AUS" & id == 3763) ~ "FET",
                             (panelapp_instance == "ENG" & id == 312 | panelapp_instance == "AUS" & id == 3894) ~ "HOT",
                             panelapp_instance == "ENG" & id == 480 ~ "HPT",
                             (panelapp_instance == "ENG" & id == 285 | panelapp_instance == "AUS" & id == 250) ~ "INT",
                             (panelapp_instance == "ENG" & id == 112 | panelapp_instance == "AUS" & id == 203) ~ "MIT",
                             (panelapp_instance == "ENG" & id == 85 | panelapp_instance == "AUS" & id == 3120) ~ "NEU",
                             panelapp_instance == "ENG" & id == 482 ~ "RKT",
                             (panelapp_instance == "ENG" & id == 548 | panelapp_instance == "AUS" & id == 199) ~ "TKD"),
         source = "PanelApp", color = tolower(color),
         disease = ifelse(grepl("CHH", disease), "CHH", disease))
panel_app_combo <- panel_app_all %>%
  distinct(disease, Gene, .keep_all = TRUE)

## article count ####

# total (PubMed queries)
nrow(pmids_1_all)
length(unique(pmids_1_all$V1))

df_all_article <- df_all %>%
  filter(source != "REL") %>%
  left_join(pmids_1_df, by = c("PMID")) %>%
  mutate(disease = ifelse(disease.x == "", disease.y, disease.x)) %>%
  bind_rows(df_all %>%
              filter(source == "REL") %>%
              separate_longer_delim(disease, delim = ";")) %>%
  filter(Genes != "") %>%
  select(-disease.x, -disease.y)

df_all_article %>%
  count(source)

df_all_article %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, FTP, API, PTC, REL)

## gene count ####

# total all (unique)
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(Genes) %>%
  filter(!is.na(Genes)) %>%
  count()

gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & locus_group == "protein-coding gene") %>%
  count()

gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(Genes, .keep_all = TRUE) %>%
  filter(grepl("rs[0-9]|[0-9]+#[a-z]\\.|[0-9]+##[a-z]\\.|HGVS|RS#:", id_var)) %>%
  count()

# total by disease
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(disease, Genes) %>%
  filter(!is.na(Genes)) %>%
  count(disease)

gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(disease, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & locus_group == "protein-coding gene") %>%
  count(disease)

gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(disease, Genes, .keep_all = TRUE) %>%
  filter(grepl("rs[0-9]|[0-9]+#[a-z]\\.|[0-9]+##[a-z]\\.|HGVS|RS#:", id_var)) %>%
  count(disease)

# total by disease and source
gene_counts_per_disease <- gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, disease, Genes) %>%
  filter(!is.na(Genes)) %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, PT3FTP, PT3API, PTC, REL, REL2)
colSums(gene_counts_per_disease[,2:ncol(gene_counts_per_disease)])

gene_counts_per_disease_hgnc <- gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & locus_group == "protein-coding gene") %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, PT3FTP, PT3API, PTC, REL, REL2)
colSums(gene_counts_per_disease_hgnc[,2:ncol(gene_counts_per_disease_hgnc)])

gene_counts_per_disease_var <- gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  filter(grepl("rs[0-9]|[0-9]+#[a-z]\\.|[0-9]+##[a-z]\\.|HGVS|RS#:", id_var)) %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, PT3FTP, PT3API, PTC, REL, REL2)
colSums(gene_counts_per_disease_var[,2:ncol(gene_counts_per_disease_var)])

# total by source (unique)
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, Genes) %>%
  filter(!is.na(Genes)) %>%
  count(source)

gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & locus_group == "protein-coding gene") %>%
  count(source)

gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & grepl("rs[0-9]|[0-9]+#[a-z]\\.|[0-9]+##[a-z]\\.|HGVS|RS#:", id_var)) %>%
  count(source)

# All genes
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes)) %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, PT3FTP, PT3API, PTC, REL, REL2)
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & !is.na(n_var)) %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, PT3FTP, PT3API, PTC, REL, REL2)
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & (!is.na(n_var) | (is.na(n_var) & n_total <= 5)) & color != "") %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, PT3FTP, PT3API, PTC, REL, REL2)

# PanelApp panel genes only
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & color != "") %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, PT3FTP, PT3API, PTC, REL, REL2)
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & !is.na(n_var) & color != "") %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, PT3FTP, PT3API, PTC, REL, REL2)
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  filter(!is.na(Genes) & (!is.na(n_var) | (is.na(n_var) & n_total <= 5)) & color != "") %>%
  count(source, disease) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  select(disease, PT3FTP, PT3API, PTC, REL, REL2)

## gene recurrence ####

# all together
total_n <- gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  filter(!is.na(Genes)) %>%
  count()

gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  filter(!is.na(Genes)) %>%
  count(n_total) %>%
  mutate(pct_total = (n/total_n$n)*100)

# by source
source_n <- gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  filter(!is.na(Genes)) %>%
  count(source)

gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  filter(!is.na(Genes)) %>%
  count(source, n_total) %>%
  mutate(pct_total = case_when(source == "PT3API" ~ (n / subset(source_n, source_n$source == "PT3API")$n)*100,
                               source == "PT3FTP" ~ (n / subset(source_n, source_n$source == "PT3FTP")$n)*100,
                               source == "PTC" ~ (n / subset(source_n, source_n$source == "PTC")$n)*100,
                               source == "REL" ~ (n / subset(source_n, source_n$source == "REL")$n)*100,
                               source == "REL2" ~ (n / subset(source_n, source_n$source == "REL2")$n)*100))

# by source and disease
source_dis_n <- gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  filter(!is.na(Genes)) %>%
  count(source, disease)

gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  filter(!is.na(Genes) & source != "REL2") %>%
  count(source, disease, n_total) %>%
  left_join(source_dis_n %>%
              rename(count = n), by = c("source", "disease")) %>%
  mutate(pct_total = (n / count)*100,
         group = factor(case_when(n_total == 1 ~ "1",
                                  n_total == 2 ~ "2",
                                  n_total >= 3 & n_total <= 5 ~ "3-5",
                                  n_total > 5 & n_total <= 10 ~ "6-10",
                                  n_total > 10 & n_total <= 25 ~ "11-25",
                                  n_total > 25 & n_total <= 50 ~ "26-50",
                                  n_total > 50 ~ ">50"),
                        levels = c("1", "2", "3-5", "6-10", "11-25", "26-50", ">50")),
         disease = case_when(disease == "CHO" ~ "CS", disease == "CKD" ~ "CyKD", disease == "DBA" ~ "DBAS",
                             disease == "FET" ~ "FETAL", disease == "HOT" ~ "HypoPT", disease == "INT" ~ "ID",
                             disease == "MIT" ~ "MITO", disease == "NEU" ~ "NEUROP", disease == "RKT" ~ "HR",
                             TRUE ~ disease),
         source = factor(case_when(source == "PT3FTP" ~ "PT3 (FTP)",
                                   source == "PT3API" ~ "PT3 (API)",
                                   source == "PTC" ~ "PTC", 
                                   source == "REL" ~ "PT3 (Relations)"), 
                         levels = c("PT3 (FTP)", "PTC", "PT3 (API)", "PT3 (Relations)"))) %>%
  ggplot(aes(x = "", y = pct_total, fill = group)) +
  geom_col() + 
  coord_polar(theta = "y") +
  scale_fill_brewer() + 
  theme_void() + 
  facet_grid(source~disease) + 
  theme(legend.title = element_blank(), legend.position = "bottom",
        strip.text.x = element_text(size = 15, vjust = 0, margin = margin(b=5)), 
        strip.text.y = element_text(size = 15, hjust = 0),
        legend.text = element_text(size = 15))

ggsave(filename = "path/to/Gene_recurrence.jpg", 
       plot = last_plot(), device = "jpg", dpi = 400, height = 8, width = 24, units = "in")

## missing ####
missing_list_venn_df <- missing_list %>%
  separate_longer_delim(source, delim = ";") %>%
  mutate(value = TRUE,
         source = case_when(source == "PT3FTP" ~ "FTP", source == "PT3API" ~ "API", 
                            source != "PT3FTP" & source != "PT3API" ~ source)) %>%
  pivot_wider(names_from = source, values_from = value, values_fill = FALSE)

legend_plot_missing <- ggplot(data.frame(x = 1:4, y = 1, group = factor(c("PT3 (FTP)", "PT3 (API)", "PTC", "PT3 (Relations)"),
                                                                        levels = c("PT3 (FTP)", "PT3 (API)", "PTC", "PT3 (Relations)"))),
                              aes(x, y, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c("PT3 (FTP)" = "#0277BD", "PT3 (API)" = "#318E34", "PTC" = "#E1C420", "PT3 (Relations)" = "#E58EDD")) +
  labs(x = "", y = "") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank(), 
        legend.text = element_text(size = 15))

(missing_venn_plot <- missing_list_venn_df %>%
    ggplot(mapping = aes(A = FTP, B = API, C = PTC, D = REL)) +
    geom_venn(show_percentage = FALSE, 
              #auto_scale = TRUE, 
              fill_color = c("#0277BD", "#318E34", "#E1C420", "#E58EDD"),
              stroke_alpha = c(0, 0, 0, 0),
              text_size = 4, 
              set_name_size = 0, 
              digits = 0) + 
    facet_wrap(~disease, ncol = 6) +
    theme_void() +
    theme(strip.text = element_text(size = 15, vjust = 0), panel.spacing = unit(0, "line")))

missing_venn_plot + legend_plot_missing +
  plot_layout(widths = c(1, 0))

## sensitivity - Venn diagram ####
euler_data <- panel_app_all %>%
  mutate(source = "PanelApp", color = tolower(color), disease = ifelse(grepl("CHH", disease), "CHH", disease)) %>%
  rename(Genes = entrez_id) %>%
  select(disease, Genes, source) %>%
  bind_rows(gene_list %>%
              filter(panelapp_instance == "COMBO" | disease %in% c("CDH", "DBA", "HPT", "RKT")) %>%
              select(disease, Genes, source) %>%
              mutate(disease = ifelse(grepl("CHH", disease), "CHH", disease))) %>%
  mutate(disease = case_when(disease == "CHO" ~ "CS", disease == "CKD" ~ "CyKD", disease == "DBA" ~ "DBAS",
                             disease == "FET" ~ "FETAL", disease == "HOT" ~ "HypoPT", disease == "INT" ~ "ID",
                             disease == "MIT" ~ "MITO", disease == "NEU" ~ "NEUROP", disease == "RKT" ~ "HR",
                             TRUE ~ disease)) %>%
  filter(!is.na(Genes)) %>%
  distinct(source, disease, Genes) %>%
  group_by(disease) %>%
  group_split()

names(euler_data) <- euler_data %>% 
  group_keys() %>% 
  pull(disease)

# total counts
counts_panelapp <- sapply(euler_data, function(x) nrow(x[x$source == "PanelApp", ]))
counts_ftp <- sapply(euler_data, function(x) nrow(x[x$source == "PT3FTP", ]))
counts_api <- sapply(euler_data, function(x) nrow(x[x$source == "PT3API", ]))
counts_ptc <- sapply(euler_data, function(x) nrow(x[x$source == "PTC", ]))
counts_rel <- sapply(euler_data, function(x) nrow(x[x$source == "REL", ]))

# difference EXTRACT from PANELAPP
setdiff_ftp_panelapp <- sapply(euler_data, function(x) length(setdiff(x[x$source == "PT3FTP", ]$Genes, x[x$source == "PanelApp", ]$Genes)))
setdiff_api_panelapp <- sapply(euler_data, function(x) length(setdiff(x[x$source == "PT3API", ]$Genes, x[x$source == "PanelApp", ]$Genes)))
setdiff_ptc_panelapp <- sapply(euler_data, function(x) length(setdiff(x[x$source == "PTC", ]$Genes, x[x$source == "PanelApp", ]$Genes)))
setdiff_rel_panelapp <- sapply(euler_data, function(x) length(setdiff(x[x$source == "REL", ]$Genes, x[x$source == "PanelApp", ]$Genes)))

# difference PANELAPP from EXTRACT
setdiff_panelpp_ftp <- sapply(euler_data, function(x) length(setdiff(x[x$source == "PanelApp", ]$Genes, x[x$source == "PT3FTP", ]$Genes)))
setdiff_panelpp_api <- sapply(euler_data, function(x) length(setdiff(x[x$source == "PanelApp", ]$Genes, x[x$source == "PT3API", ]$Genes)))
setdiff_panelpp_ptc <- sapply(euler_data, function(x) length(setdiff(x[x$source == "PanelApp", ]$Genes, x[x$source == "PTC", ]$Genes)))
setdiff_panelpp_rel <- sapply(euler_data, function(x) length(setdiff(x[x$source == "PanelApp", ]$Genes, x[x$source == "REL", ]$Genes)))

# intersect EXTRACT and PANELAPP
intersect_panelpp_ftp <- lapply(euler_data, function(x) length(base::intersect(x[x$source == "PT3FTP", ]$Genes, x[x$source == "PanelApp", ]$Genes)))
intersect_panelpp_api <- lapply(euler_data, function(x) length(base::intersect(x[x$source == "PT3API", ]$Genes, x[x$source == "PanelApp", ]$Genes)))
intersect_panelpp_ptc <- lapply(euler_data, function(x) length(base::intersect(x[x$source == "PTC", ]$Genes, x[x$source == "PanelApp", ]$Genes)))
intersect_panelpp_rel <- lapply(euler_data, function(x) length(base::intersect(x[x$source == "REL", ]$Genes, x[x$source == "PanelApp", ]$Genes)))

panelapp_parts <- data.frame(disease = names(counts_panelapp), value = as.numeric(counts_panelapp), source = "PA")

parts_list <- list(
  setdiff_ftp_panelapp = setdiff_from_panelapp[[1]],
  setdiff_api_panelapp = setdiff_from_panelapp[[2]],
  setdiff_ptc_panelapp = setdiff_from_panelapp[[3]],
  setdiff_rel_panelapp = setdiff_from_panelapp[[4]],
  intersect_panelpp_ftp = intersect_with_panelapp[[1]],
  intersect_panelpp_api = intersect_with_panelapp[[2]],
  intersect_panelpp_ptc = intersect_with_panelapp[[3]],
  intersect_panelpp_rel = intersect_with_panelapp[[4]]
)

parts <- bind_rows(
  map2(parts_list, c("FTP", "API", "PTC", "REL", "FTP+PA", "API+PA", "PTC+PA", "REL+PA"), 
       ~ {.x %>% 
           bind_rows() %>%
           pivot_longer(everything(), names_to = "disease", values_to = "value") %>%
           mutate(source = .y)}),
  panelapp_parts)

parts <- parts %>%
  mutate(label.x = case_when(source == "PA" ~ 0, source == "REL" ~ 0, source == "REL+PA" ~ 0,
                             source == "FTP" ~ 3.5, source == "FTP+PA" ~ 1.75, source == "API" ~ -3.5,
                             source == "API+PA" ~ -1.75, source == "PTC" ~ 0, source == "PTC+PA" ~ 0),
         label.y = case_when(source == "PA" ~ 0, source == "REL" ~ 3.5, source == "REL+PA" ~ 1.75,
                             source == "FTP" ~ 0, source == "FTP+PA" ~ 0, source == "API" ~ 0,
                             source == "API+PA" ~ 0, source == "PTC" ~ -3.5, source == "PTC+PA" ~ -1.75),
         label.hjust = case_when(source == "PA" ~ 0.5, source == "REL" ~ 0.5, source == "REL+PA" ~ 0.5,
                                 source == "FTP" ~ 0, source == "FTP+PA" ~ 0.25, source == "API" ~ 1,
                                 source == "API+PA" ~ 0.75, source == "PTC" ~ 0.5, source == "PTC+PA" ~ 0.5),
         label.vjust = case_when(source == "PA" ~ 0.5, source == "REL" ~ -1, source == "REL+PA" ~ 0,
                                 source == "FTP" ~ 0.5, source == "FTP+PA" ~ 0.5, source == "API" ~ 0.5,
                                 source == "API+PA" ~ 0.5, source == "PTC" ~ 1.5, source == "PTC+PA" ~ 0.5))

parts_overlap <- parts %>%
  select(disease, source, value) %>%
  filter(!source %in% c("FTP", "API", "PTC", "REL")) %>%
  pivot_wider(names_from = source, values_from = value) %>%
  mutate(overlap_ftp = ifelse(PA == `FTP+PA`, "bold", "plain"), overlap_api = ifelse(PA == `API+PA`, "bold", "plain"),
         overlap_ptc = ifelse(PA == `PTC+PA`, "bold", "plain"), overlap_rel = ifelse(PA == `REL+PA`, "bold", "plain")) %>%
  pivot_longer(cols = starts_with("overlap"), names_to = "overlap", values_to = "bold") %>%
  mutate(source = case_when(overlap == "overlap_ftp" ~ "FTP+PA", overlap == "overlap_api" ~ "API+PA",
                            overlap == "overlap_ptc" ~ "PTC+PA", overlap == "overlap_rel" ~ "REL+PA"))

parts <- parts %>%
  left_join(parts_overlap %>% 
              select(bold, source, disease)) %>%
  mutate(bold = ifelse(is.na(bold), "plain", bold),
         bold = ifelse(source == "PA", "bold", bold))

pubtator_counts <- gene_list %>%
  filter(panelapp_instance == "COMBO" | disease %in% c("CDH", "DBA", "HPT", "RKT")) %>%
  select(disease, Genes, source) %>%
  mutate(disease = ifelse(grepl("CHH", disease), "CHH", disease)) %>%
  mutate(disease = case_when(disease == "CHO" ~ "CS", disease == "CKD" ~ "CyKD", disease == "DBA" ~ "DBAS",
                             disease == "FET" ~ "FETAL", disease == "HOT" ~ "HypoPT", disease == "INT" ~ "ID",
                             disease == "MIT" ~ "MITO", disease == "NEU" ~ "NEUROP", disease == "RKT" ~ "HR",
                             TRUE ~ disease)) %>%
  filter(!is.na(Genes)) %>%
  distinct(source, disease, Genes) %>%
  count(source, disease) %>%
  mutate(N = "N1")

intersect_counts <- rbindlist(list(PT3FTP = intersect_panelpp_ftp, PT3API = intersect_panelpp_api, 
                                   PTC = intersect_panelpp_ptc, REL = intersect_panelpp_rel),
                              idcol = "source") %>%
  pivot_longer(cols = -source, names_to = "disease", values_to = "n")  %>%
  mutate(N = "Noverlap")

panelapp_counts <- panel_app_all %>% 
  mutate(disease = case_when(disease == "CHO" ~ "CS", disease == "CKD" ~ "CyKD", disease == "DBA" ~ "DBAS",
                             disease == "FET" ~ "FETAL", disease == "HOT" ~ "HypoPT", disease == "INT" ~ "ID",
                             disease == "MIT" ~ "MITO", disease == "NEU" ~ "NEUROP", disease == "RKT" ~ "HR",
                             TRUE ~ disease)) %>%
  distinct(disease, entrez_id) %>% 
  count(disease) %>% 
  mutate(source = "PT3FTP")

panelapp_counts_all <- panelapp_counts %>%
  bind_rows(panelapp_counts %>%
              mutate(source = "PT3API")) %>%
  bind_rows(panelapp_counts %>%
              mutate(source = "PTC")) %>%
  bind_rows(panelapp_counts %>%
              mutate(source = "REL"))  %>%
  mutate(N = "N2")

venn_df <- pubtator_counts %>%
  bind_rows(panelapp_counts_all) %>%
  bind_rows(intersect_counts) %>%
  pivot_wider(names_from = N, values_from = n) %>%
  group_by(disease) %>%
  group_split()
names(venn_df) <- unique(panelapp_counts_all$disease)

params <- lapply(venn_df, function(x) list(top = list(N1 = x[x$source == "REL", ]$N1,
                                                      N2 = x[x$source == "REL", ]$N2,
                                                      Noverlap = x[x$source == "REL", ]$Noverlap),
                                           bottom = list(N1 = x[x$source == "PTC", ]$N1,
                                                         N2 = x[x$source == "PTC", ]$N2,
                                                         Noverlap = x[x$source == "PTC", ]$Noverlap),
                                           left = list(N1 = x[x$source == "PT3API", ]$N1,
                                                       N2 = x[x$source == "PT3API", ]$N2,
                                                       Noverlap = x[x$source == "PT3API", ]$Noverlap),
                                           right = list(N1 = x[x$source == "PT3FTP", ]$N1,
                                                        N2 = x[x$source == "PT3FTP", ]$N2,
                                                        Noverlap = x[x$source == "PT3FTP", ]$Noverlap)))

d_list <- list(
  # CDH CHH CHO CKD DBA DSD FET
  list(D = 10), list(D = 12.5), list(D = 25), list(D = 20), list(D = 10), list(D = 20), list(D = 65),
  # HPT RKT HOT INT MIT NEU TKD
  list(D = 10), list(D = 12.5), list(D = 10), list(D = 65), list(D = 45), list(D = 25), list(D = 10))
df_all <- Map(function(p, args) {do.call(build_plus_venn, c(list(p), args))}, params, d_list)

df_all <- lapply(names(df_all), function(nm) {
  df <- df_all[[nm]]$polygons
  df$disease <- nm
  list(polygons = df, labels = df_all[[nm]]$labels)
})
names(df_all) <- unique(panelapp_counts_all$disease)

venns <- lapply(df_all, plot_plus_venn)

legend_plot_venn <- ggplot(data.frame(x = 0, y = 0,
                                      group = factor(c("PanelApp", "PT3 (FTP)", "PT3 (API)", "PTC", "PT3 (Relations)"),
                                                     levels = c("PanelApp", "PT3 (FTP)", "PTC", "PT3 (API)", "PT3 (Relations)"))),
                           aes(x, y, fill = group), alpha = 0.4) +
  geom_col() +
  theme_void() + 
  scale_fill_manual(values = c("PanelApp" = "#007A82", "PT3 (FTP)" = "#E55B4A", "PT3 (API)" = "#318E34", "PTC" = "#E1C420", "PT3 (Relations)" = "#E58EDD")) +
  labs(x = "", y = "") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank(), 
        legend.position = "right", legend.text = element_text(size = 15))

wrap_plots(venns[[1]], venns[[2]], venns[[3]], venns[[4]], venns[[5]],
           venns[[6]], venns[[7]], venns[[8]], venns[[9]], venns[[10]],
           venns[[11]], venns[[12]], venns[[13]], venns[[14]], legend_plot_venn,
           ncol = 5, nrow = 3)

ggsave(filename = "path/to/venn_25-11-2025.jpg", 
       plot = last_plot(), device = "jpg", dpi = 400, height = 7.5, width = 11.25, units = "in")

## PubTator vs LitVar ####
gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease)) %>%
  filter(!is.na(n_var) & is.na(n_human) & is.na(n_ortho)) %>%
  select(source, disease, Genes, Symbol, n_var, n_total, n_human, n_ortho) %>%
  View()
count(source, disease)

pt_lv_df <- gene_list %>%
  full_join(panel_app_all %>%
              mutate(n_panelapp = 1) %>%
              distinct(disease, entrez_id, .keep_all = TRUE) %>%
              select(disease, entrez_id, n_panelapp),
            by = join_by(disease == disease, Genes == entrez_id)) %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         n_panelapp = ifelse(is.na(n_panelapp), 0, n_panelapp),
         n_pubtator = ifelse(is.na(n_human), 0, n_human),
         n_var = ifelse(is.na(n_var), 0, n_var),
         n_ortho = ifelse(is.na(n_ortho), 0, n_ortho)) %>%
  filter(!is.na(Genes)) %>%
  select(source, disease, Genes, n_panelapp, n_pubtator, n_var, n_ortho) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  rename(PanelApp = n_panelapp, PubTator = n_pubtator, LitVar = n_var, Orthologs = n_ortho) %>%
  mutate(PubTator = ifelse(PubTator >= 1, TRUE, FALSE),
         LitVar = ifelse(LitVar >= 1, TRUE, FALSE),
         Orthologs = ifelse(Orthologs >= 1, TRUE, FALSE),
         PanelApp = ifelse(PanelApp >= 1, TRUE, FALSE)) %>%
  mutate(source = ifelse(is.na(source), "PT3FTP;PTC;PT3API;REL", source)) %>%
  separate_longer_delim(cols = source, delim = ";") %>%
  mutate(disease = case_when(disease == "CHO" ~ "CS", disease == "CKD" ~ "CyKD", disease == "DBA" ~ "DBAS",
                             disease == "FET" ~ "FETAL", disease == "HOT" ~ "HypoPT", disease == "INT" ~ "ID",
                             disease == "MIT" ~ "MITO", disease == "NEU" ~ "NEUROP", disease == "RKT" ~ "HR",
                             TRUE ~ disease),
         source = factor(case_when(source == "PT3FTP" ~ "PT3 (FTP)",
                                   source == "PT3API" ~ "PT3 (API)",
                                   source == "PTC" ~ "PTC", 
                                   source == "REL" ~ "PT3 (Relations)"), 
                         levels = c("PT3 (FTP)", "PTC", "PT3 (API)", "PT3 (Relations)")))

pt_lv_plot <- pt_lv_df %>%
  ggplot(mapping = aes(A = PanelApp, B = PubTator, C = LitVar, D = Orthologs)) +
  geom_venn(show_percentage = FALSE, 
            fill_color = c("#007A82", "#112E51", "#FF964F", "#823FFF"),
            stroke_color = "#00363A",
            text_size = 2.5, 
            stroke_alpha = c(0.75, 0, 0, 0), 
            set_name_size = 0, 
            set_name_color = "white",
            digits = 0) + 
  facet_grid(source~disease) +
  theme_void() +
  theme(strip.text.x = element_text(size = 15, vjust = 0, margin = margin(b=5)), 
        strip.text.y = element_text(size = 15, hjust = 0))

legend_plot <- ggplot(data.frame(x = 1:4, y = 1, group = factor(c("PanelApp","PubTator","LitVar","Orthologs"),
                                                                levels = c("PanelApp","PubTator","LitVar","Orthologs"))),
                      aes(x, y, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c("PanelApp" = "#007A82", "PubTator" = "#112E51", "LitVar" = "#FF964F", "Orthologs" = "#823FFF")) +
  labs(x = "", y = "") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank(), 
        legend.position = "bottom", legend.text = element_text(size = 15))

wrap_plots(pt_lv_plot, legend_plot, heights = c(1, 0))

ggsave(filename = "path/to/Venn_LitVar_22-11-2025.jpg", 
       plot = last_plot(), device = "jpg", dpi = 400, height = 8, width = 24, units = "in")

## precision ####
stats_list %>%
  filter((panelapp == "COMBO" | disease %in% c("CDH", "DBA", "HPT", "RKT")) & stat == "precision") %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         source = factor(source, levels = c("PT3FTP", "PTC", "PT3API", "REL"))) %>%
  ggplot(aes(x = disease, y = value, fill = source)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  ylim(0, 1) + 
  theme_minimal() + 
  labs(y = "Precision (%)", x = "") + 
  theme(legend.title = element_blank())

## gene ranking ####
gene_list %>%
  filter(panelapp_instance == "COMBO" | disease %in% c("CDH", "DBA", "HPT", "RKT")) %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         disease = case_when(disease == "CHO" ~ "CS", disease == "CKD" ~ "CyKD", disease == "DBA" ~ "DBAS",
                             disease == "FET" ~ "FETAL", disease == "HOT" ~ "HypoPT", disease == "INT" ~ "ID",
                             disease == "MIT" ~ "MITO", disease == "NEU" ~ "NEUROP", disease == "RKT" ~ "HR",
                             TRUE ~ disease),
         source = factor(case_when(source == "PT3FTP" ~ "PT3 (FTP)",
                                   source == "PT3API" ~ "PT3 (API)",
                                   source == "PTC" ~ "PTC", 
                                   source == "REL" ~ "PT3 (Relations)"), 
                         levels = c("PT3 (FTP)", "PTC", "PT3 (API)", "PT3 (Relations)")))%>%
  group_by(source, disease) %>%
  mutate(rank = max(rank_total, na.rm = TRUE), rank_pct = (rank_total / rank) * 100) %>%
  ungroup() %>%
  filter(color %in% c("green", "amber", "red")) %>%
  ggplot(aes(x = rank_pct, fill = disease)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 25) +
  scale_x_continuous(breaks = c(0, 50, 100), labels = c("0", "50", "100")) +
  labs(x = "Normalized Rank", y = "Normalized count", fill = "Panel") +
  facet_grid(source~disease) + 
  theme_minimal() +
  theme(strip.text = element_text(size = 15), axis.title = element_text(size = 12.5), 
        axis.text = element_text(size = 10), legend.position = "none")

ggsave(filename = "path/to/rank_22-11-2025.jpg", 
       plot = last_plot(), device = "jpg", dpi = 400, height = 8, width = 24, units = "in")

## GSEA ####
gmt_file_gobp <- "path/to/c5.go.bp.v2025.1.Hs.entrez.gmt"
gmt_pathways_gobp <- gmtPathways(gmt_file_gobp)
length(unique(unlist(gmt_pathways_gobp)))

rank_and_name <- function(x, pa = FALSE) {
  x <- x[!is.na(x$n_total), ]
  n <- x$n_total
  if (pa) {
    names(n) <- x$entrez_id
  } else {
    names(n) <- x$Genes
  }
  return(n)
}

list_name <- function(x, pa = FALSE) {
  if (pa) {
    d <- unique(x$disease)
    return(d)
  } else {
    q <- unique(x$source)
    d <- unique(x$disease)
    return(paste(d, q, sep = "_"))
  }
}

run_fgseabp <- function(x, list) {
  message(paste0("Running ", x))
  return(fgsea(pathways = gmt_pathways_gobp, stats = list[[x]], scoreType = "pos", minSize = 10, maxSize = 100))
}

# combine gene counts from all queries and remove PanelApp genes
gene_list_gsea <- gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         color = ifelse(color == "", NA, color)) %>%
  group_by(source, disease, Genes) %>%
  mutate(remove = ifelse(all(is.na(color)), 0, 1)) %>% 
  ungroup() %>%
  filter(remove == 0)

# remove genes without entrez id
dim(gene_list_gsea)
gene_list_gsea <- gene_list_gsea %>%
  filter(!is.na(Genes))
dim(gene_list_gsea)

gene_list_gsea %>% 
  count(source, disease, Genes) %>%
  filter(n > 1)

# define ranked data for GSEA
count_stats <- gene_list_gsea %>%
  filter(panelapp_instance == "COMBO" | disease %in% c("CDH", "DBA", "HPT", "RKT")) %>%
  arrange(source, disease, desc(n_total)) %>%
  group_by(source, disease) %>%
  select(source, disease, Genes, n_total) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  group_by(source, disease) %>%
  group_split()

rank_stats <- lapply(count_stats, rank_and_name)
names(rank_stats) <- unlist(lapply(count_stats, list_name))

count_stats_pa <- panel_app_all %>% 
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         n_total = case_when(color == "Green" ~ 3, color == "Amber" ~ 2, color == "Red" ~ 1)) %>%
  arrange(disease, desc(n_total)) %>%
  group_by(disease) %>%
  select(disease, entrez_id, n_total) %>%
  distinct(disease, entrez_id, .keep_all = TRUE) %>%
  group_by(disease) %>%
  group_split()

rank_stats_pa <- lapply(count_stats_pa, rank_and_name, pa = TRUE)
names(rank_stats_pa) <- unlist(lapply(count_stats_pa, list_name, pa = TRUE))

# GO - BP
gsea_gobp <- lapply(names(rank_stats), run_fgseabp, list = rank_stats)
names(gsea_gobp) <- names(rank_stats)

gsea_gobp_pa <- lapply(names(rank_stats_pa), run_fgseabp, list = rank_stats_pa)
names(gsea_gobp_pa) <- names(rank_stats_pa)
gsea_gobp_all_pa <- rbindlist(gsea_gobp_pa, idcol = "name")

gsea_gobp_all <- rbindlist(gsea_gobp, idcol = "name")
gsea_gobp_all <- gsea_gobp_all %>%
  separate_wider_delim(cols = name, delim = "_", names = c("disease", "source")) %>%
  bind_rows(gsea_gobp_all_pa %>%
              rename(disease = name) %>%
              mutate(source = "PanelApp")) %>%
  arrange(source, disease, pval)

fwrite(gsea_gobp_all, sep = "\t",
       "path/to/GSEA_21-11-2025.txt")

# top 10 for each 
gsea_gobp_all_top10 <- gsea_gobp_all %>%
  group_by(source, disease) %>%
  slice(1:10)

fwrite(gsea_gobp_all_top10, sep = "\t",
       "path/to/GSEA_top10_21-11-2025.txt")


# stats

# total pathways by source and disease
gsea_gobp_all %>% 
  count(source, disease) %>%
  pivot_wider(names_from = source, values_from = n)

# significant pathways by source and disease
gsea_gobp_all %>% 
  filter(pval < 0.05) %>%
  count(source, disease) %>%
  pivot_wider(names_from = source, values_from = n) %>%
  select(PanelApp, PT3FTP, PT3API, PTC, REL)

# pathways
overlapp_pathways_pa <- gsea_gobp_all %>%
  select(source, disease, pathway, pval) %>%
  pivot_wider(names_from = source, values_from = pval) %>%
  select(disease, pathway, PanelApp, PT3FTP, PT3API, PTC, REL)

# overlaps with PanelApp
overlapp_pathways_pa %>%
  group_by(disease) %>%
  reframe(overlap_ftp = sum(!is.na(PanelApp) & !is.na(PT3FTP)), 
          overlap_api = sum(!is.na(PanelApp) & !is.na(PT3API)),
          overlap_ptc = sum(!is.na(PanelApp) & !is.na(PTC)), 
          overlap_rel = sum(!is.na(PanelApp) & !is.na(REL))) %>% 
  fwrite("path/to/table.txt", sep = "\t")

# overlaps with PanelApp where PanelApp and at least one group is significant
gsea_gobp_all %>% 
  filter(pval < 0.05) %>%
  select(source, disease, pathway, pval) %>%
  pivot_wider(names_from = source, values_from = pval) %>%
  filter(!is.na(PanelApp)) %>%
  select(disease, pathway, PanelApp, PT3FTP, PT3API, PTC, REL) %>%
  group_by(disease) %>%
  reframe(overlap_ftp = sum(!is.na(PanelApp) & !is.na(PT3FTP)), 
          overlap_api = sum(!is.na(PanelApp) & !is.na(PT3API)),
          overlap_ptc = sum(!is.na(PanelApp) & !is.na(PTC)), 
          overlap_rel = sum(!is.na(PanelApp) & !is.na(REL))) %>%
  fwrite("path/to/table.txt", sep = "\t")

overlapp_pathways_pa %>%
  filter(!is.na(PanelApp) & PanelApp < 0.05) %>%
  group_by(disease) %>%
  reframe(overlap_ftp = sum(!is.na(PanelApp) & !is.na(PT3FTP)), overlap_api = sum(!is.na(PanelApp) & !is.na(PT3API)),
          overlap_ptc = sum(!is.na(PanelApp) & !is.na(PT3FTP)), overlap_rel = sum(!is.na(PanelApp) & !is.na(REL))) %>% 
  fwrite("path/to/table.txt", sep = "\t")

overlapp_pathways_pa %>%
  filter(!is.na(PanelApp) & PanelApp < 0.05) %>%
  arrange(disease, PanelApp) %>%
  fwrite("path/to/GSEA_signif_PanelApp_21-11-2025.txt", sep = "\t")

### HGNC protein-coding ####

# combine gene counts from all queries and remove PanelApp genes
gene_list_gsea_hgnc <- gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         color = ifelse(color == "", NA, color)) %>%
  group_by(source, disease, Genes) %>%
  mutate(remove = ifelse(all(is.na(color)), 0, 1)) %>% 
  ungroup() %>%
  filter(remove == 0 & locus_group == "protein-coding gene")

# remove genes without entrez id
dim(gene_list_gsea_hgnc)
gene_list_gsea_hgnc <- gene_list_gsea_hgnc %>%
  filter(!is.na(Genes))
dim(gene_list_gsea_hgnc)

# define ranked data for GSEA
count_stats_hgnc <- gene_list_gsea_hgnc %>%
  filter(panelapp_instance == "COMBO" | disease %in% c("CDH", "DBA", "HPT", "RKT")) %>%
  arrange(source, disease, desc(n_total)) %>%
  group_by(source, disease) %>%
  select(source, disease, Genes, n_total) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  group_by(source, disease) %>%
  group_split()

rank_stats_hgnc <- lapply(count_stats_hgnc, rank_and_name)
names(rank_stats_hgnc) <- unlist(lapply(count_stats_hgnc, list_name))

gsea_gobp_hgnc <- lapply(names(rank_stats_hgnc), run_fgseabp, list = rank_stats_hgnc)
names(gsea_gobp_hgnc) <- names(rank_stats_hgnc)

gsea_gobp_hgnc_all <- rbindlist(gsea_gobp_hgnc, idcol = "name")
gsea_gobp_hgnc_all <- gsea_gobp_hgnc_all %>%
  separate_wider_delim(cols = name, delim = "_", names = c("disease", "source")) %>%
  bind_rows(gsea_gobp_all_pa %>%
              rename(disease = name) %>%
              mutate(source = "PanelApp")) %>%
  arrange(source, disease, pval)

fwrite(gsea_gobp_hgnc_all, sep = "\t",
       "path/to/GSEA_HGNC_21-11-2025.txt")

# top 10 for each 
gsea_gobp_hgnc_all_top10 <- gsea_gobp_hgnc_all %>%
  group_by(source, disease) %>%
  slice(1:10)

fwrite(gsea_gobp_hgnc_all_top10, sep = "\t",
       "path/to/GSEA_HGNC_top10_21-11-2025.txt")


# stats

# total pathways by source and disease
gsea_gobp_hgnc_all %>% 
  count(source, disease) %>%
  pivot_wider(names_from = source, values_from = n) %>%
  select(disease, PanelApp, PT3FTP, PT3API, PTC, REL)

# significant pathways by source and disease
gsea_gobp_hgnc_all_signif_pathways <- gsea_gobp_hgnc_all %>% 
  filter(pval < 0.05) %>%
  count(source, disease) %>%
  pivot_wider(names_from = source, values_from = n) %>%
  select(PanelApp, PT3FTP, PT3API, PTC, REL)

fwrite(gsea_gobp_hgnc_all_signif_pathways, "path/to/table.txt", sep = "\t")

# pathways
overlapp_pathways_pa_hgnc <- gsea_gobp_hgnc_all %>%
  select(source, disease, pathway, pval) %>%
  pivot_wider(names_from = source, values_from = pval) %>%
  select(disease, pathway, PanelApp, PT3FTP, PT3API, PTC, REL)

# overlaps with PanelApp
overlaps_hgnc_panelapp <- overlapp_pathways_pa_hgnc %>%
  group_by(disease) %>%
  reframe(overlap_ftp = sum(!is.na(PanelApp) & !is.na(PT3FTP)), 
          overlap_api = sum(!is.na(PanelApp) & !is.na(PT3API)),
          overlap_ptc = sum(!is.na(PanelApp) & !is.na(PTC)), 
          overlap_rel = sum(!is.na(PanelApp) & !is.na(REL)))

fwrite(overlaps_hgnc_panelapp, "path/to/table.txt", sep = "\t")

# overlaps with PanelApp where PanelApp and at least one group is significant
overlaps_hgnc_panelapp_signif <- gsea_gobp_hgnc_all %>% 
  filter(pval < 0.05) %>%
  select(source, disease, pathway, pval) %>%
  pivot_wider(names_from = source, values_from = pval) %>%
  filter(!is.na(PanelApp)) %>%
  select(disease, pathway, PanelApp, PT3FTP, PT3API, PTC, REL) %>%
  group_by(disease) %>%
  reframe(overlap_ftp = sum(!is.na(PanelApp) & !is.na(PT3FTP)), 
          overlap_api = sum(!is.na(PanelApp) & !is.na(PT3API)),
          overlap_ptc = sum(!is.na(PanelApp) & !is.na(PTC)), 
          overlap_rel = sum(!is.na(PanelApp) & !is.na(REL)))

fwrite("path/to/table.txt", overlaps_hgnc_panelapp_signif, sep = "\t")

overlapp_pathways_pa_hgnc %>%
  filter(!is.na(PanelApp) & PanelApp < 0.05) %>%
  group_by(disease) %>%
  reframe(overlap_ftp = sum(!is.na(PanelApp) & !is.na(PT3FTP)), overlap_api = sum(!is.na(PanelApp) & !is.na(PT3API)),
          overlap_ptc = sum(!is.na(PanelApp) & !is.na(PT3FTP)), overlap_rel = sum(!is.na(PanelApp) & !is.na(REL))) %>% 
  fwrite("path/to/table.txt", sep = "\t")

overlapp_pathways_pa_hgnc %>%
  filter(!is.na(PanelApp) & PanelApp < 0.05) %>%
  arrange(disease, PanelApp) %>%
  fwrite("path/to/GSEA_HGNC_signif_PanelApp_21-11-2025.txt", sep = "\t")

### Genes w/ variants only ####
# combine gene counts from all queries and remove PanelApp genes
gene_list_var <- gene_list %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         color = ifelse(color == "", NA, color)) %>%
  filter(grepl("rs[0-9]|[0-9]+#[a-z]\\.|[0-9]+##[a-z]\\.|HGVS|RS#:", id_var))

gene_list_gsea_var <- gene_list_var %>%
  mutate(disease = ifelse(grepl("^CHH", disease), "CHH", disease),
         color = ifelse(color == "", NA, color)) %>%
  group_by(source, disease, Genes) %>%
  mutate(remove = ifelse(all(is.na(color)), 0, 1)) %>% 
  ungroup() %>%
  filter(remove == 0)

gene_list_var %>% 
  filter(!is.na(color)) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  count(source, disease) %>%
  left_join(panel_app_combo %>% 
              count(disease) %>%
              rename(count_panelapp = n)) %>%
  mutate(tpr = n/count_panelapp) %>%
  select(disease, source, tpr) %>%
  pivot_wider(names_from = source, values_from = tpr)

gene_list_var %>% 
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  group_by(source, disease) %>%
  reframe(n = n(), n_present = sum(!is.na(color))) %>%
  mutate(precision = n_present / n) %>%
  select(disease, source, precision) %>%
  pivot_wider(names_from = source, values_from = precision)

gene_list %>% 
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  group_by(source, disease) %>%
  reframe(n = n(), n_present = sum(color != "")) %>%
  mutate(precision = n_present / n) %>%
  select(disease, source, precision) %>%
  pivot_wider(names_from = source, values_from = precision)

gene_list_gsea_var %>%
  count(source, disease) %>%
  pivot_wider(names_from = source, values_from = n)

gene_list_gsea %>%
  count(source, disease) %>%
  pivot_wider(names_from = source, values_from = n)

# remove genes without entrez id
dim(gene_list_gsea_var)
gene_list_gsea_var <- gene_list_gsea_var %>%
  filter(!is.na(Genes))
dim(gene_list_gsea_var)

# define ranked data for GSEA
count_stats_var <- gene_list_gsea_var %>%
  filter(panelapp_instance == "COMBO" | disease %in% c("CDH", "DBA", "HPT", "RKT")) %>%
  arrange(source, disease, desc(n_total)) %>%
  group_by(source, disease) %>%
  select(source, disease, Genes, n_total) %>%
  distinct(source, disease, Genes, .keep_all = TRUE) %>%
  group_by(source, disease) %>%
  group_split()

rank_stats_var <- lapply(count_stats_var, rank_and_name)
names(rank_stats_var) <- unlist(lapply(count_stats_var, list_name))

gsea_gobp_var <- lapply(names(rank_stats_var), run_fgseabp, list = rank_stats_var)
names(gsea_gobp_var) <- names(rank_stats_var)

gsea_gobp_var_all <- rbindlist(gsea_gobp_var, idcol = "name")
gsea_gobp_var_all <- gsea_gobp_var_all %>%
  separate_wider_delim(cols = name, delim = "_", names = c("disease", "source")) %>%
  bind_rows(gsea_gobp_all_pa %>%
              rename(disease = name) %>%
              mutate(source = "PanelApp")) %>%
  arrange(source, disease, pval)

fwrite(gsea_gobp_var_all, sep = "\t",
       "path/to/GSEA_var_21-11-2025.txt")

# top 10 for each 
gsea_gobp_var_all_top10 <- gsea_gobp_var_all %>%
  group_by(source, disease) %>%
  slice(1:10)

fwrite(gsea_gobp_var_all_top10, sep = "\t",
       "path/to/GSEA_var_top10_21-11-2025.txt")


# stats

# total pathways by source and disease
gsea_gobp_var_all %>% 
  count(source, disease) %>%
  pivot_wider(names_from = source, values_from = n) %>%
  select(disease, PanelApp, PT3FTP, PT3API, PTC, REL) %>%
  mutate(disease = case_when(disease == "RKT" ~ "HR", disease == "HOT" ~ "HypoPT", TRUE ~ disease)) %>%
  arrange(disease)

# significant pathways by source and disease
gsea_gobp_var_all %>% 
  filter(pval < 0.05) %>%
  count(source, disease) %>%
  pivot_wider(names_from = source, values_from = n) %>%
  select(disease, PanelApp, PT3FTP, PT3API, PTC, REL)  %>%
  mutate(disease = case_when(disease == "RKT" ~ "HR", disease == "HOT" ~ "HypoPT", TRUE ~ disease)) %>%
  arrange(disease)

# pathways
overlapp_pathways_pa_var <- gsea_gobp_var_all %>%
  select(source, disease, pathway, pval) %>%
  pivot_wider(names_from = source, values_from = pval) %>%
  select(disease, pathway, PanelApp, PT3FTP, PT3API, PTC, REL) %>%
  mutate(disease = case_when(disease == "RKT" ~ "HR", disease == "HOT" ~ "HypoPT", TRUE ~ disease)) %>%
  arrange(disease)

# overlaps with PanelApp
overlapp_pathways_pa_var %>%
  group_by(disease) %>%
  reframe(overlap_ftp = sum(!is.na(PanelApp) & !is.na(PT3FTP)), 
          overlap_api = sum(!is.na(PanelApp) & !is.na(PT3API)),
          overlap_ptc = sum(!is.na(PanelApp) & !is.na(PTC)), 
          overlap_rel = sum(!is.na(PanelApp) & !is.na(REL))) %>% 
  fwrite("path/to/table.txt", sep = "\t")

# overlaps with PanelApp where PanelApp and at least one group is significant
gsea_gobp_var_all %>% 
  filter(pval < 0.05) %>%
  select(source, disease, pathway, pval) %>%
  pivot_wider(names_from = source, values_from = pval) %>%
  filter(!is.na(PanelApp)) %>%
  select(disease, pathway, PanelApp, PT3FTP, PT3API, PTC, REL) %>%
  group_by(disease) %>%
  reframe(overlap_ftp = sum(!is.na(PanelApp) & !is.na(PT3FTP)),
          overlap_api = sum(!is.na(PanelApp) & !is.na(PT3API)),
          overlap_ptc = sum(!is.na(PanelApp) & !is.na(PTC)), 
          overlap_rel = sum(!is.na(PanelApp) & !is.na(REL))) %>%
  fwrite("path/to/table.txt", sep = "\t")

overlapp_pathways_pa_var %>%
  filter(!is.na(PanelApp) & PanelApp < 0.05) %>%
  group_by(disease) %>%
  reframe(overlap_ftp = sum(!is.na(PanelApp) & !is.na(PT3FTP)), overlap_api = sum(!is.na(PanelApp) & !is.na(PT3API)),
          overlap_ptc = sum(!is.na(PanelApp) & !is.na(PT3FTP)), overlap_rel = sum(!is.na(PanelApp) & !is.na(REL))) %>% 
  fwrite("path/to/table.txt", sep = "\t")

overlapp_pathways_pa_var %>%
  filter(!is.na(PanelApp) & PanelApp < 0.05) %>%
  arrange(disease, PanelApp) %>%
  fwrite("path/to/GSEA_var_signif_PanelApp_21-11-2025.txt", sep = "\t")



overlapp_pathways_All <- overlapp_pathways_pa %>%
  filter(!is.na(PanelApp) & PanelApp < 0.05) %>%
  arrange(disease, PanelApp) %>%
  full_join(overlapp_pathways_pa_hgnc %>%
              filter(!is.na(PanelApp) & PanelApp < 0.05) %>%
              arrange(disease, PanelApp) %>%
              rename(PanelApp_HGNC = PanelApp, PT3FTP_HGNC = PT3FTP, PT3API_HGNC = PT3API, PTC_HGNC = PTC, REL_HGNC = REL),
            by = c("disease", "pathway")) %>%
  full_join(overlapp_pathways_pa_var %>%
              filter(!is.na(PanelApp) & PanelApp < 0.05) %>%
              arrange(disease, PanelApp) %>%
              rename(PanelApp_VAR = PanelApp, PT3FTP_VAR = PT3FTP, PT3API_VAR = PT3API, PTC_VAR = PTC, REL_VAR = REL),
            by = c("disease", "pathway")) %>%
  select(-PanelApp_HGNC, -PanelApp_VAR)

fwrite(overlapp_pathways_All, sep = "\t",
       "path/to/table_S12.txt")


