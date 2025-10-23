# ANALYSIS ####
## READ IN DATA ####
stats_list_1 <- fread("path/to/Query1_stats_df.txt", sep = "\t")
stats_list_2 <- fread("path/to/Query2_stats_df.txt", sep = "\t")
stats_list_PT <- fread("path/to/QueryPT1_stats.txt", sep = "\t")

gene_list_1 <- fread("path/to/Query1_genes_df.txt", sep = "\t")
gene_list_2 <- fread("path/to/Query2_genes_df.txt", sep = "\t")
gene_list_PT <- fread("path/to/QueryPT1_genes_df.txt", sep = "\t")

missing_list_1 <- fread("path/to/Query1_missing_df.txt", sep = "\t")
missing_list_2 <- fread("path/to/Query2_missing_df.txt", sep = "\t")
missing_list_PT <- fread("path/to/QueryPT1_missing_df.txt", sep = "\t")

complete_query1 <- stats_list_1 %>% 
  filter(grepl("^tpr$", stat)) %>% 
  filter(value == 1) %>%
  mutate(keep = paste(disease, panelapp, sep = "_"))

# PanelApp
panel_app_all <- fread("path/to/PanelAppENG-AUS_panels_PUB-EDIT_NO-DUPLICATES.txt", sep = "\t")
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
                             (panelapp_instance == "ENG" & id == 548 | panelapp_instance == "AUS" & id == 199) ~ "TKD"))
  
## Article count ####
stats_list_compare_art <- stats_list_1 %>%
  filter(grepl("^n_art", stat)) %>% 
  distinct(value, .keep_all = TRUE) %>%
  left_join(stats_list_2 %>% 
              filter(grepl("^n_art", stat)) %>%
              distinct(value, .keep_all = TRUE) %>%
              rename(value2 = value),
            by = join_by(stat, disease, panelapp)) %>%
  left_join(stats_list_PT %>% 
              filter(grepl("^n_art", stat)) %>%
              distinct(value, .keep_all = TRUE) %>%
              rename(valuePT = value),
            by = join_by(stat, disease, panelapp)) %>%
  pivot_longer(cols = c(value, value2, valuePT), names_to = "query", values_to = "value") %>%
  mutate(query = gsub("value", "query", query))

(article_plot <- stats_list_compare_art %>%
  mutate(query = factor(case_when(query == "query" ~ "First",
                                  query == "query2" ~ "Updated",
                                  query == "queryPT" ~ "Relations"),
                               levels = c("First", "Updated", "Relations"))) %>%
  ggplot(aes(x = disease, y = value, fill = query, group = query)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  # geom_text(aes(x = disease, y = value, group = query, label = round(value, 2)),
  #           position = position_dodge(width = 0.9), size = 2.25, vjust = -0.25, inherit.aes = FALSE) + 
  scale_fill_manual(values = c("#E55B4A", "#E1C420", "#E58EDD"),
                    guide = guide_legend(override.aes = list(pattern = 'none'))) +
  scale_y_continuous(breaks = seq(0, 900000, by = 100000)) +
  labs(x = "", y = "Article count", fill = "Query", tag = "A") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12.5), plot.margin = unit(c(0,0,0,0), "mm")))

## Gene count ####
stats_list_compare_gene <- stats_list_1 %>%
  filter(grepl("^n_gene", stat)) %>% 
  distinct(value, .keep_all = TRUE) %>%
  left_join(stats_list_2 %>% 
                     filter(grepl("^n_gene", stat)) %>%
                     distinct(value, .keep_all = TRUE) %>%
                     rename(value2 = value),
                   by = join_by(stat, disease, panelapp)) %>%
  left_join(stats_list_PT %>% 
                     filter(grepl("^n_gene", stat)) %>%
                     distinct(value, .keep_all = TRUE) %>%
                     rename(valuePT = value),
                   by = join_by(stat, disease, panelapp)) %>%
  pivot_longer(cols = c(value, value2, valuePT), names_to = "query", values_to = "value") %>%
  mutate(query = gsub("value", "query", query))

(gene_plot <- stats_list_compare_gene %>%
    mutate(query = factor(case_when(query == "query" ~ "First",
                                                  query == "query2" ~ "Updated",
                                                  query == "queryPT" ~ "Relations"),
                                 levels = c("First", "Updated", "Relations"))) %>%
    ggplot(aes(x = disease, y = value, fill = query, group = query)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("#E55B4A", "#E1C420", "#E58EDD"),
                      guide = guide_legend(override.aes = list(pattern = 'none'))) +
    # geom_text(aes(x = disease, y = value, group = query, label = round(value, 2)),
    #           position = position_dodge(width = 0.9), size = 2.25, vjust = -0.25, inherit.aes = FALSE) + 
    scale_y_continuous(breaks = seq(0, 35000, by = 5000)) +
    labs(x = "", y = "Gene count", fill = "Query", tag = "B") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12.5), plot.margin = unit(c(0,0,0,0), "mm")))

# ggsave(filename = "path/to/Article_Gene_Precision_Queries1-2-3.jpg", 
#        plot = last_plot(), device = "jpg", dpi = 400, height = 4, width = 14, units = "in")

### Article + Gene plot ####
article_gene_plot <- article_plot + gene_plot +
  plot_layout(guides = "collect")

ggsave(filename = "path/to/Paper_gene_plot_16-10-2025.jpg",
       plot = last_plot(), device = "jpg", dpi = 400, height = 4, width = 11.5, units = "in")

## Missing ####
missing_all <- missing_list_1 %>%
  mutate(query = "query1") %>%
  bind_rows(missing_list_2 %>%
              mutate(query = "query2")) %>%
  bind_rows(missing_list_PT %>%
              mutate(query = "queryPT"))
  
missing_all %>%
  distinct(query, disease, missing)

missing_all %>%
  distinct(query, disease, missing) %>%
  count(query, disease) %>%
  pivot_wider(names_from = query, values_from = n) %>%
  arrange(disease)

# missing Relations query, not PubMed query
missing_all %>%
  distinct(query, disease, missing) %>%
  mutate(n = 1) %>%
  pivot_wider(names_from = query, values_from = n, values_fill = 0) %>%
  arrange(disease) %>%
  filter(queryPT > 0 & (query1== 0 | query2 == 0)) %>%
  count(disease)

missing_relations_not_pubmed <- missing_all %>%
  distinct(query, disease, missing) %>%
  mutate(n = 1) %>%
  pivot_wider(names_from = query, values_from = n, values_fill = 0) %>%
  arrange(disease) %>%
  filter(queryPT > 0 & (query1== 0 | query2 == 0)) %>%
  mutate(context = paste(disease, missing, sep = "_"))

missing_relations_not_pubmed_PMID <- gene_list_2 %>%
  bind_rows(gene_list_1 %>%
              mutate(context = paste(disease, panelapp_instance, sep = "_")) %>%
              filter(context %in% complete_query1$keep)) %>%
  mutate(context = paste(disease, Gene, sep = "_")) %>%
  filter(context %in% missing_relations_not_pubmed$context) %>%
  rowwise() %>%
  mutate(PMIDs = paste_unique(c(PMID, PMID_ortho, PMID_var))) %>%
  distinct(disease, Gene, PMIDs) %>%
  mutate(PMIDs = gsub("^;|;$", "", PMIDs), PMIDs = gsub(";;", "", PMIDs))

fwrite(missing_relations_not_pubmed_PMID, sep = "\t", file = "path/to/missingRE_notPubMed.txt")

missing_relations_pmids <- missing_relations_not_pubmed_PMID %>%
  separate_longer_delim(PMIDs, delim = ";") %>%
  filter(PMIDs != "") %>%
  distinct(disease, PMIDs) %>%
  pull(PMIDs)
missing_relations_pubtator <- fread("path/to/PT_Master.txt", sep = "\t") %>%
  filter(PMID %in% missing_relations_pmids)

missing_relations_not_pubmed_panelapp <- missing_relations_not_pubmed_PMID %>%
  separate_longer_delim(PMIDs, delim = ";") %>%
  filter(PMIDs != "") %>%
  distinct(disease, PMIDs, .keep_all = TRUE) %>%
  left_join(missing_relations_pubtator, by = join_by(PMIDs == PMID)) %>%
  select(disease, Gene, PMIDs, Relations) %>%
  separate_longer_delim(Relations, delim = "~") %>%
  rowwise() %>%
  filter(grepl(Gene, Relations) | Relations == "") %>%
  group_by(disease, Gene) %>%
  reframe(PMIDs = paste_unique(PMIDs), Relations = paste_unique(Relations, c = "~")) %>%
  left_join(panel_app_all[,c("disease", "Gene", "publications2")], by = c("disease", "Gene")) %>%
  mutate(publications2 = ifelse(is.na(publications2), "", publications2)) %>%
  select(disease, Gene, PMIDs, Relations, publications2) %>%
  group_by(disease, Gene) %>%
  reframe(PMIDs = paste_unique(PMIDs), Relations = paste_unique(Relations, c = "~"),
          publications2 = paste_unique(publications2)) %>%
  mutate(Relations = gsub("^~|~$", "", Relations)) %>%
  rowwise() %>%
  mutate(PMIDs_PM_PA = paste_unique(c(unlist(strsplit(PMIDs, ";")), unlist(strsplit(publications2, ";"))))) %>%
  distinct(disease, Gene, PMIDs_PM_PA, Relations)

fwrite(missing_relations_not_pubmed_panelapp, sep = "\t", file = "path/to/missingRE_notPubMed_PanelApp.txt")

## Sensitivity ####
### Barplot ####
stats_list_compare <- stats_list_1 %>%
  filter(grepl("^tpr", stat)) %>% 
  left_join(stats_list_2 %>%
              filter(grepl("^tpr", stat)) %>%
              dplyr::rename(value2 = value),
            by = join_by(stat, disease, panelapp)) %>%
  pivot_longer(cols = c(value, value2), names_to = "query", values_to = "value") %>%
  mutate(query = gsub("value", "query", query),
         panelapp = if_else(disease == "CHH_GMS", "ENG", panelapp),
         panelapp = if_else(panelapp == "AUS", "Australia", "Genomics England"),
         disease = if_else(disease == "CHH_GMS", "CHH\n(GMS)", disease))

(tpr_plot <- stats_list_compare %>%
  filter(grepl("tpr", stat)) %>%
  mutate(stat = factor(case_when(stat == "tpr" ~ "All",
                                 stat == "tpr_green" ~ "Green",
                                 stat == "tpr_amber" ~ "Amber",
                                 stat == "tpr_red" ~ "Red"),
                       levels = c("All", "Green", "Amber", "Red")),
         panelapp = factor(panelapp),
         query = case_when(query == "query" ~ "First", query == "query2" ~ "Updated")) %>%
  ggplot(aes(x = disease, y = value, fill = stat, alpha = query, group = stat)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(mapping = aes(x = disease, y = value, group = stat, alpha = query, label = round(value, 2)),
            inherit.aes = FALSE, show.legend = FALSE, size = 2.5, vjust = -0.25,
            position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#515151", "#4EA72E", "#FFBF00", "#F8766D"),
                    guide = guide_legend(override.aes = list(pattern = 'none'))) +
  scale_alpha_manual(values = c(1, 0.5)) +
  facet_wrap(~panelapp, ncol = 1) +
  labs(x = "", y = "Sensitivity", fill = "PanelApp rank", alpha = "Query", title = "PubMed queries") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12.5), plot.margin = unit(c(0,0,0,0), "mm")))

stats_list_compare_PT <- stats_list_2 %>%
  bind_rows(stats_list_1 %>%
              mutate(keep = paste(disease, panelapp, sep = "_")) %>%
              filter(keep %in% complete_query1$keep)) %>%
  select(-keep) %>%
  filter(grepl("^tpr", stat)) %>% 
  left_join(stats_list_PT %>%
              filter(grepl("^tpr", stat)) %>%
              dplyr::rename(valuePT = value),
            by = join_by(stat, disease, panelapp)) %>%
  pivot_longer(cols = c(value, valuePT), names_to = "query", values_to = "value") %>%
  mutate(query = gsub("value", "query", query),
        panelapp = if_else(disease == "CHH_GMS", "ENG", panelapp),
        panelapp = if_else(panelapp == "AUS", "Australia", "Genomics England"),
        disease = if_else(disease == "CHH_GMS", "CHH\n(GMS)", disease))

(tpr_plot_PT <- stats_list_compare_PT %>%
  filter(grepl("tpr", stat)) %>%
  filter(query == "queryPT") %>%
  mutate(stat = factor(case_when(stat == "tpr" ~ "All",
                                 stat == "tpr_green" ~ "Green",
                                 stat == "tpr_amber" ~ "Amber",
                                 stat == "tpr_red" ~ "Red"),
                              levels = c("All", "Green", "Amber", "Red")),
                panelapp = factor(panelapp),
                query = case_when(query == "query" ~ "PubMed query", 
                                  query == "queryPT" ~ "PubTator Relations")) %>%
  
  ggplot(aes(x = disease, y = value, fill = stat, group = stat)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(mapping = aes(x = disease, y = value, group = stat, label = round(value, 2)),
            inherit.aes = FALSE, show.legend = FALSE, size = 2.5, vjust = -0.25,
            position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#515151", "#4EA72E", "#FFBF00", "#F8766D"),
                    guide = guide_legend(override.aes = list(pattern = 'none'))) +
  scale_alpha_manual(values = c(1, 0.5, 0.25)) +
  facet_wrap(~panelapp, ncol = 1) +
  labs(x = "", y = "Sensitivity", fill = "PanelApp rank", title = "Relations query") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12.5), plot.margin = unit(c(0,0,0,0), "mm")))

tprs_plot <- (tpr_plot / tpr_plot_PT) +
  plot_layout(guides = "collect")

ggsave(filename = "path/to/TPR_queries_16-10-2025.jpg",
       plot = last_plot(), device = "jpg", dpi = 400, height = 12, width = 18, units = "in")

### Venn diagram ####
# ggvenn
tpr_venn_df <- panel_app_all %>%
  mutate(data = "PanelApp", color = tolower(color)) %>%
  full_join(gene_list_1 %>%
              mutate(data = "First", 
                     panelapp_instance = if_else(disease == "CHH_GMS", "ENG", panelapp_instance)) %>%
              distinct(disease, Gene, .keep_all = TRUE)) %>%
  full_join(gene_list_2 %>%
              bind_rows(gene_list_1 %>%
                          mutate(context = paste(disease, panelapp_instance, sep = "_")) %>%
                          filter(context %in% complete_query1$keep)) %>%
              mutate(data = "Updated", 
                     panelapp_instance = if_else(disease == "CHH_GMS", "ENG", panelapp_instance)) %>%
              distinct(disease, Gene, .keep_all = TRUE)) %>%
  full_join(gene_list_PT %>%
              mutate(data = "Relations", 
                     panelapp_instance = if_else(disease == "CHH_GMS", "ENG", panelapp_instance)) %>%
              distinct(disease, Gene, .keep_all = TRUE)) %>%
  mutate(disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%
  filter(Gene != "") %>%
  distinct(data, disease, Gene) %>%
  count(data, disease, Gene) %>%
  pivot_wider(names_from = data, values_from = n, values_fill = 0) %>%
  mutate(First = if_else(First >= 1, TRUE, FALSE),
         Updated = if_else(Updated >= 1, TRUE, FALSE),
         Relations = if_else(Relations >= 1, TRUE, FALSE),
         PanelApp = if_else(PanelApp >= 1, TRUE, FALSE))

#### 4-way ####
legend_plot <- ggplot(data.frame(x = 1:4, y = 1, group = factor(c("PanelApp","First","Updated","Relations"),
                                                                levels = c("PanelApp","First","Updated","Relations"))),
                      aes(x, y, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c("PanelApp" = "#007A82", "First" = "#E55B4A", "Updated" = "#E1C420", "Relations" = "#E58EDD")) +
  labs(x = "", y = "") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank())

(tpr_venn_plot_left <- tpr_venn_df %>%
    # filter(grepl("^[A-F]", disease)) %>%
    filter(disease %in% c("CDH", "CHO", "DSD", "HOT", "INT", "NEU", "TKD")) %>%
    ggplot(mapping = aes(A = PanelApp, B = First, C = Updated, D = Relations)) +
    geom_venn(show_percentage = FALSE, 
              #fill_color = c("#007A82", "#5688A5", "#AC95C7", "#E58EDD"),
              fill_color = c("#007A82", "#E55B4A", "#E1C420", "#E58EDD"),
              #auto_scale = TRUE, 
              stroke_color = "#00363A",
              text_size = 4, 
              stroke_alpha = c(0.75, 0, 0, 0), 
              set_name_size = 0, 
              digits = 0) + 
    facet_wrap(~disease, ncol = 1) +
    theme_void() +
    theme(strip.text = element_text(size = 15, vjust = 0), panel.spacing = unit(0, "line")))

(tpr_venn_plot_right <- tpr_venn_df %>%
    # filter(grepl("^[A-F]", disease)) %>%
    filter(disease %in% c("CHH", "CKD", "FET", "HPT", "MIT", "RKT")) %>%
    ggplot(mapping = aes(A = PanelApp, B = First, C = Updated, D = Relations)) +
    geom_venn(show_percentage = FALSE, 
              #fill_color = c("#007A82", "#5688A5", "#AC95C7", "#E58EDD"),
              fill_color = c("#007A82", "#E55B4A", "#E1C420", "#E58EDD"),
              #auto_scale = TRUE, 
              stroke_color = "#00363A",
              text_size = 4, 
              stroke_alpha = c(0.75, 0, 0, 0), 
              set_name_size = 0, 
              digits = 0) + 
    facet_wrap(~disease, ncol = 1) +
    theme_void() +
    theme(strip.text = element_text(size = 15, vjust = 0), panel.spacing = unit(0, "line")))

tpr_venn_plot_left 
tpr_venn_plot_right_full <- (plot_spacer() / tpr_venn_plot_right / plot_spacer() + plot_layout(heights = c(1, 50, 2)))

tpr_venn_plot_left + tpr_venn_plot_right_full + legend_plot +
  plot_layout(widths = c(1, 1, 0))

ggsave(filename = "path/to/Venn_PanelApp-combined_Queries1-2-3.jpg", 
       plot = last_plot(), device = "jpg", dpi = 400, height = 12, width = 7, units = "in")

#### 3 separate ####
legend_plot <- ggplot(data.frame(x = 1:4, y = 1, group = factor(c("PanelApp","First","Updated","Relations"),
                                                                levels = c("PanelApp","First","Updated","Relations"))),
                      aes(x, y, fill = group)) +
  geom_col() +
  theme_void() +
  scale_fill_manual(values = c("PanelApp" = "#007A82", "First" = "#E55B4A", "Updated" = "#E1C420", "Relations" = "#E58EDD")) +
  labs(x = "", y = "") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank())

circle_text_size = 3
(plotA <- tpr_venn_df %>%
    filter(disease %in% c("CDH", "CHH", "CHO", "CKD", "DBA", "DSD", "FET")) %>%
    filter(PanelApp == TRUE | First == TRUE) %>%
    ggplot(mapping = aes(A = PanelApp, B = First)) +
    geom_venn(show_percentage = FALSE, 
              hjust = 0,
              fill_color = c("#007A82", "#E55B4A"),
              position = position_dodge(2.5),
              text_size = circle_text_size, 
              stroke_alpha = c(0.75, 0), 
              stroke_color = "#00363A",
              stroke_size = 0.75,
              set_name_color = "white") + 
    facet_wrap(~disease, ncol = 13) +
    theme_void() +
    theme(strip.text = element_text(size = 15, vjust = 0), panel.spacing = unit(0, "line")))

(plotB <- tpr_venn_df %>%
    filter(disease %in% c("CDH", "CHH", "CHO", "CKD", "DBA", "DSD", "FET")) %>%
    filter(PanelApp == TRUE | Updated == TRUE) %>%
    ggplot(mapping = aes(A = PanelApp, B = Updated)) +
    geom_venn(show_percentage = FALSE, 
              fill_color = c("#007A82", "#E1C420"),
              position = position_dodge(2.5),
              text_size = circle_text_size, 
              stroke_alpha = c(0.75, 0), 
              stroke_color = "#00363A",
              stroke_size = 0.75,
              set_name_color = "white") + 
    facet_wrap(~disease, ncol = 13) +
    theme_void() +
    theme(strip.background = element_blank(), strip.text = element_blank(), 
          panel.spacing = unit(0, "line")))

(plotC <- tpr_venn_df %>%
    filter(disease %in% c("CDH", "CHH", "CHO", "CKD", "DBA", "DSD", "FET")) %>%
    filter(PanelApp == TRUE | Relations == TRUE) %>%
    ggplot(mapping = aes(A = PanelApp, B = Relations)) +
    geom_venn(show_percentage = FALSE, 
              fill_color = c("#007A82", "#E58EDD"),
              position = position_dodge(2.5),
              text_size = circle_text_size, 
              stroke_alpha = c(0.75, 0), 
              stroke_color = "#00363A",
              stroke_size = 0.75,
              set_name_color = "white") + 
    facet_wrap(~disease, ncol = 13) +
    theme_void() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          panel.spacing = unit(0, "line")))

(plotD <- tpr_venn_df %>%
    filter(disease %in% c("HOT", "HPT", "INT", "MIT", "NEU", "RKT", "TKD")) %>%
    filter(PanelApp == TRUE | First == TRUE) %>%
    ggplot(mapping = aes(A = PanelApp, B = First)) +
    geom_venn(show_percentage = FALSE, 
              fill_color = c("#007A82", "#E55B4A"),
              position = position_dodge(2.5),
              text_size = circle_text_size, 
              stroke_alpha = c(0.75, 0), 
              stroke_color = "#00363A",
              stroke_size = 0.75,
              set_name_color = "white") + 
    facet_wrap(~disease, ncol = 13) +
    theme_void() +
    theme(strip.text = element_text(size = 15, vjust = 0), panel.spacing = unit(0, "line")))

(plotE <- tpr_venn_df %>%
    filter(disease %in% c("HOT", "HPT", "INT", "MIT", "NEU", "RKT", "TKD")) %>%
    filter(PanelApp == TRUE | Updated == TRUE) %>%
    ggplot(mapping = aes(A = PanelApp, B = Updated)) +
    geom_venn(show_percentage = FALSE, 
              fill_color = c("#007A82", "#E1C420"),
              position = position_dodge(2.5),
              text_size = circle_text_size, 
              stroke_alpha = c(0.75, 0), 
              stroke_color = "#00363A",
              stroke_size = 0.75,
              set_name_color = "white") + 
    facet_wrap(~disease, ncol = 13) +
    theme_void() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          panel.spacing = unit(0, "line")))

(plotF <- tpr_venn_df %>%
    filter(disease %in% c("HOT", "HPT", "INT", "MIT", "NEU", "RKT", "TKD")) %>%
    filter(PanelApp == TRUE | Relations == TRUE) %>%
    ggplot(mapping = aes(A = PanelApp, B = Relations)) +
    geom_venn(show_percentage = FALSE, 
              fill_color = c("#007A82", "#E58EDD"),
              position = position_dodge(2.5),
              text_size = circle_text_size, 
              stroke_alpha = c(0.75, 0), 
              stroke_color = "#00363A",
              stroke_size = 0.75,
              set_name_color = "white") + 
    facet_wrap(~disease, ncol = 13) +
    theme_void() +
    theme(strip.background = element_blank(), strip.text = element_blank(), 
          panel.spacing = unit(0, "line")))

venns <- plotA / plotB / plotC / plotD / plotE / plotF

(final <- wrap_plots(venns, legend_plot, widths = c(1, 0)))

ggsave(filename = "path/to/Venn_PanelApp-separate_Queries_16-10-2025.jpg", 
       plot = last_plot(), device = "jpg", dpi = 400, height = 14, width = 14, units = "in")

## PubTator vs LitVar ####
pt_lv_df <- gene_list_1 %>%
  distinct(disease, Gene, .keep_all = TRUE) %>%
  full_join(panel_app_all %>%
              mutate(n_panelapp = 1, disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%
              distinct(disease, Gene, n_panelapp)) %>%
  mutate(query = "First") %>%
  bind_rows(gene_list_2 %>%
              bind_rows(gene_list_1 %>%
                          mutate(context = paste(disease, panelapp_instance, sep = "_")) %>%
                          filter(context %in% complete_query1$keep)) %>%
              distinct(disease, Gene, .keep_all = TRUE) %>%
              full_join(panel_app_all %>%
                          mutate(n_panelapp = 1, disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%
                          distinct(disease, Gene, n_panelapp)) %>%
              mutate(query = "Updated")) %>%
  bind_rows(gene_list_PT %>%
              distinct(disease, Gene, .keep_all = TRUE) %>%
              full_join(panel_app_all %>%
                          mutate(n_panelapp = 1, disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%
                          distinct(disease, Gene, n_panelapp)) %>%
              mutate(query = "Relations")) %>%
  mutate(disease = if_else(disease == "CHH_GMS", "CHH", disease),
         n_panelapp = if_else(is.na(n_panelapp), 0, n_panelapp),
         n_pubtator = if_else(is.na(n), 0, n),
         n_var = if_else(is.na(n_var), 0, n_var),
         n_ortho = if_else(is.na(n_ortho), 0, n_ortho)) %>%
  filter(Gene != "") %>%
  select(query, disease, Gene, n_panelapp, n_pubtator, n_var, n_ortho) %>%
  distinct(query, disease, Gene, .keep_all = TRUE) %>%
  dplyr::rename(PanelApp = n_panelapp, PubTator = n_pubtator, LitVar = n_var, Orthologs = n_ortho) %>%
  mutate(PubTator = if_else(PubTator >= 1, TRUE, FALSE),
         LitVar = if_else(LitVar >= 1, TRUE, FALSE),
         Orthologs = if_else(Orthologs >= 1, TRUE, FALSE),
         PanelApp = if_else(PanelApp >= 1, TRUE, FALSE))

pt_lv_plot <- pt_lv_df %>%
  mutate(query = factor(query, levels = c("First", "Updated", "Relations"))) %>%
  ggplot(mapping = aes(A = PanelApp, B = PubTator, C = LitVar, D = Orthologs)) +
  geom_venn(show_percentage = FALSE, 
            #fill_color = c("#007A82", "#5688A5", "#AC95C7", "#E58EDD"),
            fill_color = c("#007A82", "#112E51", "#FF964F", "#823FFF"),
            #auto_scale = TRUE, 
            stroke_color = "#00363A",
            text_size = 2.5, 
            stroke_alpha = c(0.75, 0, 0, 0), 
            set_name_size = 0, 
            set_name_color = "white",
            digits = 0) + 
  facet_grid(query~disease) +
  theme_void() +
  theme(strip.text.x = element_text(size = 15, vjust = 0), 
        strip.text.y = element_text(size = 15, hjust = 0), 
        panel.spacing = unit(0, "line"))

legend_plot <- ggplot(data.frame(x = 1:4, y = 1, group = factor(c("PanelApp","PubTator","LitVar","Orthologs"),
                                                                levels = c("PanelApp","PubTator","LitVar","Orthologs"))),
                      aes(x, y, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c("PanelApp" = "#007A82", "PubTator" = "#112E51", "LitVar" = "#FF964F", "Orthologs" = "#823FFF")) +
  labs(x = "", y = "") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank())

wrap_plots(pt_lv_plot, legend_plot, widths = c(1, 0))

ggsave(filename = "path/to/Venn_PubTator-LitVar_16-10-2025.jpg", 
       plot = last_plot(), device = "jpg", dpi = 400, height = 6, width = 24, units = "in")



pt_lv_df_instance <- gene_list_1 %>%
  distinct(disease, panelapp_instance, Gene, .keep_all = TRUE) %>%
  full_join(panel_app_all %>%
              mutate(n_panelapp = 1, disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%
              distinct(disease, Gene, n_panelapp)) %>%
  mutate(query = "First") %>%
  bind_rows(gene_list_2 %>%
              bind_rows(gene_list_1 %>%
                          mutate(context = paste(disease, panelapp_instance, sep = "_")) %>%
                          filter(context %in% complete_query1$keep)) %>%
              distinct(disease, panelapp_instance, Gene, .keep_all = TRUE) %>%
              full_join(panel_app_all %>%
                          mutate(n_panelapp = 1, disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%
                          distinct(disease, Gene, n_panelapp)) %>%
              mutate(query = "Updated")) %>%
  bind_rows(gene_list_PT %>%
              distinct(disease, Gene, .keep_all = TRUE) %>%
              full_join(panel_app_all %>%
                          mutate(n_panelapp = 1, disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%
                          distinct(disease, panelapp_instance, Gene, n_panelapp)) %>%
              mutate(query = "Relations")) %>%
  mutate(disease = if_else(disease == "CHH_GMS", "CHH", disease),
         n_panelapp = if_else(is.na(n_panelapp), 0, n_panelapp),
         n_pubtator = if_else(is.na(n), 0, n),
         n_var = if_else(is.na(n_var), 0, n_var),
         n_ortho = if_else(is.na(n_ortho), 0, n_ortho)) %>%
  filter(Gene != "") %>%
  select(query, panelapp_instance, disease, Gene, n_panelapp, n_pubtator, n_var, n_ortho) %>%
  distinct(query, panelapp_instance, disease, Gene, .keep_all = TRUE) %>%
  dplyr::rename(PanelApp = n_panelapp, PubTator = n_pubtator, LitVar = n_var, Orthologs = n_ortho) %>%
  mutate(PubTator = if_else(PubTator >= 1, TRUE, FALSE),
         LitVar = if_else(LitVar >= 1, TRUE, FALSE),
         Orthologs = if_else(Orthologs >= 1, TRUE, FALSE),
         PanelApp = if_else(PanelApp >= 1, TRUE, FALSE))

## Gene ranking ####
gene_list_compare <- gene_list_1 %>%
  mutate(query = "First") %>%
  bind_rows(gene_list_2 %>%
              mutate(query = "Updated")) %>%
  bind_rows(gene_list_PT %>%
              mutate(query = "Relations")) %>%
  mutate(query = factor(query, levels = c("First", "Updated", "Relations")),
         panelapp_instance = if_else(disease == "CHH_GMS", "ENG", panelapp_instance),
         panelapp_instance = if_else(panelapp_instance == "AUS", "Australia", "Genomics England"),
         disease = if_else(disease == "CHH_GMS", "CHH (GMS)", disease))

(rank_plot <- gene_list_compare %>%
  group_by(query, disease, panelapp_instance) %>%
  mutate(rank_total = max(rank, na.rm = TRUE), rank_pct = (rank / rank_total) * 100) %>%
  ungroup() %>%
  filter(color %in% c("green", "amber", "red")) %>%
  
  ggplot(aes(x = rank_pct, fill = disease)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", bins = 25) +
  # scale_fill_manual(guide = guide_legend(override.aes = list(pattern = 'none'))) +
  scale_alpha_manual(values = c(1, 0.5, 0.25)) +
  labs(x = "Normalized Rank", y = "Normalized count", fill = "Panel") +
  facet_grid(query~panelapp_instance) + #, labeller = as_labeller(c(`green` = "Green", `amber` = "Amber", `red` = "Red"))) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12.5)))

ggsave(filename = "path/to/Rank_Queries_16-10-2025.jpg", 
       plot = last_plot(), device = "jpg", dpi = 400, height = 6, width = 10, units = "in")

## GSEA ####
gmt_file_gobp <- "path/to/gene_sets/c5.go.bp.v2023.2.Hs.symbols.gmt"
gmt_file_gomf <- "path/to/gene_sets/c5.go.mf.v2023.2.Hs.symbols.gmt"

gmt_pathways_gobp <- gmtPathways(gmt_file_gobp)
gmt_pathways_gomf <- gmtPathways(gmt_file_gomf)

length(unique(unlist(gmt_pathways_gobp)))
length(unique(unlist(gmt_pathways_gomf)))

rank_and_name <- function(x) {
  x <- x[!is.na(x$n_all), ]
  n <- x$n_all
  names(n) <- x$Gene
  return(n)
}

list_name <- function(x) {
  q <- unique(x$query)
  d <- unique(x$disease)
  return(paste(d, q, sep = "_"))
}

run_fgseabp <- function(x, list) {
  message(paste0("Running ", x))
  return(fgsea(pathways = gmt_pathways_gobp, stats = list[[x]], scoreType = "pos", minSize = 10, maxSize = 100))
}

run_fgseamf <- function(x, list) {
  message(paste0("Running ", x))
  return(fgsea(pathways = gmt_pathways_gomf, stats = list[[x]], scoreType = "pos", minSize = 10, maxSize = 100))
}

### Extra genes ####
# combine gene counts from all queries and remove PanelApp genes
gene_list_gsea <- gene_list_1 %>%
  mutate(query = "query1",
         disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%
  group_by(disease, Gene_ID) %>%
  mutate(remove = if_else(all(panelapp == 0), 0, 1)) %>% 
  ungroup() %>%
  filter(remove == 0) %>%
  bind_rows(gene_list_2 %>%
              mutate(query = "query2",
                     disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%              
              group_by(disease, Gene_ID) %>%
              mutate(remove = if_else(all(panelapp == 0), 0, 1)) %>% 
              ungroup() %>%
              filter(remove == 0)) %>%
  bind_rows(gene_list_PT %>%
              mutate(query = "queryPT",
                     disease = if_else(disease == "CHH_GMS", "CHH", disease)) %>%
              group_by(disease, Gene_ID) %>%
              mutate(remove = if_else(all(panelapp == 0), 0, 1)) %>% 
              ungroup() %>%
              filter(remove == 0)) %>%
  left_join(hs_genes_ncbi %>% 
              mutate(Gene = if_else(Synonyms != "", paste(Symbol, Synonyms, sep = ","), Symbol),
                     Gene_names = Gene) %>%
              select(-Synonyms, -Symbol) %>%
              separate_longer_delim(Gene, delim = ",") %>%
              group_by(GeneID) %>%
              mutate(source = c("main", rep("synonym", n() - 1))) %>%
              ungroup()) %>%
  group_by(query, disease, Gene) %>%
  filter(case_when(any(source == "main") ~ source == "main",
                   TRUE ~ source == "synonym")) %>%
  slice(1) %>%
  ungroup()

# remove empty Gene IDs and no-name Genes
dim(gene_list_gsea)
gene_list_gsea <- gene_list_gsea %>%
  filter(!(is.na(Gene_ID) & !is.na(Gene)))
dim(gene_list_gsea)

gene_list_gsea %>% 
  count(query, disease, Gene) %>%
  filter(n > 1)

# define ranked data for GSEA
count_stats <- gene_list_gsea %>%
  mutate(n = if_else(is.na(n), 0, n),
         n_ortho = if_else(is.na(n_ortho), 0, n_ortho),
         n_var = if_else(is.na(n_var), 0, n_var)) %>%
  filter(Gene_ID %in% hs_genes_ncbi$GeneID) %>%
  rowwise() %>%
  mutate(n_all = sum(n, n_ortho, n_var)) %>%
  arrange(query, disease, desc(n_all)) %>%
  group_by(query, disease) %>%
  mutate(rank_gsea = rank(n_all, ties.method = "random"),
         rank_total = max(rank_gsea, na.rm = TRUE),
         rank_pct = rank_gsea / rank_total,
         n_pct = n_all * rank_pct) %>%
  ungroup() %>%
  select(query, disease, Gene_ID, Gene, n_all, n_pct) %>%
  distinct(query, disease, Gene_ID, .keep_all = TRUE) %>%
  group_by(query, disease) %>%
  group_split()

rank_stats <- lapply(count_stats, rank_and_name)
names(rank_stats) <- unlist(lapply(count_stats, list_name))

# GO - BP
gsea_gobp <- lapply(names(rank_stats), run_fgseabp, list = rank_stats)
names(gsea_gobp) <- names(rank_stats)

gsea_gobp_all <- rbindlist(gsea_gobp, idcol = "name")
gsea_gobp_all <- arrange(gsea_gobp_all, name, pval)

fwrite(gsea_gobp_all, sep = "\t", "path/to/GSEA_GOBP_setSize10-100_Queries1-2-3.txt")

# top 10 for each 
gsea_gobp_all_top10 <- gsea_gobp_all %>%
  separate_wider_delim(cols = name, delim = "_", names = c("disease", "query")) %>%
  group_by(query, disease) %>%
  slice(1:10)

fwrite(gsea_gobp_all_top10, sep = "\t", "path/to/GSEA_GOBP_setSize10-100_Queries1-2-3_TOP10.txt")

# top 5  
gsea_gobp_all_top5 <- gsea_gobp_all %>%
  separate_wider_delim(cols = name, delim = "_", names = c("disease", "query")) %>%
  group_by(query, disease) %>%
  slice(1:5)

fwrite(gsea_gobp_all_top5, sep = "\t", "path/to/GSEA_GOBP_setSize10-100_Queries_TOP5.txt")



# GO - MF 
gsea_gomf <- lapply(names(rank_stats), run_fgseamf, list = rank_stats)
names(gsea_gomf) <- names(rank_stats)

gsea_gomf_all <- rbindlist(gsea_gomf, idcol = "name")
gsea_gomf_all <- arrange(gsea_gomf_all, name, pval)

fwrite(gsea_gomf_all, sep = "\t", "path/to/GSEA_GOMF_setSize10-100_Queries1-2-3.txt")

# top 10 for each 
gsea_gomf_all_top10 <- gsea_gomf_all %>%
  separate_wider_delim(cols = name, delim = "_", names = c("disease", "query")) %>%
  group_by(query, disease) %>%
  slice(1:10)

fwrite(gsea_gomf_all_top10, sep = "\t", "path/to/GSEA_GOMF_setSize10-100_Queries1-2-3_TOP10.txt")

### PanelApp genes ####
# combine gene counts from all queries and remove PanelApp genes
gene_list_gsea <- panel_app_all %>%
  left_join(hs_genes_ncbi %>% 
              mutate(Gene = if_else(Synonyms != "", paste(Symbol, Synonyms, sep = ","), Symbol),
                     Gene_names = Gene) %>%
              select(-Synonyms, -Symbol) %>%
              separate_longer_delim(Gene, delim = ",") %>%
              group_by(GeneID) %>%
              mutate(source = c("main", rep("synonym", n() - 1))) %>%
              ungroup()) %>%
  mutate(disease = if_else(disease == "CHH_GMS", "CHH", disease), color = tolower(color)) %>%
  group_by(disease, Gene) %>%
  filter(case_when(any(source == "main") ~ source == "main",
                   TRUE ~ source == "synonym")) %>%
  slice(1) %>%
  ungroup() %>%
  select(disease, Gene, color)

gene_list_gsea %>%
  count(disease)

# define ranked data for GSEA
count_stats <- gene_list_gsea %>%
  mutate(n_all = case_when(color == "green" ~ 3, color == "amber" ~ 2, color == "red" ~ 1)) %>%
  group_by(disease) %>%
  mutate(rank_gsea = rank(n_all, ties.method = "random")) %>%
  ungroup() %>%
  select(disease, Gene, n_all, rank_gsea) %>%
  distinct(disease, Gene, .keep_all = TRUE) %>%
  group_by(disease) %>%
  group_split()

rank_and_name_pa_color <- function(x) {
  x <- x[!is.na(x$n_all), ]
  n <- x$n_all
  names(n) <- x$Gene
  return(n)
}

rank_and_name_pa_random <- function(x) {
  x <- x[!is.na(x$rank_gsea), ]
  n <- x$rank_gsea
  names(n) <- x$Gene
  return(n)
}

rank_stats <- lapply(count_stats, rank_and_name_pa_color)
names(rank_stats) <- unlist(lapply(count_stats, function(x) {unique(x$disease)}))

gsea_gobp <- lapply(names(rank_stats), run_fgseabp, list = rank_stats)
names(gsea_gobp) <- names(rank_stats)

gsea_gobp_all <- rbindlist(gsea_gobp, idcol = "name")
gsea_gobp_all <- arrange(gsea_gobp_all, name, pval)

fwrite(gsea_gobp_all, sep = "\t", "path/to/PANELAPP_Green3Amber2Red1_GSEA_GOBP_setSize10-100.txt")

### PanelApp / Query comparison #### 
# read in data
pa_gsea <- fread("path/to/PANELAPP_Green3Amber2Red1_GSEA_GOBP_setSize10-100.txt", sep = "\t")
dim(pa_gsea)

queries_gsea <- fread("path/to/GSEA_GOBP_setSize10-100_Queries1-2-3.txt", sep = "\t")

# number of pathways per query and disease
pathway_count <- queries_gsea %>%
  separate_wider_delim(name, delim = "_", names = c("disease", "query")) %>%
  count(query, disease) %>%
  pivot_wider(names_from = "query", values_from = "n") %>%
  left_join(pa_gsea %>%
              count(name) %>%
              dplyr::rename(panelapp = n),
            by = join_by("disease" == "name"))
fwrite(pathway_count, sep = "\t", "path/to/pathway_counts.txt")

# number of significant pathways per query and disease
sig_pathway_count <- queries_gsea %>%
  separate_wider_delim(name, delim = "_", names = c("disease", "query")) %>%
  filter(pval < 0.05) %>%
  count(query, disease) %>%
  pivot_wider(names_from = "query", values_from = "n") %>%
  left_join(pa_gsea %>%
              filter(pval < 0.05) %>%
              count(name) %>%
              dplyr::rename(panelapp = n),
            by = join_by("disease" == "name"))
fwrite(sig_pathway_count, sep = "\t", "path/to/sig_pathway_counts.txt")

# proportion significant pathways
prop_sig_pathway_count <- cbind(disease = sig_pathway_count$disease, 
                                (sig_pathway_count[,-1] / pathway_count[,-1])*100)
fwrite(prop_sig_pathway_count, sep = "\t", "path/to/prop_sig_pathway_counts.txt")

queries_gsea_rank <- queries_gsea %>%
  separate_wider_delim(name, delim = "_", names = c("disease", "query")) %>%
  group_by(disease, query) %>%
  mutate(rank = dense_rank(pval)) %>%
  select(query, disease, pathway, rank) %>%
  pivot_wider(names_from = "query", values_from = "rank")

queries_gsea_pval <- queries_gsea %>%
  separate_wider_delim(name, delim = "_", names = c("disease", "query")) %>%
  select(query, disease, pathway, pval) %>%
  pivot_wider(names_from = "query", values_from = "pval")

combined_gsea_rank <- queries_gsea_rank %>%
  full_join(pa_gsea %>%
              group_by(name) %>%
              mutate(panelapp = dense_rank(pval)) %>%
              select(name, pathway, panelapp), 
            by = join_by(disease == name, pathway == pathway))

combined_gsea_pval <- queries_gsea_pval %>%
  full_join(pa_gsea %>%
              dplyr::rename(panelapp = pval) %>%
              select(name, pathway, panelapp), 
            by = join_by(disease == name, pathway == pathway))

combined_gsea_filter <- combined_gsea_pval %>%
  # filter(!is.na(query1) & !is.na(query2) & !is.na(queryPT) & !is.na(panelapp)) %>%
  left_join(combined_gsea_rank %>%
              dplyr::rename(rank_query1 = query1, rank_query2 = query2,
                            rank_queryPT = queryPT, rank_panelapp = panelapp), 
            by = c("disease", "pathway")) %>%
  select(disease, pathway, panelapp, rank_panelapp, query1, rank_query1, query2, rank_query2, queryPT, rank_queryPT) %>%
  filter(panelapp < 0.05) %>%
  arrange(disease, panelapp)
fwrite(combined_gsea_filter, sep = "\t", "path/to/gsea_panelapp_queries.txt")

combined_gsea_filter_overlap_panelapp <- queries_gsea %>%
  separate_wider_delim(name, delim = "_", names = c("disease", "query")) %>%
  select(query, disease, pathway, pval) %>%
  filter(pval < 0.05) %>%
  full_join(pa_gsea %>%
              select(name, pathway, pval) %>%
              filter(pval < 0.05) %>%
              dplyr::rename(panelapp = pval),
            by = join_by(disease == name, pathway == pathway)) %>%
  filter(panelapp < 0.05 & pval < 0.05) %>%
  count(query, disease) %>%
  pivot_wider(names_from = "query", values_from = "n")
fwrite(combined_gsea_filter_overlap_panelapp, sep = "\t", "path/to/sig_overlap_panelapp.txt")

dim(combined_gsea_filter)
unique(combined_gsea_filter$disease)
unique(pa_gsea$name)
