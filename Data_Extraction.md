Start by downloading PanelApp data

```R
# DATA EXTRACTION ####
source(PUBPANEL_FUNCTIONS.R)

# PanelApp Genomics England
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

fwrite(panelapp_genes_eng_df_tabreplace, sep = "\t", file = "path/to/PanelAppEngland_panels_DATE.txt")

# PanelApp Australia
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
fwrite(panelapp_genes_aus_df_tabreplace, sep = "\t", file = "path/to/PanelAppAustralia_panels_DATE.txt"))

# Combined
panel_app_all <- panelapp_genes_eng_df_tabreplace %>%
  mutate(panelapp = "ENG") %>%
  bind_rows(panelapp_genes_aus_df_tabreplace %>%
              mutate(panelapp = "AUS"))
dim(panel_app_all)

fwrite(panel_app_all, sep = "\t", file = "path/to/PanelAppENG-AUS_panels_DATE.txt")
```
