# DEFINE FUNCTIONS ####
## HELPER ####

chunk <- function(x, n) {split(seq_along(1:x), ceiling(seq_along(1:x)/n))}

MeSH_to_PubTator <- function(x) {  return(paste0("@DISEASE_", gsub(",\\s|-|\\s", "_", x)))}

paste_unique <- function(x, c = ";") {return(paste(unique(x), collapse = c))}

pubtator_df_clean <- function(lines) {
  lines <- str_replace(lines, "(;\\d{8})$", "\\1;")
  is_continuation <- str_detect(lines, "\\\t$")
  result <- {
    group <- integer(length(lines))
    g <- 0
    
    for (i in seq_along(lines)) {
      if (!is_continuation[i]) {
        g <- g + 1
      }
      group[i] <- g
    }
    
    tapply(lines, group, function(grp) {
      str_c(str_remove_all(grp, "\\\t$"), collapse = "")
    }) %>% 
      unname()
  }
  
  result_df <- read.delim(textConnection(result), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  return(result_df)
}

tidy_variant <- function(x) {
  
  require(stringr)
  
  vars <- paste_unique(x, c = "|")
  
  vars <- unlist(str_split(x, "\\|"))
  vars_edit <- c()
  for (var in 1:length(vars)) {
    if (!is.null(vars[var])) {
      if (grepl("^rs[0-9]", vars[var])) { # RS number
        vars_edit <- c(vars_edit, str_extract(vars[var], "rs[0-9]+"))
      } else if (grepl("#[0-9]+#[cgp].+", vars[var])) { # gene then HGVS
        vars_edit <- c(vars_edit, str_extract(vars[var], "(?<=#)[0-9]+#.+(?=~)"))
      } else if (grepl("^##[cgp].+", vars[var])) { # HGVS
        vars_edit <- c(vars_edit, str_extract(vars[var], "(?<=##)[cgp].+(?=~)"))
      } else if (grepl("^#[0-9]+#", vars[var])) { # gene
        vars_edit <- c(vars_edit, str_extract(vars[var], "#[0-9]+#(?=~)"))
      } else if (grepl("NA", vars[var])) { # NA
        vars_edit <- c(vars_edit, NA)
      } else {
        vars_edit <- c(vars_edit, vars[var])
      }
    } else {
      vars_edit <- c(vars_edit, NULL)
    }
  }
  vars_edit <- vars_edit[vars_edit != ""]
  return(paste_unique(vars_edit))
}

## API ####
### PanelApp ####
PanelIDs_PanelApp <- function(x) {
  lapply(c("httr", "jsonlite"), require, character.only = TRUE)
  
  if (!x %in% c("ENG", "AUS")) {
    stop("Invalid input: must be either 'ENG' or 'AUS'")
  }
  
  if (x == "ENG") {
    get <- RETRY("GET", times = 10, "https://panelapp.genomicsengland.co.uk/api/v1/panels/?format=json")
  } else {
    get <- RETRY("GET", times = 10, "https://panelapp-aus.org/api/v1/panels/?format=json")
  }
  data <- fromJSON(rawToChar(get$content))
  panel_ids <- data$results$id
  
  while (!is.null(data$`next`)) {
    get <- RETRY("GET", data$`next`, times = 10)
    data <- fromJSON(rawToChar(get$content))
    panel_ids <- c(panel_ids, data$results$id)
  }
  return(panel_ids)
}

PanelGenes_PanelApp <- function(x, panelapp = c("ENG", "AUS")) {
  lapply(c("httr", "jsonlite"), require, character.only = TRUE)
  
  if (!panelapp %in% c("ENG", "AUS")) {
    stop("Invalid input: must be either 'ENG' or 'AUS'")
  }
  
  if (panelapp == "ENG") {
    base_url <- "https://panelapp.genomicsengland.co.uk/api/v1/panels/"
  } else {
    base_url <- "https://panelapp-aus.org/api/v1/panels/"
  }
  
  get <- RETRY(verb = "GET", times = 10, pause_cap = 1, quiet = FALSE, terminate_on = c(400),
               url = paste0(base_url, x, "/?format=json"))
  
  data <- fromJSON(rawToChar(get$content))

  if (length(data$genes) > 0) {
    df <- data$genes[,2:ncol(data$genes)]

    df$publications <- sapply(df$publications, paste, collapse = "~")
    df$evidence <- sapply(df$evidence, paste, collapse = "~")
    df$phenotypes <- sapply(df$phenotypes, paste, collapse = "~")
    df$transcript <- sapply(df$transcript, paste, collapse = "~")
    df$tags <- sapply(df$tags, paste, collapse = "~")
    
    df$panel <- data$name
    df$id <- data$id
    df$disease_group <- data$disease_group
    df$disease_sub_group <- data$disease_sub_group
    df$version <- data$version 
  } else {
    df <- data.frame(panel = data$name, id = data$id, disease_group = data$disease_group,
                     disease_sub_group = data$disease_sub_group, version = data$version)
  }
  return(df)
}

### PubTator ####
pubtator_function_JSON_UPDATE <- function (x) {
  lapply(c("httr", "jsonlite"), require, character.only = TRUE)
  
  get <- RETRY(verb = "GET", times = 10, pause_cap = 1, quiet = FALSE, terminate_on = c(400),
               url = paste("https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson?pmids=", x, sep = ""))

  data <- fromJSON(rawToChar(get$content))
  
  gene <- disease <- chemical <- variant <- species <- cellline <- NULL

  temp_annot <- data$PubTator3$passages[[1]]$annotations

  if (length(temp_annot) != 0) {
    for (i in 1:length(temp_annot)) {
      temp_annot_infon <- temp_annot[[i]]$infons
      if (!is.null(temp_annot_infon)) {
        for (j in 1:nrow(temp_annot_infon)) {
          if (temp_annot_infon[j,]$type == "Gene") {
            gene <- c(gene, paste(temp_annot_infon[j,]$normalized_id, temp_annot_infon[j,]$name, sep = "~"))
          } else if (temp_annot_infon[j,]$type == "Disease") { 
            disease <- c(disease, paste(temp_annot_infon[j,]$normalized_id, temp_annot_infon[j,]$name, sep = "~"))
          } else if (temp_annot_infon[j,]$type == "Chemical") {
            chemical <- c(chemical, paste(temp_annot_infon[j,]$normalized_id, temp_annot_infon[j,]$name, sep = "~"))
          } else if (temp_annot_infon[j,]$type == "Variant") {
            variant <- c(variant, paste(temp_annot_infon[j,]$normalized_id, temp_annot_infon[j,]$name, sep = "~"))
          } else if (temp_annot_infon[j,]$type == "Species") {
            species <- c(species, paste(temp_annot_infon[j,]$normalized_id, temp_annot_infon[j,]$name, sep = "~"))
          } else if (temp_annot_infon[j,]$type == "CellLine")  {
            cellline <- c(cellline, paste(temp_annot_infon[j,]$normalized_id, temp_annot_infon[j,]$name, sep = "~"))
          }
        }
      }
    }
    gene <- union(gene, gene)
    disease <- union(disease, disease)
    variant <- union(variant, variant)
    chemical <- union(chemical, chemical)
    species <- union(species, species)
    cellline <- union(cellline, cellline)
  }

  temp_relat <- data$PubTator3$relations_display[[1]]
  if(length(temp_relat) > 0) {
    relat <- paste(gsub(";", "", unlist(temp_relat)), collapse = "~")
  } else {
    relat <- NULL
  }
      
  return(list(PMID = x, Genes = gene, Diseases = disease, Chemicals = chemical, 
              Variants = variant, Species = species, CellLines = cellline, 
              Relations = relat))
}

pubtator_to_table_UPDATE <- function (x, file = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/PT_Master.txt", write = TRUE) {

  # if wanting to write to file
  if (write) {
    if (!endsWith(file, ".txt")) {
      message("Please make file '.txt'")
    }
    
    if (!file.exists(file)) {
      file.create(file)
      write(c("PMID", "Genes", "Diseases", "Chemicals", "Variants", "Species", "CellLines", "Relations"), 
            file = file, append = TRUE, sep = "\t", ncolumns = 8)
    }
    
    # if writing, will write line-by-line: lapply(PMIDs, function(x) _table(_PubTator))
    for (i in 1:length(x)) {
      write.table(cbind(x[[i]]$PMID,
                        paste(x[[i]]$Genes, collapse = "|"),
                        paste(x[[i]]$Diseases, collapse = "|"),
                        paste(x[[i]]$Chemicals, collapse = "|"),
                        paste(x[[i]]$Variants, collapse = "|"),
                        paste(x[[i]]$Species, collapse = "|"),
                        paste(x[[i]]$CellLines, collapse = "|"),
                        paste(x[[i]]$Relations, collapse = "|")), 
                  file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  } else { # if NOT writing, will collapse all in one go: _table(lapply(PMIDS, _PubTator))
    df <- data.frame(PMID = numeric(), Genes = character(), Diseases = character(), Chemicals = character(),
                     Variants = character(), Species = character(), CellLines = character(), Relations = character())
    for (i in 1:length(x)) {
      df <- rbind(df, cbind(x[[i]]$PMID, 
                            paste(x[[i]]$Genes, collapse = "|"), paste(x[[i]]$Diseases, collapse = "|"),
                            paste(x[[i]]$Chemicals, collapse = "|"), paste(x[[i]]$Variants, collapse = "|"),
                            paste(x[[i]]$Species, collapse = "|"), paste(x[[i]]$CellLines, collapse = "|"),
                            paste(x[[i]]$Relations, collapse = ",")))
    }
    colnames(df) <- c("PMID", "Genes", "Diseases", "Chemicals", "Variants", "Species", "CellLines", "Relations")
    return(df)
  }
}

PMIDs_pubtator_query <- function (q, date = "today", write = TRUE, type = "genetic", force = FALSE,
                                  master = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/RE_Master.txt") {
  lapply(c("httr", "jsonlite"), require, character.only = TRUE)
  
  # if wanting to write to file
  if (write) {
    if (!endsWith(master, ".txt")) {
      stop("Please make file '.txt'")
    }
  }
  
  if (!file.exists(master)) {
    file.create(master)
    write(c("date", "query", "pmids"), 
          file = master, append = TRUE, sep = "\t", ncolumns = 3)
  }
  
  if (!type %in% c("genetic", "disease")) {stop("Argument 'type' must be 'genetic' or 'disease'")}
  
  if (date == "today") {date <- gsub("-", "", as.character(Sys.Date()))} 
  
  pubtator_query_base <- "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?format=json&text="
  
  if (type == "genetic") {
    gene_q <- unlist(lapply(q, function(x) paste("relations:ANY|", x, "|GENE", sep = "")))
    variant_q <- unlist(lapply(q, function(x) paste("relations:ANY|", x, "|VARIANT", sep = "")))
    
    queries <- c(gene_q, variant_q)

  } else if (type == "disease") {
    disease_q <- unlist(lapply(q, function(x) paste("relations:ANY|", x, "|DISEASE", sep = "")))
    
    queries <- c(disease_q)
  }

  if (!force) {
    queries_date <- paste(date, queries, sep = "~")
    queries_pass <- find_in_master(queries_date, file = master, master_type = "RE")

    if (length(queries_pass) < 1) {
      stop("All queries already in master file")
    }
  } else {
    queries_pass <- queries
  }

  queries_full <- paste0(pubtator_query_base, queries_pass)
  
  pmid_df <- data.frame(matrix(ncol = 3, nrow = length(queries_full)))
  names(pmid_df) <- c("date", "query", "pmids")
  
  for (q in 1:length(queries_full)) {
    message("Processing query ", q, " of ", length(queries_full))
    get <- RETRY(verb = "GET", url = queries_full[q], times = 10, pause_cap = 1, quiet = FALSE, terminate_on = c(400))

    pmid <- NULL
    data <- fromJSON(rawToChar(get$content))
    test_pages <- data$total_pages
    
    if (test_pages == 0) {
      pmid <- "No publications"
    } else {
      message(paste0("Extracting from ", test_pages, " pages"))
      for (page in 1:test_pages) {
        test_page <- fromJSON(paste(queries_full[q], "&page=", page, sep = ""))
        pmid <- c(pmid, test_page$results$pmid)
        Sys.sleep(0.25)
      }
    }
    pmid_df[q,]$date <- date
    pmid_df[q,]$query <- queries_pass[q]
    pmid_df[q,]$pmids <- paste(pmid, collapse = ";")
    
    if (write) {
      write.table(cbind(date, queries_pass[q], paste(pmid, collapse = ";")), 
                  file = master,  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
  return(pmid_df)
}

### LitVar ####
litvaR <- function(x, fn = c("s", "q", "p")) {
  lapply(c("httr", "jsonlite", "urltools"), require, character.only = TRUE)
  
  fn <- match.arg(fn)
  
  if (fn == "s") {
    if (grepl("^litvar@", x)) {
      x <- gsub("^litvar@", "", x)
    }
    
    if (!grepl("#", x)) {
      message(paste0(x, " must be VarID (e.g., rs## or #GeneID#HGVS)"))
      return(data.frame(id = x))
    }
    
    x_encode <- url_encode(x)

    get <- RETRY(verb = "GET", times = 10, pause_cap = 1, quiet = FALSE, terminate_on = c(400),
                 url = paste0("https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/litvar@", x_encode))

    if (get$status_code == 200) {
      data <- fromJSON(rawToChar(get$content))
      data <- as.data.frame(t(unlist(data)))
      colnames(data)[colnames(data) == "_id"] <- "id"
      data$input <- x
      return(data)
    } else {
      return(data.frame(id = x, input = x))
    }
  } else if (fn == "q") {
    x_encode <- url_encode(x)
    
    get <- RETRY(verb = "GET", times = 10, pause_cap = 1, quiet = FALSE, terminate_on = c(400),
                 url = paste0("https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/?query=", x_encode))
    
    if (get$status_code == 200) {
      data <- fromJSON(rawToChar(get$content))
      return(data)
    } else {
      stop("Error in 'x' or accession")
    }
  } else if (fn == "p") {

    if (grepl("^litvar@", x)) {
      x <- gsub("^litvar@", "", x)
    }
    
    x_encode <- url_encode(x)
    
    get <- RETRY("GET", times = 10, pause_cap = 1, quiet = FALSE, terminate_on = c(400),
                 paste0("https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/litvar@", x_encode, "/publications"))
    
    if (get$status_code == 200) {
      data <- fromJSON(rawToChar(get$content))
      if ("pmids" %in% names(data)) {
        return(data$pmids)
      } else {
        return(NA)
      }
    } else {
      return <- fromJSON(rawToChar(get$content))
      if (grepl("Variant not found", return$detail)) {
        return(NA)
      }
    }
  }
}

litvaR_to_table <- function (x, file = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/LV_Master.txt") {

  if (!endsWith(file, ".txt")) {
    message("Please make file '.txt'")
  }
  
  if (!file.exists(file)) {
    file.create(file)
    write(c("id", "id_pubtator", "genes", "species", "clinical"), 
          file = file, append = TRUE, sep = "\t", ncolumns = 8)
  }
  
  for (i in 1:length(x)) {
    a <- gsub("litvar@", "", x[[i]]$id)
    b <- x[[i]]$input
    write.table(cbind(x[[i]]$id,
                      ifelse(a == b, a, b),
                      ifelse(length(x[[i]][,grepl("^gene", colnames(x[[i]]))]) > 0, paste(x[[i]][,grepl("^gene", colnames(x[[i]]))], collapse = "|"), ""),
                      ifelse(length(x[[i]][,grepl("^data_species", colnames(x[[i]]))]) > 0, paste(x[[i]][,grepl("^data_species", colnames(x[[i]]))], x[[i]][,grepl("^data_tax_id", colnames(x[[i]]))], sep = "~"), ""),
                      ifelse(length(x[[i]][,grepl("^data_clinical_significance", colnames(x[[i]]))]) > 0, paste(x[[i]][,grepl("^data_clinical_significance", colnames(x[[i]]))], collapse = "|"), "")
                      ),
                file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

## ANALYSIS ####
find_in_master <- function(x, file = NULL, master_type = "PT") {
  
  if (master_type == "PT") {
    if (is.null(file)) {
      master_file <- "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/PT_Master.txt"
    } else {
      master_file = file
    }
    master_pt <- fread(master_file, sep = "\t", fill = TRUE)
    x <- x[!x %in% master_pt$PMID]
    rm(master_pt)
    return(x)
    
  } else if (master_type == "LV") {
    if (is.null(file)) {
      master_file <- "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/LV_Master.txt"
    } else {
      master_file = file
    }
    master_lv <- fread(master_file, sep = "\t", fill = TRUE)
    x <- x[!x %in% master_lv$id_pubtator]
    rm(master_lv)
    return(x)
    
  } else if (master_type == "RE") {
    if (is.null(file)) {
      master_file <- "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/RE_Master.txt"
    } else {
      master_file = file
    }
    master_re <- readLines(master_file)
    master_re <- pubtator_df_clean(master_re)
    master_re <- master_re %>% 
      separate_longer_delim(cols = pmids, delim = ";") %>% 
      group_by(date, query) %>% 
      reframe(pmids = paste_unique(pmids)) %>%
      mutate(date_query = paste(date, query, sep = "~"))
    x <- x[!x %in% master_re$date_query]
    x <- gsub(pattern = ".*~", "", x)
    rm(master_re)
    return(x)
    
  } else {
    stop("master_type must be either 'PT' (PubTator), 'LV' (LitVar), or 'RE' (PubTator Relations)")
  }
}

annotatoR <- function(ids, chunk_size_pt = 50, chunk_start_pt = 1, chunk_size_lv = 50, chunk_start_lv = 1,
                      pubtator = TRUE, litvar = TRUE, 
                      pt_master = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/PT_Master.txt") {
  lapply(c("data.table", "dplyr", "tidyr"), require, character.only = TRUE)
  
  if (pubtator) {
    pmids <- ids
    
    if (is.numeric(pmids)) {
      pmids <- as.character(pmids)
    }
    
    message(paste0("Checking ", length(pmids), " PMIDs in Master"))
    pmids_annot <- find_in_master(pmids, master_type = "PT")

    if (length(pmids_annot) > 0) {
      message("Retrieving PubTator annotations for ", length(pmids_annot), " PMIDs")
      div_pmids <- chunk(length(pmids_annot), chunk_size_pt)
      for (block in chunk_start_pt:length(div_pmids)) {
        message(paste0("Working on ", block, "/", length(div_pmids)))
        pubtator_to_table_UPDATE(lapply(pmids_annot[div_pmids[[block]]], pubtator_function_JSON_UPDATE))
      }
    }
  }

  if (litvar) {
    if (pubtator) {
      df <- fread(pt_master, sep = "\t")
      vars <- df %>%
        separate_longer_delim(cols = "Variants", delim = "|") %>%
        filter(Variants != "" & Variants != "NA~NA" & Variants != "Could not retrieve publications") %>%
        separate_wider_delim(cols = "Variants", delim = "~", names = c("varID", "var")) %>%
        filter(grepl("#", varID)) %>%
        pull(varID) %>%
        unique()
    } else {
      vars <- ids
    }
    
    if (!all(grepl("#", vars))) {
      stop(paste0("'x' contains: ", vars[!grepl("#", vars)]), "; must contain only VarIDs (e.g., rs## or #GeneID#HGVS)")
    }
    vars_annot <- find_in_master(vars, master_type = "LV")

    if (length(vars_annot) > 0) {
      message(paste0("Retrieving LitVar annotations for ", length(vars_annot), " variants"))
      div_vars <- chunk(length(vars_annot), chunk_size_lv)
      for (block in chunk_start_lv:length(div_vars)) {
        message(paste0("Working on ", block, "/", length(div_vars)))
        litvaR_to_table(lapply(vars_annot[div_vars[[block]]], litvaR, fn = "s"))
      }
    }
  }
}

gene_countR <- function(pubtator, ortholog = TRUE, variant = TRUE, pa = NULL, pa_id = NULL,
                        hs_genes = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/ncbi_dataset_SYN.tsv", 
                        ortholog_db = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/gene_orthologs_MOD.txt",
                        variant_db = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/LV_Master.txt") {
  # load dependencies
  lapply(c("data.table", "dplyr", "tidyr"), require, character.only = TRUE)
  
  hs_genes_ncbi <- fread(hs_genes, sep = "\t")
  colnames(hs_genes_ncbi)[1] <- gsub("NCBI ", "", colnames(hs_genes_ncbi)[1])
  
  if (is.character(pubtator)) {
    df <- fread(pubtator, sep = "\t")
  } else {
    df <- pubtator
  }
  
  message("computing gene counts")
  df_gene <- df %>%
    separate_longer_delim(cols = "Genes", delim = "|") %>% 
    filter(Genes != "" & Genes != "Could not retrieve publications") %>%
    separate_wider_delim(cols = "Genes", delim = "~", names = c("Gene_ID", "Gene_name")) %>%
    mutate(Gene_ID = suppressWarnings(as.numeric(Gene_ID))) %>%
    filter(Gene_ID %in% hs_genes_ncbi$GeneID) %>%
    select(Gene_ID, Gene_name, PMID, Variants) %>%
    distinct(Gene_ID, Gene_name, PMID, .keep_all = TRUE) %>%
    group_by(Gene_ID, Gene_name) %>%
    summarize(PMID = paste_unique(PMID), n = n(), Variants = tidy_variant(Variants), .groups = "drop") %>%
    distinct() %>%
    arrange(desc(n)) %>%
    mutate(rank = dense_rank(desc(n))) %>%
    relocate(PMID, .after = last_col())
  
  # gene ortholog counts
  if (ortholog) {
    if (!exists("ortholog_df")) {
      ortho <- fread(ortholog_db)
    }
    
    message("computing orthologous gene counts")
    df_ortholog <- df %>%
      separate_longer_delim(cols = "Genes", delim = "|") %>% 
      filter(Genes != "" & Genes != "Could not retrieve publications") %>%
      separate_wider_delim(cols = "Genes", delim = "~", names = c("Gene_ID", "Gene_name")) %>%
      mutate(Gene_ID = suppressWarnings(as.numeric(Gene_ID))) %>%
      filter(!Gene_ID %in% hs_genes_ncbi$GeneID) %>%
      select(Gene_ID, Gene_name, PMID) %>%
      distinct() %>%
      left_join(ortho, by = join_by(Gene_ID == ortholog_gene)) %>%
      left_join(hs_genes_ncbi, by = join_by(hs_gene == GeneID)) %>%
      mutate(Gene_join = if_else(!is.na(hs_gene), hs_gene, Gene_ID),
             Symbol_join = if_else(!is.na(Symbol), Symbol, Gene_name)) %>%
      group_by(Gene_join, Symbol_join) %>%
      summarize(n_ortho = n(), genes_ortho = paste_unique(paste0(Gene_ID, "=", Gene_name, "(", ortholog_tax, ")")),
                PMID_ortho = paste_unique(PMID), .groups = "drop") %>%
      select(Gene_join, Symbol_join, n_ortho, genes_ortho, PMID_ortho) %>%
      rename(Gene_ID = Gene_join, Gene_name = Symbol_join) %>%
      distinct() %>%
      filter(!is.na(Gene_name)) %>%
      mutate(rank_ortho = dense_rank(desc(n_ortho)))
    
    df_gene <- df_gene %>%
      full_join(df_ortholog, by = c("Gene_ID", "Gene_name"))
  }
  
  # variant gene counts
  if (variant) {
    if (!exists("variant_df")) {
      variant_df <- fread(variant_db, sep = "\t", fill = TRUE)
    }
    
    message("computing variant gene counts")
    df_vars <- df %>%
      separate_longer_delim(cols = "Variants", delim = "|") %>%
      filter(Variants != "" & Variants != "NA~NA" & Variants != "Could not retrieve publications") %>%
      separate_wider_delim(cols = "Variants", delim = "~", names = c("varID", "var")) %>%
      left_join(variant_df, by = join_by(varID == id_pubtator)) %>%
      separate_longer_delim(cols = "genes", delim = "|") %>%
      filter(!is.na(genes) & genes != "") %>%
      select(PMID, id, genes) %>%
      distinct() %>%
      group_by(genes) %>%
      filter(!is.na(id)) %>%
      summarize(PMID_var = paste_unique(PMID), id_var = paste_unique(id), n_var = n(), .groups = "drop") %>%
      distinct() %>%
      arrange(desc(n_var)) %>%
      mutate(rank_var = dense_rank(desc(n_var)))
    
    df_gene <- df_gene %>%
      full_join(df_vars, by = join_by(Gene_name == genes))
  }
  
  # relations
  df_gene_relations <- df %>%
    separate_longer_delim(cols = "Genes", delim = "|") %>% 
    filter(Genes != "" & Genes != "Could not retrieve publications") %>%
    separate_wider_delim(cols = "Genes", delim = "~", names = c("Gene_ID", "Gene_name")) %>%
    separate_longer_delim(cols = "Relations", delim = "~") %>%
    filter(Relations != "" & Relations != "Could not retrieve publications") %>%
    separate_wider_delim(cols = "Relations", delim = "|", names = c("Relation", "Entity1", "Entity2"), cols_remove = FALSE) %>%
    select(PMID, Relations) %>%
    distinct() %>%
    filter(!is.na(Relations)) %>%
    group_by(Relations) %>%
    summarize(n = n(), PMID_relations = paste_unique(PMID), .groups = "drop") %>%
    arrange(desc(n)) %>%
    separate_wider_delim(cols = Relations, delim = "|", names = c("Relation", "Entity1", "Entity2"))

  # panel app comparison
  if (is.null(pa)) {
    return(list(gene_count = df_gene, relations_count = df_gene_relations))
  } else {
    panel <- subset(pa, pa$id == pa_id)
    panel_genes <- panel$entity_name
    
    panel_genes_green <- subset(panel, color == "Green")$entity_name
    panel_genes_amber <- subset(panel, color == "Amber")$entity_name
    panel_genes_red <- subset(panel, color == "Red")$entity_name
    
    df_gene_syn <- df_gene %>%
      mutate(Gene_ID = as.numeric(Gene_ID)) %>% 
      left_join(hs_genes_ncbi, by = join_by(Gene_ID == GeneID)) %>% 
      mutate(Gene_names = if_else(Gene_name == Symbol, 
                                  if_else(Synonyms != "", paste(Gene_name, Synonyms, sep = ","), Gene_name),
                                  if_else(Synonyms != "", paste(Symbol, Gene_name, Synonyms, sep = ","), paste(Symbol, Gene_name, sep = ","))),
             Gene = Gene_names) %>%
      separate_longer_delim(cols = Gene, delim = ",") %>%
      select(-c(Symbol, Synonyms)) %>%
      group_by(Gene_ID) %>%
      mutate(source = c("main", rep("synonym", n() - 1))) %>%
      ungroup()

    present <- panel_genes[panel_genes %in% unique(df_gene_syn$Gene)]
    present_green <- panel_genes_green[panel_genes_green %in% unique(df_gene_syn$Gene)]
    present_amber <- panel_genes_amber[panel_genes_amber %in% unique(df_gene_syn$Gene)]
    present_red <- panel_genes_red[panel_genes_red %in% unique(df_gene_syn$Gene)]
    
    df_gene <- df_gene_syn %>%
      left_join(panel[,c("entity_name", "color")], by = join_by(Gene == entity_name)) %>%
      group_by(Gene_ID) %>%
      filter(case_when(any(source == "main" & !is.na(color)) ~ source == "main",
                       any(source == "synonym" & !is.na(color)) ~ source == "synonym" & !is.na(color),
                       TRUE ~ source == "main")) %>%
      ungroup() %>%
      mutate(panelapp = if_else(!is.na(color), 1, 0), color = tolower(color)) %>%
      select(Gene_ID, Gene, Gene_names, panelapp, color, Variants,
             n, rank, PMID, n_ortho, rank_ortho, genes_ortho, PMID_ortho, n_var, rank_var, id_var, PMID_var)

    missing <- panel_genes[!panel_genes %in% unique(df_gene$Gene)]
    missing_green <- panel_genes_green[!panel_genes_green %in% unique(df_gene$Gene)]
    missing_amber <- panel_genes_amber[!panel_genes_amber %in% unique(df_gene$Gene)]
    missing_red <- panel_genes_red[!panel_genes_red %in% unique(df_gene$Gene)]
    
    extra <- df_gene[!df_gene$Gene %in% panel_genes, ]$Gene
    
    missing_genes <- paste(sort(missing), collapse = ";")
    missing_genes_green <- paste(sort(missing_green), collapse = ";")
    missing_genes_amber <- paste(sort(missing_amber), collapse = ";")
    missing_genes_red <- paste(sort(missing_red), collapse = ";")
    
    pct_present <- length(present)/length(panel_genes)
    pct_present_green <- length(present_green)/length(panel_genes_green)
    pct_present_amber <- length(present_amber)/length(panel_genes_amber)
    pct_present_red <- length(present_red)/length(panel_genes_red)
    
    precision <- length(present)/nrow(df_gene)
    
    df_stats <- data.frame(stat = c("n_panel", "n_genes", "n_present", "tpr", "tpr_green", 
                                    "tpr_amber", "tpr_red", "n_extra", "precision", "n_articles"),
                           value = c(length(panel_genes), nrow(df_gene), length(present), 
                                     pct_present, pct_present_green, pct_present_amber, 
                                     pct_present_red, length(extra), precision, length(df$PMID)))
    
    return(list(panel = unique(panel$panel), panelapp_genes = sort(panel_genes), gene_count = df_gene, 
                relations_count = df_gene_relations, stats = df_stats, missing = paste(missing, collapse = ";"),
                missing_green = paste(missing_green, collapse = ";"), missing_amber = paste(missing_amber, collapse = ";"),
                missing_red = paste(missing_red, collapse = ";")))
  }
}
