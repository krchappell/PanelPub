
# PanelPub functions ####

annotatoR <- function(ids, chunk_size_pt = 50, chunk_start_pt = 1, chunk_size_lv = 50, chunk_start_lv = 1,
                      pubtator = TRUE, litvar = TRUE, 
                      pt_master = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/PT_Master.txt") {
  require(dplyr)
  
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

build_plus_venn <- function(params, D=15, scale_radius=1){
  center_left   <- c(-D, 0)
  center_right  <- c( D, 0)
  center_top    <- c(0,  D)
  center_bottom <- c(0, -D)
  
  df_left <- build_venn_df(params$left$N1, params$left$N2, params$left$Noverlap,
                           "left", "left", scale_radius, center_left)
  df_right <- build_venn_df(params$right$N1, params$right$N2, params$right$Noverlap,
                            "right", "right", scale_radius, center_right)
  df_top <- build_venn_df(params$top$N1, params$top$N2, params$top$Noverlap,
                          "top", "top", scale_radius, center_top)
  df_bottom <- build_venn_df(params$bottom$N1, params$bottom$N2, params$bottom$Noverlap,
                             "bottom", "bottom", scale_radius, center_bottom)
  
  labels_left <- compute_labels(df_left, "left", params$left$N1, params$left$N2,
                                params$left$Noverlap, "left")
  labels_right <- compute_labels(df_right, "right", params$right$N1, params$right$N2,
                                 params$right$Noverlap, "right")
  labels_top <- compute_labels(df_top, "top", params$top$N1, params$top$N2,
                               params$top$Noverlap, "top")
  labels_bottom <- compute_labels(df_bottom, "bottom", params$bottom$N1,
                                  params$bottom$N2, params$bottom$Noverlap, "bottom")
  
  list(
    polygons = as.data.frame(rbind(df_left, df_right, df_top, df_bottom)),
    labels = as.data.frame(rbind(labels_left, labels_right, labels_top, labels_bottom))
  )
}


build_venn_df <- function(N1, N2, Noverlap, id, layout, scale_radius,
                          inner_center) {
  
  R1 <- sqrt(N1/pi) * scale_radius
  R2 <- sqrt(N2/pi) * scale_radius
  
  # Identify inner vs outer
  if (R1 < R2) {
    inner_set <- "A"; outer_set <- "B"; R_inner <- R1; R_outer <- R2
  } else if (R2 < R1) {
    inner_set <- "B"; outer_set <- "A"; R_inner <- R2; R_outer <- R1
  } else {
    inner_set <- "A"; outer_set <- "B"; R_inner <- R1; R_outer <- R2
  }
  
  # Compute distance for desired overlap
  center_distance <- if(Noverlap >= min(N1,N2)) abs(R_outer - R_inner) else
    uniroot(function(d) intersection_area(d, R1, R2) - Noverlap,
            interval = c(0,R1+R2))$root
  
  dir_vec <- switch(layout, "left"   = c(-1,0), "right"  = c(1,0), "top"    = c(0,1), "bottom" = c(0,-1))
  
  center_outer <- inner_center + dir_vec * center_distance
  
  # Inner circle
  df_inner <- circle_df(inner_center, R_inner)
  df_inner$set <- paste0(id,"_",inner_set)
  df_inner$venn_id <- id
  df_inner$R <- R_inner
  df_inner$cx <- inner_center[1]
  df_inner$cy <- inner_center[2]
  df_inner$is_inner <- TRUE
  
  # Outer circle
  df_outer <- circle_df(center_outer, R_outer)
  df_outer$set <- paste0(id,"_",outer_set)
  df_outer$venn_id <- id
  df_outer$R <- R_outer
  df_outer$cx <- center_outer[1]
  df_outer$cy <- center_outer[2]
  df_outer$is_inner <- FALSE
  
  rbind(df_inner, df_outer)
}

chunk <- function(x, n) {split(seq_along(1:x), ceiling(seq_along(1:x)/n))}

circle_df <- function(center = c(0,0), radius = 1, npoints = 100) {
  theta <- seq(0, 2*pi, length.out = npoints)
  data.frame(x = center[1] + radius * cos(theta),
             y = center[2] + radius * sin(theta))
}

compute_labels <- function(df, id, N1, N2, Noverlap, layout){
  
  require(dplyr)
  
  df_inner <- df[df$is_inner, ]
  df_outer <- df[!df$is_inner, ]
  
  cx_inner <- unique(df_inner$cx)
  cy_inner <- unique(df_inner$cy)
  R_inner <- unique(df_inner$R)
  
  cx_outer <- unique(df_outer$cx)
  cy_outer <- unique(df_outer$cy)
  R_outer <- unique(df_outer$R)
  
  # Inner circle labels
  theta_intersection <- switch(layout, "left" = pi/2, "right" = pi/2, "top" = pi/2, "bottom" = 3*pi/2)
  theta_inner_only <- switch(layout, "left" = 3*pi/2, "right" = 3*pi/2, "top" = 3*pi/2, "bottom" = pi/2)
  
  x_inter <- cx_inner
  y_inter <- cy_inner
  
  x_inner_only <- cx_inner + R_inner * cos(theta_inner_only)
  y_inner_only <- cy_inner + R_inner * sin(theta_inner_only)
  
  # Outer-only label
  theta_outer <- switch(layout, "left" = pi, "right" = 0, "top" = pi/2, "bottom" = 3*pi/2)
  x_outer <- cx_outer 
  y_outer <- cy_outer + R_outer
  y_outer_bottom <- cy_outer - R_outer
  
  hjust <- switch(layout, "left"=0.5, "right"=0.5, "top"=0.5, "bottom"=0.5)
  vjust <- switch(layout, "left"=1.5, "right"=1.5, "top"=1.5, "bottom"=-0.5)
  
  center <- data.frame(venn_id = "center", region = "center", value = N2, x = 0, y = 0, 
                       hjust = 0.5, vjust = 0.5, label = as.character(N2))
  
  # Assemble
  labs <- data.frame(
    venn_id = id,
    region = c("inner_only","outer_only","intersection"),
    value = c(min(N1,N2)-Noverlap, max(N1,N2)-Noverlap, Noverlap),
    x = c(x_inner_only, x_outer, x_inter),
    y = c(y_inner_only, y_outer, y_inter),
    y_bottom = c(y_inner_only, y_outer_bottom, y_inter),
    hjust = c(0.5, hjust, 0.5),
    vjust = c(0.5, vjust, 0.5))
  
  # Remove zero labels
  labs$label <- as.character(labs$value)
  labs$label[labs$value == "0"] <- ""

  labs <- labs %>%
    filter(region != "inner_only") %>%
    bind_rows(center) %>%
    mutate(size = ifelse(region == "intersection", 3, 4)) %>%
    mutate(y = ifelse(venn_id == "bottom", y_bottom, y)) %>%
    select(-y_bottom)
  
  labs
}

find_in_master <- function(x, file = NULL, master_type = "PT") {
  
  require(dplyr)
  
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

gene_countR <- function(pubtator, ortholog = TRUE, variant = TRUE, pa = NULL, pa_id = NULL,
                        hs_genes = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/Nov_5_2025/gene_info_human.txt", 
                        ortholog_db = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/gene_orthologs_MOD.txt",
                        variant_db = "C:/Users/kenny/Desktop/MOODS/projects/team/GMPH/PubTator/Nov_5_2025/LV_data_1_RS.txt") {
  # load dependencies
  lapply(c("data.table", "dplyr", "tidyr"), require, character.only = TRUE)
  
  hs_genes_ncbi <- fread(hs_genes, sep = "\t") %>%
    select(GeneID, Symbol, Synonyms, chromosome)
  
  if (is.character(pubtator)) {
    df <- fread(pubtator, sep = "\t")
  } else {
    df <- pubtator
  }
  
  message("computing gene counts")
  df_gene <- df %>%
    separate_longer_delim(cols = "Genes", delim = "~") %>% 
    filter(!is.na(Genes)) %>%
    mutate(Genes = suppressWarnings(as.integer(Genes))) %>%
    filter(Genes %in% hs_genes_ncbi$GeneID) %>%
    left_join(hs_genes_ncbi, by = join_by(Genes == GeneID)) %>%
    select(Genes, Symbol, Synonyms, PMID) %>%
    distinct(Genes, PMID, .keep_all = TRUE) %>%
    group_by(Genes) %>%
    summarize(Symbol = paste_unique(Symbol), Synonyms = paste_unique(Synonyms), 
              PMID_human = paste_unique(PMID), n_human = n(), .groups = "drop") %>%
    distinct() %>%
    arrange(desc(n_human)) %>%
    mutate(rank_human = dense_rank(desc(n_human))) %>%
    relocate(PMID_human, .after = last_col())
  
  if (variant) {
    variant_df <- fread(variant_db, sep = "\t") %>%
      select(rsid, gene, hgvs, data_chromosome_base_position) %>%
      separate_wider_delim(data_chromosome_base_position, delim = ":", names = c("chr", "pos"), too_few = "align_start") %>%
      filter(!is.na(rsid) & rsid != "" & !is.na(gene) & gene != "" & !is.na(chr) & chr != "") %>% 
      group_by(rsid) %>%
      reframe(gene = paste_unique(gene), hgvs = paste_unique(hgvs), chr = paste_unique(chr))
    
    message("computing variant gene counts")
    df_vars <- df %>%
      separate_longer_delim(cols = "Variants", delim = "~") %>%
      filter(!is.na(Variants)) %>%
      separate_wider_delim(cols = "Variants", delim = "$", names = c("varID", "var"), too_few = "align_start") %>%
      mutate(Genes = as.integer(case_when(grepl("CorrespondingGene:[0-9]", varID) ~  str_extract(varID, "(?<=CorrespondingGene:)[0-9]+"),
                                          grepl("#[0-9]+#", varID) ~ str_extract(varID, "(?<=#)[0-9]+(?=#)"), 
                                          TRUE ~ NA)),
             varID = ifelse(grepl("HGVS", varID),  
                            paste0(sub(".*Gene:([^;]+).*", "\\1", varID), "#",
                                   sub(".*HGVS:([^;]+).*", "\\1", varID)), 
                            varID)) %>%
      left_join(variant_df, by = join_by(varID == rsid)) %>%
      rename(Symbol_LitVar = gene) %>%
      separate_longer_delim(cols = "Symbol_LitVar", delim = ";") %>%
      left_join(hs_genes_ncbi, by = join_by(Genes == GeneID)) %>%
      left_join(hs_genes_ncbi, by = join_by(Symbol_LitVar == Symbol, chr == chromosome)) %>%
      mutate(Genes = ifelse(is.na(Genes), GeneID, Genes), 
             Symbol = ifelse(is.na(Symbol), Symbol_LitVar, Symbol),
             Synonyms = ifelse(is.na(Synonyms.x), Synonyms.y, Synonyms.x)) %>% 
      select(PMID, Genes, Symbol, Synonyms, varID) %>% 
      filter(!is.na(Genes)) %>%      
      distinct() %>%
      group_by(Genes) %>%
      summarize(n_var = n(), PMID_var = paste_unique(PMID, ), id_var = paste_unique(varID), 
                Symbol = paste_unique(Symbol), Synonyms = paste_unique(Synonyms),
                .groups = "drop") %>%
      arrange(desc(n_var)) %>%
      mutate(rank_var = dense_rank(desc(n_var)))
    
    df_gene <- df_gene %>%
      full_join(df_vars, by = c("Genes", "Symbol", "Synonyms"))
  }
  
  # gene ortholog counts
  if (ortholog) {
    if (!exists("ortholog_df")) {
      ortho <- fread(ortholog_db)
    }
    
    message("computing orthologous gene counts")
    df_ortholog <- df %>%
      separate_longer_delim(cols = "Genes", delim = "~") %>% 
      filter(!is.na(Genes)) %>%
      mutate(Genes = suppressWarnings(as.integer(Genes))) %>%
      filter(!Genes %in% hs_genes_ncbi$GeneID) %>%
      left_join(ortho, by = join_by(Genes == ortholog_gene)) %>%
      left_join(hs_genes_ncbi, by = join_by(hs_gene == GeneID)) %>%
      mutate(Genes = if_else(!is.na(hs_gene), hs_gene, Genes)) %>%
      group_by(Genes) %>%
      summarize(n_ortho = n(), genes_ortho = paste_unique(paste0(Genes, "(", ortholog_tax, ")")),
                PMID_ortho = paste_unique(PMID), Symbol = paste_unique(Symbol), Synonyms = paste_unique(Synonyms),
                .groups = "drop") %>%
      select(Genes, Symbol, Synonyms, n_ortho, genes_ortho, PMID_ortho) %>%
      distinct() %>%
      mutate(rank_ortho = dense_rank(desc(n_ortho)))
    
    df_gene <- df_gene %>%
      full_join(df_ortholog, by = c("Genes", "Symbol", "Synonyms"))
  }
  
  df_gene <- df_gene %>%
    rowwise() %>%
    mutate(PMID_all = case_when(!is.na(PMID_human) & !is.na(PMID_ortho) & !is.na(PMID_var) ~ 
                                  paste_unique(union(union(strsplit(PMID_human, ";")[[1]],
                                                           strsplit(PMID_ortho, ";")[[1]]),
                                                     strsplit(PMID_var, ";")[[1]])),
                                !is.na(PMID_human) & !is.na(PMID_ortho) & is.na(PMID_var) ~ 
                                  paste_unique(union(strsplit(PMID_human, ";")[[1]], 
                                                     strsplit(PMID_ortho, ";")[[1]])),
                                !is.na(PMID_human) & is.na(PMID_ortho) & !is.na(PMID_var) ~ 
                                  paste_unique(union(strsplit(PMID_human, ";")[[1]], 
                                                     strsplit(PMID_var, ";")[[1]])),
                                is.na(PMID_human) & !is.na(PMID_ortho) & !is.na(PMID_var) ~ 
                                  paste_unique(union(strsplit(PMID_ortho, ";")[[1]], 
                                                     strsplit(PMID_var, ";")[[1]])),
                                !is.na(PMID_human) & is.na(PMID_ortho) & is.na(PMID_var) ~ PMID_human,
                                is.na(PMID_human) & !is.na(PMID_ortho) & is.na(PMID_var) ~ PMID_ortho,
                                is.na(PMID_human) & is.na(PMID_ortho) & !is.na(PMID_var) ~ PMID_var),
           n_total = case_when(!is.na(PMID_human) & !is.na(PMID_ortho) & !is.na(PMID_var) ~ 
                                 length(union(union(strsplit(PMID_human, ";")[[1]],
                                                    strsplit(PMID_ortho, ";")[[1]]),
                                              strsplit(PMID_var, ";")[[1]])),
                               !is.na(PMID_human) & !is.na(PMID_ortho) & is.na(PMID_var) ~ 
                                 length(union(strsplit(PMID_human, ";")[[1]], 
                                              strsplit(PMID_ortho, ";")[[1]])),
                               !is.na(PMID_human) & is.na(PMID_ortho) & !is.na(PMID_var) ~ 
                                 length(union(strsplit(PMID_human, ";")[[1]], 
                                              strsplit(PMID_var, ";")[[1]])),
                               is.na(PMID_human) & !is.na(PMID_ortho) & !is.na(PMID_var) ~ 
                                 length(union(strsplit(PMID_ortho, ";")[[1]], 
                                              strsplit(PMID_var, ";")[[1]])),
                               !is.na(PMID_human) & is.na(PMID_ortho) & is.na(PMID_var) ~ n_human,
                               is.na(PMID_human) & !is.na(PMID_ortho) & is.na(PMID_var) ~ n_ortho,
                               is.na(PMID_human) & is.na(PMID_ortho) & !is.na(PMID_var) ~ n_var)) %>%
    ungroup() %>%
    arrange(desc(n_total)) %>%
    mutate(rank_total = dense_rank(desc(n_total)))
  
  # panel app comparison
  if (is.null(pa)) {
    return(list(gene_count = df_gene))
  } else {
    panel <- subset(pa, pa$id == pa_id)
    panel_genes <- panel$entrez_id
    
    panel_genes_green <- subset(panel, color == "Green")$entrez_id
    panel_genes_amber <- subset(panel, color == "Amber")$entrez_id
    panel_genes_red <- subset(panel, color == "Red")$entrez_id
    
    present <- panel_genes[panel_genes %in% unique(df_gene$Genes)]
    present_green <- panel_genes_green[panel_genes_green %in% unique(df_gene$Genes)]
    present_amber <- panel_genes_amber[panel_genes_amber %in% unique(df_gene$Genes)]
    present_red <- panel_genes_red[panel_genes_red %in% unique(df_gene$Genes)]
    
    df_gene <- df_gene %>%
      left_join(panel[,c("entrez_id", "color")], by = join_by(Genes == entrez_id)) %>%
      mutate(color = tolower(color)) %>%
      select(Genes, Symbol, Synonyms, color, n_total, rank_total, 
             n_human, rank_human, PMID_human,
             n_var, rank_var, id_var, PMID_var,
             n_ortho, rank_ortho, genes_ortho, PMID_ortho)
    
    missing <- panel %>%
      filter(entrez_id %in% panel_genes[!panel_genes %in% unique(df_gene$Genes)]) %>%
      pull(entity_name)
    missing_green <- panel %>%
      filter(entrez_id %in% panel_genes_green[!panel_genes_green %in% unique(df_gene$Genes)]) %>%
      pull(entity_name)
    missing_amber <- panel %>%
      filter(entrez_id %in% panel_genes_amber[!panel_genes_amber %in% unique(df_gene$Genes)]) %>%
      pull(entity_name)
    missing_red <- panel %>%
      filter(entrez_id %in% panel_genes_red[!panel_genes_red %in% unique(df_gene$Genes)]) %>%
      pull(entity_name)
    
    extra <- df_gene[!df_gene$Genes %in% panel_genes, ]$Symbol
    
    missing_genes <- paste(sort(missing), collapse = ";")
    missing_genes_green <- paste(sort(missing_green), collapse = ";")
    missing_genes_amber <- paste(sort(missing_amber), collapse = ";")
    missing_genes_red <- paste(sort(missing_red), collapse = ";")
    
    pct_present <- length(present)/length(panel_genes)
    pct_present_green <- length(present_green)/length(panel_genes_green)
    pct_present_amber <- length(present_amber)/length(panel_genes_amber)
    pct_present_red <- length(present_red)/length(panel_genes_red)
    
    precision <- length(present)/nrow(df_gene)
    
    df_stats <- data.frame(stat = c("n_panel", "n_genes", "n_present", "tpr", "tpr_green", "tpr_amber", "tpr_red", 
                                    "n_extra", "precision", "n_articles"),
                           value = c(length(panel_genes), nrow(df_gene), length(present), pct_present, pct_present_green, 
                                     pct_present_amber, pct_present_red, length(extra), precision, length(df$PMID)))
    
    return(list(panel = unique(panel$panel), panelapp_genes = sort(panel_genes),  gene_count = df_gene, 
                stats = df_stats, missing = missing, missing_green = missing_green, missing_amber = missing_amber, 
                missing_red = missing_red))
  }
}

intersection_area <- function(d, R1, R2) {
  if(d >= R1 + R2) return(0)
  if(d <= abs(R1 - R2)) return(pi * min(R1,R2)^2)
  r1sq <- R1^2
  r2sq <- R2^2
  part1 <- r1sq * acos((d^2 + r1sq - r2sq)/(2*d*R1))
  part2 <- r2sq * acos((d^2 + r2sq - r1sq)/(2*d*R2))
  part3 <- 0.5 * sqrt((-d+R1+R2)*(d+R1-R2)*(d-R1+R2)*(d+R1+R2))
  part1 + part2 - part3
}

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

paste_unique <- function(x, sep = ";") {
  if (all(is.na(x)) | all(x == "")) {
    return(NA)
  } else {
    return(paste(unique(x[!is.na(x)]), collapse = sep))
  }
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
    
    df_gene <- data$genes$gene_data
    df$alias <- sapply(df_gene$alias, paste, collapse = "~")
    df$biotype <- df_gene$biotype
    df$hgnc_id <- df_gene$hgnc_id
    df$gene_name <- df_gene$gene_name
    df$omim_gene <- sapply(df_gene$omim_gene, paste, collapse = "~")
    df$alias_name <- sapply(df_gene$alias_name, paste, collapse = "~")
    df$gene_symbol <- df_gene$gene_symbol
    df$hgnc_release <- df_gene$hgnc_release
    
  } else {
    df <- data.frame(panel = data$name, id = data$id, disease_group = data$disease_group,
                     disease_sub_group = data$disease_sub_group, version = data$version)
  }
  return(df)
}

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

plot_plus_venn <- function(res) {
  
  require(ggplot2)
  
  ggplot() +
    geom_polygon(data = res$polygons,
                 aes(x, y, group = interaction(venn_id, set), fill = set), alpha = 0.4) +
    geom_text(data = res$labels,
              aes(x, y, label = label, hjust = hjust, vjust = vjust, size = 3)) +
    scale_size_identity() +
    coord_equal() +
    theme_void() +
    facet_wrap(~disease) + 
    scale_fill_manual(values=c("top_A"="#E58EDD","top_B"="#007A82",
                               "bottom_A"="#E1C420","bottom_B"="#007A82",
                               "left_A"="#318E34","left_B"="#007A82",
                               "right_A"="#E55B4A","right_B"="#007A82")) +
    theme(legend.position = "none", strip.text = element_text(size = 17.5))
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
