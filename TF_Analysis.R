###############################################################################
# tf_enrichment_lib.R
# A small R helper library for discovering transcription factors (TFs) that
# regulate an arbitrary set of “hub” genes and for comparing their expression
# across two bulk-RNA-seq conditions.
#
# Author: Abtin Tondar
###############################################################################

# Required packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(httr)          # API calls (CheA3)
  library(jsonlite)      # API response parsing
  library(data.table)    # Efficient tables
  library(ggplot2)       # Visualisation
})

###############################################################################
# 1.  CHEA3 ENRICHMENT
###############################################################################
run_chea3 <- function(gene_vec) {
  url <- "https://amp.pharm.mssm.edu/chea3/api/enrich"
  body <- list(query_name = "hub_gene_query", gene_set = gene_vec)
  res  <- httr::POST(url, body = body, encode = "json")
  if (res$status_code != 200) stop("CHEA3 API call failed")
  parsed <- jsonlite::fromJSON(rawToChar(res$content))
  # We keep the Integrated_Rank sheet, which combines all evidence sources
  tf_tab <- as.data.table(parsed[["Integrated_Rank"]])
  setnames(tf_tab, c("tf", "rank", "p_value", "z_score"))
  return(tf_tab[])
}

###############################################################################
# 2.  SELECT TOP TFs
###############################################################################
select_top_tfs <- function(tf_dt, top_percent = 10) {
  keep_n <- ceiling(nrow(tf_dt) * top_percent / 100)
  tf_dt[rank <= keep_n]
}

###############################################################################
# 3.  MATCH WITH DEG RESULTS
###############################################################################
intersect_degs_with_tfs <- function(deg_dt,
                                    tf_vec,
                                    gene_col = "symbol",
                                    padj_col = "padj",
                                    padj_cut = 0.05) {
  setnames(deg_dt, gene_col, "symbol_temp")
  res <- deg_dt[symbol_temp %in% tf_vec]
  res <- res[get(padj_col) < padj_cut]
  setnames(res, "symbol_temp", "symbol")
  return(res[])
}

###############################################################################
# 4.  COMPARE TWO CONDITIONS
###############################################################################
compare_two_conditions <- function(cond1_dt, cond2_dt,
                                   logfc1_col, logfc2_col,
                                   padj1_col, padj2_col,
                                   logfc_diff_cut = 0.2) {

  merged <- merge(cond1_dt, cond2_dt,
                  by = "symbol",
                  suffixes = c("_c1", "_c2"))
  merged[, logfc_gap := abs(get(logfc1_col)) - abs(get(logfc2_col))]
  merged[, signif_both :=
           get(padj1_col) < 0.05 & get(padj2_col) < 0.05]
  merged[, flag :=
           fifelse(signif_both & abs(logfc_gap) > logfc_diff_cut,
                   "Stronger in cond1",
                   "Similar or stronger in cond2")]
  return(merged[])
}

###############################################################################
# 5.  BAR PLOT FOR SHARED TFs
###############################################################################
plot_shared_tfs <- function(shared_dt,
                            tf_col = "symbol",
                            logfc1_col, logfc2_col,
                            cond_names = c("Condition-1", "Condition-2"),
                            outfile = NULL) {

  long <- rbind(
    shared_dt[, .(symbol = get(tf_col),
                  logFC = get(logfc1_col),
                  condition = cond_names[1])],
    shared_dt[, .(symbol = get(tf_col),
                  logFC = get(logfc2_col),
                  condition = cond_names[2])]
  )

  gp <- ggplot(long,
               aes(x = reorder(symbol, logFC),
                   y = logFC,
                   fill = condition)) +
        geom_col(position = position_dodge(0.8),
                 colour = "black") +
        coord_flip() +
        labs(x = "Transcription factors",
             y = "log2 fold change",
             fill = NULL) +
        theme_minimal(base_size = 14)

  if (!is.null(outfile))
    ggsave(outfile, gp, width = 10, height = 6, dpi = 400)
  return(gp)
}

###############################################################################
# 6.  ONE-LINER PIPELINE FOR QUICK USE
###############################################################################
tf_enrichment_workflow <- function(hub_genes,
                                   deg1_dt, deg2_dt,
                                   deg1_name = "Condition1",
                                   deg2_name = "Condition2",
                                   logfc1_col, logfc2_col,
                                   padj1_col = "padj",
                                   padj2_col = "padj",
                                   top_percent = 10) {

  message("Calling CHEA3 …")
  tf_raw   <- run_chea3(hub_genes)
  tf_top   <- select_top_tfs(tf_raw, top_percent)
  tf_list  <- tf_top$tf

  message("Intersecting with DEG tables …")
  sig1 <- intersect_degs_with_tfs(deg1_dt, tf_list,
                                  padj_col = padj1_col)
  sig2 <- intersect_degs_with_tfs(deg2_dt, tf_list,
                                  padj_col = padj2_col)

  message("Comparing shared TFs …")
  shared <- compare_two_conditions(sig1, sig2,
                                   logfc1_col, logfc2_col,
                                   padj1_col, padj2_col)

  gp <- plot_shared_tfs(shared,
                        logfc1_col = logfc1_col,
                        logfc2_col = logfc2_col,
                        cond_names = c(deg1_name, deg2_name))

  list(raw_enrichment = tf_raw,
       top_tfs        = tf_top,
       significant_1  = sig1,
       significant_2  = sig2,
       shared_compare = shared,
       plot           = gp)
}

###############################################################################
# End of library
###############################################################################
