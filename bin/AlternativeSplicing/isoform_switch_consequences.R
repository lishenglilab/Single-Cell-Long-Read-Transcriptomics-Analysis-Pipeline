#!/usr/bin/env Rscript
# Isoform Switch Consequences Analysis

suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(ggplot2)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("-i", "--switchlist_file"), type="character", default=NULL,
              help="Path to switchList RDS file", metavar="FILE"),
  make_option(c("-o", "--output_dir"), type="character", default="./IsoformSwitchConsequences",
              help="Output directory [default: %default]", metavar="DIR"),
  make_option(c("--dIF_cutoff"), type="numeric", default=0.1,
              help="dIF cutoff for significant isoform switches [default: %default]", metavar="NUM"),
  make_option(c("--alpha"), type="numeric", default=0.05,
              help="Alpha level for significance [default: %default]", metavar="NUM")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$switchlist_file)) {
  print_help(opt_parser)
  stop("Missing required argument: --switchlist_file")
}

# Assign variables
switchlist_file <- opt$switchlist_file
output_dir <- opt$output_dir
dIF_cutoff <- opt$dIF_cutoff
alpha <- opt$alpha

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to log messages
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(paste0("[", timestamp, "] ", msg))
  write(paste0("[", timestamp, "] ", msg), 
        file = file.path(output_dir, "consequences_analysis.log"), append = TRUE)
}

log_message("Starting Isoform Switch Consequences Analysis...")
log_message(paste("Switch list file:", switchlist_file))
log_message(paste("Output directory:", output_dir))
log_message(paste("dIF cutoff:", dIF_cutoff))
log_message(paste("Alpha:", alpha))

# ----------------------------------------------------------------------
# Load switch list
# ----------------------------------------------------------------------
log_message("Loading switch list...")

if (!file.exists(switchlist_file)) {
  stop(paste("Switch list file not found:", switchlist_file))
}

aSwitchList <- readRDS(switchlist_file)
log_message(paste("Loaded switch list with", nrow(aSwitchList$isoformFeatures), "isoform features"))

# ----------------------------------------------------------------------
# Analyze consequences
# ----------------------------------------------------------------------
log_message("Analyzing consequences...")

aSwitchList <- analyzeIntronRetention(aSwitchList, onlySwitchingGenes = FALSE)
log_message("Intron retention analysis completed")

aSwitchList <- isoformSwitchTestDEXSeq(aSwitchList, dIFcutoff = dIF_cutoff, alpha = alpha)
log_message("Isoform switch test completed")

# Analyze switch consequences
aSwitchList <- analyzeSwitchConsequences(aSwitchList)
log_message("Switch consequences analysis completed")

# Save results
log_message("Saving results...")
saveRDS(aSwitchList, file.path(output_dir, paste0("switchList_consequences.rds")))

# Save switch consequences table
write.table(aSwitchList$switchConsequence,
            file.path(output_dir, "switch_consequences.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------------------------------------------------
# Generate plots
# ----------------------------------------------------------------------
log_message("Generating plots...")

# Consequence summary plot
p_consequence_summary <- extractConsequenceSummary(
  aSwitchList,
  consequencesToAnalyze = "all",
  plotGenes = FALSE,
  asFractionTotal = FALSE
)

ggsave(plot = p_consequence_summary,
       filename = file.path(output_dir, "consequence_summary.pdf"),
       width = 10, height = 6)
log_message("Saved consequence summary plot")

# Consequence enrichment plot
p_consequence_enrichment <- extractConsequenceEnrichment(
  aSwitchList,
  consequencesToAnalyze = "all",
  analysisOppositeConsequence = TRUE,
  localTheme = theme_bw(base_size = 14),
  returnResult = FALSE,
  minEventsForPlotting = 1
)

ggsave(plot = p_consequence_enrichment,
       filename = file.path(output_dir, "consequence_enrichment.pdf"),
       width = 10, height = 12)
log_message("Saved consequence enrichment plot")

# Volcano plot
p_volcano <- ggplot(data = aSwitchList$isoformFeatures, 
                    aes(x = dIF, y = -log10(isoform_switch_q_value))) +
  geom_point(
    aes(color = abs(dIF) > dIF_cutoff & isoform_switch_q_value < alpha),
    size = 1
  ) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  geom_vline(xintercept = c(-dIF_cutoff, dIF_cutoff), linetype = "dashed") +
  facet_wrap(~ condition_2) +
  scale_color_manual("Significant\nIsoform Switch", values = c("black", "red")) +
  labs(x = "dIF", y = "-Log10 (Isoform Switch Q Value)") +
  theme_bw()

ggsave(plot = p_volcano,
       filename = file.path(output_dir, "volcano_plot.pdf"),
       width = 15, height = 8)
log_message("Saved volcano plot")

# Example switch plot for top gene
if (nrow(aSwitchList$switchConsequence) > 0) {
  top_gene <- unique(aSwitchList$switchConsequence$gene_name)[1]
  pdf(file.path(output_dir, "example_switch_plot.pdf"),width = 20, height = 15)
  switchPlot(
    aSwitchList,
    gene = top_gene,
    localTheme = theme_bw(base_size = 13),
    plotTopology = F
  )
  dev.off()
  
  log_message(paste("Saved example switch plot for gene:", top_gene))
} else {
  log_message("No switch consequences found, skipping example switch plot")
}

# ----------------------------------------------------------------------
# Create summary statistics
# ----------------------------------------------------------------------
log_message("Creating summary statistics...")

summary_stats <- list(
  Total_Isoforms_Analyzed = nrow(aSwitchList$isoformFeatures),
  Total_Genes_Analyzed = length(unique(aSwitchList$isoformFeatures$gene_id)),
  Significant_Isoform_Switches = sum(aSwitchList$isoformFeatures$isoform_switch_q_value < alpha & 
                                       abs(aSwitchList$isoformFeatures$dIF) > dIF_cutoff, na.rm = TRUE),
  Genes_with_Switches = length(unique(aSwitchList$isoformFeatures$gene_id[
    aSwitchList$isoformFeatures$isoform_switch_q_value < alpha & 
      abs(aSwitchList$isoformFeatures$dIF) > dIF_cutoff
  ])),
  dIF_Cutoff = dIF_cutoff,
  Alpha = alpha,
  Analysis_Date = as.character(Sys.time())
)

# Write summary to file
summary_file <- file.path(output_dir, "consequences_summary.txt")
sink(summary_file)
cat("Isoform Switch Consequences Analysis Summary\n")
cat("============================================\n\n")
for (name in names(summary_stats)) {
  cat(paste0(name, ": ", summary_stats[[name]], "\n"))
}
sink()

# Save session info
writeLines(capture.output(sessionInfo()), 
           file.path(output_dir, "consequences_session_info.txt"))

writeLines(capture.output(opt), 
           file.path(output_dir, "consequences_command_args.txt"))

log_message("Isoform switch consequences analysis completed successfully!")
log_message(paste("Results saved in:", output_dir))