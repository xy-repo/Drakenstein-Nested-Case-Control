args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop("Usage: run_maaslin.R features.tsv metadata.tsv output fixed_effects min_prevalence reference_variable reference_level")
}

suppressPackageStartupMessages(library(Maaslin2))
features <- read.delim(args[[1]], row.names = 1, check.names = FALSE)
metadata <- read.delim(args[[2]], row.names = 1, check.names = FALSE)
metadata$case <- factor(metadata$case, levels = c("Case", "Control"))
if ("age_days_npp_collection" %in% names(metadata)) {
    metadata$age_days_npp_collection <- as.numeric(metadata$age_days_npp_collection)
}
if ("antibiotics_use_at_collection_analysis" %in% names(metadata)) {
    metadata$antibiotics_use_at_collection_analysis <- factor(
        metadata$antibiotics_use_at_collection_analysis
    )
}

fixed_effects <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
reference <- NULL
if (nzchar(args[[6]]) && nzchar(args[[7]])) {
    reference <- c(args[[6]], args[[7]])
}

result <- Maaslin2(
    input_data = features,
    input_metadata = metadata,
    output = args[[3]],
    normalization = "TSS",
    transform = "LOG",
    min_prevalence = as.numeric(args[[5]]),
    fixed_effects = fixed_effects,
    reference = reference,
    plot_scatter = FALSE,
    plot_heatmap = FALSE
)
invisible(result)
