# ============================================================================
# STANDALONE TEST SCRIPT PER server-annotation.R
# ============================================================================

# === LIBRARIES ===
library(shiny)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(IRanges)
library(Rsamtools)
library(condformat)
# Load required libraries
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(Rsamtools)
library(condformat)
# library(uncoverappLib) # If you have this package

# Mock Shiny inputs
input <- list(
  coverage_co = 30,
  UCSC_Genome = "hg19", # or "hg38" to test both
  Chromosome = "chr1",
  query_Database = "1:1902000-1902200" # Example query region
)

# File paths - UPDATE THESE PATHS
bed_file_path <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_fake191.bed.gz"
annotation_file_path <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_hg19.bed.gz"

# Mock mysample() function - replace with your actual data loading
mysample <- function() {
  # Read your BED file here
  # Assuming it has columns: chromosome, start, end, coverage, etc.
  bed_data <- read.table(bed_file_path, header = FALSE, sep = "\t")
  # Adjust column names based on your actual file structure
  colnames(bed_data) <- c("chromosome", "start", "end", "coverage") # Add more columns as needed
  return(bed_data)
}

# Mock getAnnotationFiles() function
getAnnotationFiles <- function() {
  return(annotation_file_path)
}

# Convert reactive functions to regular functions
filtered_low_nucl <- function() {
  sample_data <- mysample()
  if (is.null(sample_data))
    return(NULL)
  
  sample_data %>%
    dplyr::filter(chromosome == input$Chromosome,
                  coverage <= as.numeric(input$coverage_co))
}

intBED <- function() {
  filtered_data <- filtered_low_nucl()
  if (is.null(filtered_data))
    return(NULL)
  
  bedA <- filtered_data
  file.name <- getAnnotationFiles()
  
  query <- input$query_Database
  query.regions <- read.table(text = gsub("[:-]+", " ", query, perl = TRUE),
                             header = FALSE, col.names = c("chr", "start", "end"))
  
  if (is.null(query.regions))
    return(NULL)
  
  print("Query regions:")
  print(query.regions)
  
  result <- try({
    fq <- GenomicRanges::makeGRangesFromDataFrame(query.regions, keep.extra.columns = TRUE)
    res <- Rsamtools::scanTabix(file.name, param = fq)
    sapply(res, length)
    dff <- Map(function(elt) {
      read.csv(textConnection(elt), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    }, res)
    bedB <- as.data.frame(dff)
  })
  
  if ("try-error" %in% class(result)) {
    err_msg <- 'no coordinates recognized'
    cat("Error:", err_msg, "\n")
    return(NULL)
  }
  
  print("Head of bedB:")
  print(head(bedB))
  
  # Set column names (same for both genome builds)
  colnames(bedB) <- c('Chromo', 'start', 'end', 'REF', 'ALT',
                      'dbsnp', 'GENENAME', 'PROTEIN_ensembl',
                      'MutationAssessor', 'SIFT', 'Polyphen2',
                      'M_CAP', 'CADD_PHED', 'AF_gnomAD', 'ClinVar',
                      'clinvar_MedGen_id', 'clinvar_OMIM_id', 'HGVSc_VEP',
                      'HGVSp_VEP')
  
  str(bedB)
  
  # Add chromosome prefix
  for (i in bedB[1]) {
    Chromosome <- paste("chr", i, sep = "")
    bedB <- cbind(Chromosome, bedB)
    bedB[,2] <- NULL
  }
  
  # Convert data types
  bedB$Chromosome <- as.character(bedB$Chromosome)
  bedB$AF_gnomAD <- suppressWarnings(as.numeric(bedB$AF_gnomAD))
  bedB$CADD_PHED <- suppressWarnings(as.numeric(bedB$CADD_PHED))
  
  # Function to intersect bed files
  intersectBedFiles.GR <- function(bed1, bed2) {
    bed1_gr <- makeGRangesFromDataFrame(bed1, ignore.strand = TRUE,
                                       keep.extra.columns = TRUE)
    bed2_gr <- makeGRangesFromDataFrame(bed2, ignore.strand = TRUE,
                                       keep.extra.columns = TRUE)
    tp <- findOverlaps(query = bed1_gr, subject = bed2_gr, type = "any")
    intersect_df <- data.frame(bed1_gr[queryHits(tp),], bed2_gr[subjectHits(tp),])
    return(intersect_df)
  }
  
  intersect_df <- intersectBedFiles.GR(bedA, bedB)
  return(intersect_df)
}

condform_table <- function() {
  intersect_data <- intBED()
  if (is.null(intersect_data))
    return(NULL)
  
  cat("Creating formatted table...\n")
  
  formatted_table <- condformat(intersect_data) %>%
    rule_fill_discrete(ClinVar, expression = ClinVar != ".",
                       colours = c("TRUE" = "red", "FALSE" = "green")) %>%
    rule_fill_discrete(CADD_PHED, expression = CADD_PHED > 20,
                       colours = c("TRUE" = "red", "FALSE" = "green")) %>%
    rule_fill_discrete(MutationAssessor, expression = MutationAssessor == 'H',
                       colours = c("TRUE" = "red", "FALSE" = "green")) %>%
    rule_fill_discrete(M_CAP, expression = M_CAP == 'D',
                       colours = c("TRUE" = "red", "FALSE" = "green")) %>%
    rule_fill_discrete(AF_gnomAD, expression =
                         ifelse(is.na(AF_gnomAD) | AF_gnomAD < 0.5,
                                'TRUE', 'FALSE'),
                       colours = c("TRUE" = "red", "FALSE" = "green")) %>%
    rule_fill_discrete(c(start, end),
                       expression = grepl("H|M", MutationAssessor) &
                         ClinVar != "." & AF_gnomAD < 0.5,
                       colours = c("TRUE" = "yellow", "FALSE" = "")) %>%
    rule_css(c(start, end),
             expression = ifelse(grepl("H|M", MutationAssessor) &
                                   ClinVar != "." & AF_gnomAD < 0.5,
                                 "red", "green"),
             css_field = "color")
  
  return(formatted_table)
}

# Run the analysis
cat("Starting analysis...\n")

# Get the results
final_table <- condform_table()

# Save the results
if (!is.null(final_table)) {
  # Save as HTML (condformat creates HTML tables)
  output_file <- paste0("analysis_results_", Sys.Date(), ".html")
  
  # You can also save the raw data
  raw_data <- intBED()
  if (!is.null(raw_data)) {
    write.csv(raw_data, paste0("raw_analysis_results_", Sys.Date(), ".csv"), row.names = FALSE)
    cat("Raw data saved to:", paste0("raw_analysis_results_", Sys.Date(), ".csv"), "\n")
  }
  
  cat("Analysis completed successfully!\n")
  cat("Formatted table object created and available as 'final_table'\n")
  
  # Print summary
  cat("Summary:\n")
  cat("- Rows in final result:", nrow(raw_data), "\n")
  cat("- Columns:", ncol(raw_data), "\n")
  
} else {
  cat("No results found. Check your input parameters and file paths.\n")
}

# To view the formatted table (if in RStudio or similar environment):
# final_table
