##### Generating a synthetic gene-expression data using simphony #####
library(simphony)
library(data.table)

set.seed(42)

# defining the characteristics
num_patients <- 5
num_genes <- 10
num_samples <- 7
window_hours <- 24

# create 7 timepoints evenly distributed across a 24-hour window (0, 4, 8, 12, 16, 20, 24)
t_points <- seq(0, window_hours, length.out = num_samples)

# define patients as independent conditions
patient_list <- list()

# expression base range
main_baselines <- runif(num_genes, min = 3, max = 11)

# feature df for patients
for (i in 1:num_patients) {

  # first 5 genes: log2 amplitude between 0.1 and 0.49
  amp_low <- sample(seq(0.1, 0.49, length.out = 5))

  # next 5 genes: log2 amplitude between 0.5 and 2.0
  amp_high <- sample(seq(0.5, 2.0, length.out = 5))

  # combine them to make the full 10-gene amplitude vector
  amplitudes <- c(amp_low, amp_high)

  # build the feature matrix for this patient
  patient_list[[paste0("Patient_", i)]] <- data.table(
    fracFeatures = rep(1 / num_genes, num_genes),
    base = main_baselines,
    amp = amplitudes,
    period = window_hours
  )
}

# regulating count noise
low_noise_disp <- function(x) {0.05 * defaultDispFunc(x)}

# simulate the dataset
sim_data <- simphony(
  featureGroupsList = patient_list,
  nFeatures = num_genes,
  timepointsType = "specified",
  timepoints = t_points,
  nReps = 1,
  family = "negbinom",
  dispFunc = low_noise_disp
)

# the generated read counts
expr_matrix <- sim_data$abundData
# log2 transform the counts
simdf <- log2(expr_matrix)

# sample mapping
sample_meta <- sim_data$sampleMetadata

# gene properties
feature_meta <- sim_data$featureMetadata

# defining the column names according to what rhythmCI expects
colnames(simdf) <- c('P1_0','P1_4','P1_8','P1_12','P1_16','P1_20','P1_24',
                     'P2_0','P2_4','P2_8','P2_12','P2_16','P2_20','P2_24',
                     'P3_0','P3_4','P3_8','P3_12','P3_16','P3_20','P3_24',
                     'P4_0','P4_4','P4_8','P4_12','P4_16','P4_20','P4_24',
                     'P5_0','P5_4','P5_8','P5_12','P5_16','P5_20','P5_24')

simdf <- t(simdf)

usethis::use_data(simdf, overwrite = TRUE)

##### RCC files from GEO #####
library(GEOquery)

# GEO accession
ftp_base <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146530/suppl/"
geo_id <- "GSE146530"

# Destination directory under inst/extdata
dest_dir <- file.path("inst", "extdata", geo_id)

file <- "GSE146530_RAW.tar"

# Create directory if it doesn't exist
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

# Download each file
for (f in file) {
  url <- paste0(ftp_base, f)
  dest_file <- file.path(dest_dir, f)
  download.file(url, dest_file, mode = "wb")
}

# List all .gz files
gz_files <- utils::untar(file.path(dest_dir, file), list = TRUE)
# gz_files <- list.files(file.path(dest_dir, file), pattern = "\\.gz$", full.names = TRUE)

# Extract the contents to the destination folder
utils::untar(file.path(dest_dir, file), exdir = dest_dir)

# Delete the .tar file after successful extraction
file.remove(file.path(dest_dir, file))

# Loop over each .gz file, extract, and delete the original
for (gz_file in gz_files) {
  R.utils::gunzip(file.path(dest_dir, gz_file), overwrite = TRUE, remove = TRUE)
}

# List all .RCC files and delete unwanted
rcc_files <- list.files(dest_dir, pattern = "\\.RCC$", full.names = TRUE)
sel_rcc_files <- rcc_files[166:189]
rcc_to_delete <-rcc_files[!rcc_files %in% sel_rcc_files]
file.remove(rcc_to_delete)

# List the file names of remaining .RCC files without path
sel_rcc_file_name <- list.files(dest_dir, pattern = "\\.RCC$", full.names = FALSE)

# Create the meta data
meta_df <- data.frame(files.RCC=sel_rcc_file_name, stringsAsFactors = FALSE)

meta_df$id <- rep(c("f1", "f2", "f3"), length.out = nrow(meta_df))
meta_df$time <- rep(seq(0, 28, by = 4), each = 3)
meta_df$ext_time <- rep((seq(8, by = 4, length.out = 8) %% 24), each = 3)
meta_df$group <- 'group1'

write.csv(meta_df, file = "inst/extdata/meta_GSE146530.csv", row.names = FALSE)

