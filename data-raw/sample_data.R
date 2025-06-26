##### Generating a synthetic gene-expression data using simphony #####
library(simphony)
library(data.table)

set.seed(7)

# create a synthetic data using varying amplitude
featureGroups <- list(data.table::data.table(amp = c(0,0.05,0.1,0.2,0.4,0.5,0.9,1.2,1.6,2)),
                      data.table::data.table(amp = c(0,0.05,0.1,0.2,0.4,0.5,0.9,1.2,1.6,1.8)))
simData <- simphony::simphony(featureGroups, nFeatures = 20, interval = 4, nReps = 1, family = 'negbinom')

# assigning data to variables
simdf <- as.data.frame(simData$abundData)
simdf_sam_meta <- as.data.frame(simData$sampleMetadata)
simdf_fea_meta <- as.data.frame(simData$featureMetadata[,c(1:3,5:7)])
simdf <- log2(simdf)

colnames(simdf) <- c('P1_0','P1_4','P1_8','P1_12','P1_16','P1_20','P1_24','P1_28','P1_32','P1_36','P1_40','P1_44',
                     'P2_0','P2_4','P2_8','P2_12','P2_16','P2_20','P2_24','P2_28','P2_32','P2_36','P2_40','P2_44')

simdf_img <- simdf #creating a df for image
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

