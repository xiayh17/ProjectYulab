library(httr)

# Data fetch --------------------------------------------------------------
# Set stringDB base URL
address <- "https://string-db.org"

# Validate the address
httr::stop_for_status(httr::GET(address))

# Version info
current_version <- read.table(url(paste(address, "/api/tsv/version", sep = "")), header = TRUE)
available_version <- read.table(url(paste(address, "/api/tsv/available_api_versions", sep = "")), header = TRUE)

#' networkParamsParser
#' 
#' parameters parser for [Getting the STRING network interactions](https://string-db.org/cgi/help.pl?sessionId=btsvnCeNrBk7).
#'
#' @param identifiers required parameter for multiple items, e.g. `c("PTCH1", "TP53", "BRCA1", "BRCA2")`
#' @param species NCBI taxon identifiers (e.g. Human is 9606, see: [STRING organisms](https://string-db.org/cgi/input.pl?input_page_active_form=organisms).
#' @param required_score threshold of significance to include a interaction, a number between 0 and 1000 (default depends on the network)
#' @param network_type network type: functional (default), physical
#' @param add_nodes adds a number of proteins with to the network based on their confidence score (default:1)
#' @param show_query_node_labels when available use submitted names in the preferredName column when (0 or 1) (default:0)
#' @param caller_identity your identifier for us.
#'
#' @return a list contain parameters for query
networkParamsParser <- function(
    identifiers,
    species,
    required_score = NULL,
    network_type = "functional",
    add_nodes = 1,
    show_query_node_labels = 0,
    caller_identity = NULL
) {
  # Format the identifiers
  identifiers <- paste(identifiers, collapse = "\n")
  
  # Check parameters
  if (missing(species)) {
    stop("Please provide an NCBI taxon identifier for the species.")
  }
  
  # Create parameters list
  params <- list(
    identifiers = identifiers,
    species = species,
    required_score = required_score,
    network_type = network_type,
    add_nodes = add_nodes,
    show_query_node_labels = show_query_node_labels,
    caller_identity = caller_identity
  )
  
  # Remove NULL elements from the list
  filtered_params <- Filter(Negate(is.null), params)
  return(filtered_params)
}

networkParams <- networkParamsParser(
  identifiers = c("PTCH1", "TP53", "BRCA1", "BRCA2"),
  species = 9606,
  network_type = "functional",
  add_nodes = 1,
  show_query_node_labels = 0
)

# Define a cache directory
cache_dir <- tempdir()

# Define a function to read the file with caching
read_tsv_with_cache <- memoise::memoise(function(file_url,header = FALSE) {
  # Generate a unique cache filename based on the file URL
  cache_filename <- fs::path_join(c(cache_dir, paste0(digest::digest(file_url), ".rds")))
  
  # Check if the cached file exists
  if (file.exists(cache_filename)) {
    # If cached file exists, load and return the cached data
    cached_data <- readRDS(cache_filename)
    return(cached_data)
  } else {
    # If cached file does not exist, download and cache the data
    data <- read.delim(url(file_url), sep = "\t", header = header)
    saveRDS(data, cache_filename)
    return(data)
  }
})

# read data from stringDB api
response <- httr::GET(paste(address, "/api/tsv/network", sep = ""), query = networkParams)
res <- read_tsv_with_cache(response$url, header = TRUE)

# read data from species code list
# Define the URLs for the files
kegg_species_url <- 'https://rest.kegg.jp/list/organism'
stringdb_species_url <- "https://stringdb-static.org/download/species.v11.5.txt"

# Read the files with caching
kegg_species <- read_tsv_with_cache(kegg_species_url)
stringdb_species <- read_tsv_with_cache(stringdb_species_url,header = TRUE)
