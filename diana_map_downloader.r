# Setup the absolute path to sratoolkit
sratoolkit_path = "/home/NGS/Progs/sratoolkit.2.9.6-1-centos_linux64/"


# Checking if needed packages are installed
packages <- c("HelpersMG","stringi","stringr","devtools")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos='http://cran.us.r-project.org')  
}

bioconductor_packages <- c("GEOquery","ENCODExplorer","BiocParallel")
if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
  BiocManager::install(bioconductor_packages)
}

cat("\nLoading required packages...\n\n")
if (!all(invisible(sapply(c(packages,bioconductor_packages), function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE)))))) {
  cat("\nCould not load packages succesfully.\n")
}

# Check and load docopt package to display documenation
if(!require("docopt", quietly = T)) {
  install_github("docopt/docopt.R")
  library("docopt", logical.return = TRUE, quietly = T)
}

# Setting the documentation of the script console call
doc = " 
DIANA-mAP script to automatically download biological data sets from various online databases. Tested using R version 3.6.2 (2019-12-12).
Current databases supported: NCBI GEO, NCBI SRA, ENCODE.

Note: Absolute path to sratoolkit needs to be set manually through line 2 in the script.

Usage: 
diana_map_downloader.r [-d <output_path>] [--geofilter <filter_regex>] [--encodefilter <filter_type>] [--srafilter <filter_type>] (-i <input>...)

Options: 
-d <output_path>, --directory <output_path>     Output saving path [default: ./]
--geofilter <filter_regex>                      Filter downloaded GEO files by regular expression ex. BED [default: NULL] (only works for R version > 3.5.0)
--encodefilter <filter_type>                    Filter downloaded ENCODE files ex. fastq [default: all]
--srafilter <filter_type>                       File type to be download ex.fastq or sam to call the corresponding 
sratoolkit dumper. Only fastq and sam files supported. [default: fastq]
-h --help                                       Show this screen

Arguments:
-i <input>...                                   One or more file accession numbers to download. Also, a file path may be provided including one accession number per line. You need to include the -i before each argument.
"

# Picking the arguments from the script call
opts = docopt(doc)

# Default null value to geofilter (was string)
if (stri_cmp_eq("NULL", opts$geofilter))
  opts$geofilter = NULL

# Create a history file keeping the status of the downloaded files 
if (!file.exists(paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""))) {
  
  write(paste("AccessionID", "\tStatus", "\tLink"), file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""),
        append = TRUE)
  
}

# Function that implements GEO download
geoDownload = function(i, opts) {
  
  cat(paste("\nDownloading ", i, "...\n", sep = ""))
  
  # Download the data
  
  # This works in R.3.5.0+ version only (the filter argument)
  #geoFileDownload = getGEOSuppFiles(i, makeDirectory = TRUE, baseDir = opts$directory,
  #                                  fetch_files = TRUE, filter_regex = opts$geofilter)
  
  # Checking if a file with the given accession number exists
  geoFile = tryCatch(
    {
      getGEO(i, getGPL = FALSE)
    },
    error=function(e){NULL}
  )
  
  if (is.null(geoFile)) {
    
    cat(paste("\nDownload was empty. No file with accession number", i, "was found.\n"))
    
    # Keep record of the non-found file
    write(paste(i,"\tnot found", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", i, sep = ""), 
          file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
    return()
    
  }
  
  # GSM case
  if (grepl("^GSM", i)) {
    
    # Get supplementary of GSM
    getGEOSuppFiles(i, makeDirectory = FALSE, baseDir = paste0(opts$directory, i), filter_regex = opts$geofilter)
    
    # If the given GSM is associated with SRA sequencing data, go to SRA
    if (geoFile@header$type == "SRA") {
      
      sra_i = strsplit(geoFile@header[["relation"]][grepl("SRA", geoFile@header[["relation"]])], "term=")[[1]][2]
      
      # Keep record of the found SRA file
      write(paste(i,"->", sra_i, "\tfound", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", i, sep = ""), 
            file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
      
      sraDownload(sra_i, opts)
      return()
      
    }
    
    # Keep record of the GSM found file
    write(paste(i,"\tfound", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", i, sep = ""), 
          file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
    
    # Create a folder to store the files, if it does not exist already
    if (!file.exists(i))
      dir.create(file.path(opts$directory, i))
    
    # Sample intensity data
    write.table(Table(geoFile), paste(opts$directory, i, "/", i, "_intensities.txt", sep = ""), 
                row.names = FALSE, col.names = TRUE, na = "", sep = "\t", quote = FALSE)
    
    # Descriptions for the sample data
    write.table(Columns(geoFile), paste(opts$directory, i, "/", i, "_intensity_descriptions.txt", sep = ""), 
                row.names = FALSE, col.names = TRUE, na = "", sep = "\t", quote = FALSE)
    
    # Experiment metadata
    fileMeta = plyr::ldply(Meta(geoFile), rbind)
    
    write.table(fileMeta, paste(opts$directory, i, "/", i, "_experiment_metadata.txt", sep = ""), 
                row.names = FALSE, col.names = FALSE, na = "", sep = "\t", quote = FALSE)
    
  }
  
  # GDS case
  else if (grepl("^GDS", i)) {
    
    # Keep record of the found GDS file
    write(paste(i,"\tfound", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", i, sep = ""), 
          file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
    
    # For each platform included in the given dataset
    for (gpl in 1:length(geoFile@header$platform)) {
      
      # Get GSM list annotated with the GPL
      gsm_list = geoFile@dataTable@columns$sample
      gpl_annot = geoFile@header$platform[gpl]
      
      # Download supplementary of GPL
      getGEOSuppFiles(gpl_annot, makeDirectory = TRUE, baseDir = paste0(opts$directory, i), filter_regex = opts$geofilter)
      
      # Download each GSM
      lapply(gsm_list[1:2], function(gsm) {
        
        # Create directory
        if (!file.exists(paste0(opts$directory, i, "/", gpl_annot))) {
          dir.create(file.path(paste0(opts$directory, i, "/", gpl_annot)), recursive = TRUE)
        }
        
        gsm_data = getGEO(gsm, getGPL = FALSE, destdir = paste0(opts$directory, i, "/", gpl_annot))
        # Download supplementary of GSM
        getGEOSuppFiles(gsm, makeDirectory = FALSE, baseDir =paste0(opts$directory, i, "/", gpl_annot), filter_regex = opts$geofilter)
        
        # If GSM is associated with SRA sequencing data, go to SRA 
        if (gsm_data@header$type == "SRA") {
          
          sra_i = strsplit(gsm_data@header[["relation"]][grepl("SRA", gsm_data@header[["relation"]])], "term=")[[1]][2] 
          
          # Keep record of the found SRA file
          write(paste(gsm_data@header$geo_accession,"->", sra_i, "\tfound", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_data@header$geo_accession, sep = ""), 
                file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
          
          sraDownload(sra_i, opts)
          
        }
        
        # Keep record of the found GSM file
        write(paste(gsm_data@header$geo_accession,"\tfound", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_data@header$geo_accession, sep = ""), 
              file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
        
      })
      
      # GPL experiment metadata
      write.table(geoFile@dataTable@columns, paste(opts$directory, i, "/", gpl_annot, "/", gpl_annot, "_metadata.tsv", sep = ""), sep='\t', row.names = FALSE )
      
    }
  }
  
  # GSE case
  else if (grepl("^GSE", i)) {
    
    cat(opts$geofilter)
    # Get supplementary of GSE
    getGEOSuppFiles(i, makeDirectory = TRUE, baseDir = opts$directory, filter_regex = opts$geofilter)
    
    # For each platform included in the given series
    for (gpl in 1:length(geoFile)) {
      
      # Get GSM list annotated with the GPL
      gsm_list = geoFile[[gpl]]$geo_accession
      gpl_annot = geoFile[[gpl]]@annotation
      
      # Download supplementary of GPL
      getGEOSuppFiles(gpl_annot, makeDirectory = TRUE, baseDir = paste0(opts$directory, i), filter_regex = opts$geofilter)
      
      # Download each GSM
      lapply(gsm_list[1:2], function(gsm) {
        
        # Create directory
        if (!file.exists(paste0(opts$directory, i, "/", gpl_annot))) {
          dir.create(file.path(paste0(opts$directory, i, "/", gpl_annot)), recursive = TRUE)
        }
        
        gsm_data = getGEO(gsm, getGPL = FALSE, destdir = paste0(opts$directory, i, "/", gpl_annot))
        # Download supplementary of GSM
        getGEOSuppFiles(gsm, makeDirectory = FALSE, baseDir =paste0(opts$directory, i, "/", gpl_annot), filter_regex = opts$geofilter)
        
        # If GSM is associated with SRA sequencing data, go to SRA 
        if (gsm_data@header$type == "SRA") {
          
          sra_i = strsplit(gsm_data@header[["relation"]][grepl("SRA", gsm_data@header[["relation"]])], "term=")[[1]][2]
          
          # Keep record of the found SRA file
          write(paste(gsm_data@header$geo_accession,"->", sra_i, "\tfound", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_data@header$geo_accession, sep = ""), 
                file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
          
          sraDownload(sra_i, opts)
          
        }
        
        # Keep record of the found GSM file
        write(paste(gsm_data@header$geo_accession,"\tfound", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_data@header$geo_accession, sep = ""), 
              file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
        
      })
      
      # GPL Experiment metadata
      write.table(geoFile[[gpl]]@phenoData@data, paste(opts$directory, i, "/", gpl_annot, "/", gpl_annot, "_metadata.tsv", sep = ""), sep='\t', row.names = FALSE )
      
    }
  }
  
  # GPL case
  else {
    
    # Download supplementary of GPL
    getGEOSuppFiles(gpl_annot, makeDirectory = TRUE, baseDir = paste0(opts$directory, i), filter_regex = opts$geofilter)
    
    # Create directory
    if (!file.exists(paste0(opts$directory, i))) {
      dir.create(file.path(paste0(opts$directory, i)), recursive = TRUE)
    }
    
    # Get GSM list annotated with the GPL
    gsm_list = geoFile@header$sample_id
    
    # Download each GSM
    lapply(gsm_list[1:2], function(gsm) {
      
      gsm_data = getGEO(gsm, getGPL = FALSE, destdir = paste0(opts$directory, i))
      # Download supplementary of GSM
      getGEOSuppFiles(gsm, makeDirectory = FALSE, baseDir =paste0(opts$directory, i), filter_regex = opts$geofilter)
      
      # If GSM is associated with SRA sequencing data, go to SRA 
      if (gsm_data@header$type == "SRA") {
        
        sra_i = strsplit(gsm_data@header[["relation"]][grepl("SRA", gsm_data@header[["relation"]])], "term=")[[1]][2]
        
        # Keep record of the found SRA file
        write(paste(gsm_data@header$geo_accession,"->", sra_i, "\tfound", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_data@header$geo_accession, sep = ""), 
              file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
        
        sraDownload(sra_i, opts)
        
      }
      
      # Keep record of the found GSM file
      write(paste(gsm_data@header$geo_accession,"\tfound", "\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_data@header$geo_accession, sep = ""), 
            file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
      
    })
    
    # GPL experiment metadata
    write.table(data.frame(unlist(geoFile@header)), paste(opts$directory, i, "/", i, "_metadata.tsv", sep = ""), sep='\t', row.names = TRUE, col.names = F)
    
  }
  
  cat("\n")
  return()
  
}

# Function that implements SRA download
sraDownload = function(i, opts) {
  
  # If the accession number matches a fastq file, call fastq-dump
  if (opts$`--srafilter` == "fastq") {
    
    cat(paste("\nDownloading ", i, "...\n", sep = ""))
    
    # Download the file and convert it (synchronised) to fastq
    sraDownloadStatus = system(paste(sratoolkit_path,"bin/fasterq-dump --split-files --outdir ", opts$directory, i, " ", i, sep = ""))
    
    # If a file was not found, store it to record and continue to the next one
    if (sraDownloadStatus != 0) {
      
      cat(paste("Failed to download ", i, "...\nCheck ", "https://www.ncbi.nlm.nih.gov/sra/", i, "\n", sep = ""))
      
      # Remove the folder created by fastqdump, as it does not contain valid files
      unlink(i, recursive = TRUE)
      # Save to the record file
      write(paste(i,"\tnot found", "\thttps://www.ncbi.nlm.nih.gov/sra/", i, sep = ""), 
            file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
      return()
      
    }
    
    cat(paste("Download complete.\nFile stored at ", opts$directory, i, "\n", sep = ""))
    
    # Keep record of the found SRA file
    write(paste(i,"\tfound", "\thttps://www.ncbi.nlm.nih.gov/sra/", i, sep = ""), 
          file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
    
    # Also get the metadata
    system(paste("wget -O ", opts$directory, i, "/", i, "_metadata.csv ",
                 " 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=", i, "'", sep = ""))
  }
  # If the accession number matches a SAM file, call sam-dump
  else if (opts$`--srafilter` == "sam") {
    
    cat(paste("\nDownloading ", i, "...\n", sep = ""))
    
    # Download the file and convert it (synchronised) to SAM
    sraDownloadStatus = system(paste(sratoolkit_path,"bin/sam-dump --output-file ", 
                                     opts$directory, i, "/", i, ".sam ", i,  sep = ""))
    
    # If a file was not found, store it to record and continue to the next one
    if (sraDownloadStatus != 0) {
      
      cat(paste("Failed to download ", i, "...\nCheck ", "https://www.ncbi.nlm.nih.gov/sra/", i, "\n", sep = ""))
      
      # Keep record of the found SRA file
      write(paste(i,"\tnot found", "\thttps://www.ncbi.nlm.nih.gov/sra/", i, sep = ""), 
            file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
      return()
      
    }
    
    # Create a folder to store the SAM file and the metadata
    if (!file.exists(i))
      dir.create(file.path(opts$directory, i))
    
    cat(paste("Download complete.\nFile stored at ", opts$directory, i, "\n", sep = ""))
    
    # Save to the record file
    write(paste(i,"\tfound", "\thttps://www.ncbi.nlm.nih.gov/sra/", i, sep = ""), 
          file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
    
    # Also get the metadata
    system(paste("wget -nv -O ", opts$directory, i, "/", i, "_metadata.csv ",
                 " 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=", i, "'", sep = ""))
  }
  
  cat("\n")
  return()
  
}

encodeDownload = function(i, opts) {
  
  cat(paste("\nDownloading ", i, "...\n", sep = ""))
  fuzzy_result = fuzzySearch(i, get_encode_df(), filterVector=c("accession","file_accession"))
  
  # If no file with the given accession number was found within encode_df, try the API
  if (length(fuzzy_result$file_accession) == 0){
    
    cat(paste("\nSearching through the ENCODE API...\n"))
    encFile = tryCatch(
      {
        system(paste("wget -nv --content-disposition -P ", opts$directory, i, 
                     " https://www.encodeproject.org/files/", i, "/@@download", sep = ""))
      },
      error=function(e){NULL}
    )
    
    # If file was found with the API, try to create a metadata file and keep record
    if (encFile == 0) {
      
      # Read the JSON given by the API (it includes many "thrash" information)
      encMetaFile = jsonlite::fromJSON(paste("https://www.encodeproject.org/files/", i, "/?format=json", sep = ""))
      # Flatten the levels created by the lists of lists so it can be written in tsv file
      encMetaFile = data.frame(unlist(encMetaFile))
      
      write.table(encMetaFile, paste(opts$directory, i, "/", i, "_metadata.tsv", sep = ""),
                  row.names = TRUE, col.names = FALSE, na = "", sep = "\t", quote = FALSE)
      write(paste(i,"\tfound", "\thttps://www.encodeproject.org/files/", i, "/", sep = ""), 
            file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
      
      return()
      
    }
    # If no file with the given accession number was found using the API, keep record
    # and continue to the next one
    else {    
      
      cat(paste("No results through the API\n"))
      write(paste(i,"\tnot found", "\thttps://www.encodeproject.org/files/", i, "/", sep = ""), file = "diana_mAP_file_download_status.tsv",
            append = TRUE)
      return()
      
    }
    
  }
  # If file was found with the API, keep record of its status
  else {
    
    write(paste(i,"\tfound", "\thttps://www.encodeproject.org/files/", i, "/", sep = ""), 
          file = paste(opts$directory, "diana_mAP_file_download_status.tsv", sep = ""), append = TRUE)
    
    if (!file.exists(i))
      dir.create(file.path(opts$directory, i))
    
    # Download the asked file
    downloadEncode(fuzzy_result, format = opts$`--encodefilter`, dir = paste(opts$directory, i, sep = ""), force = TRUE)
    # Download the metadata tied to the file found
    write.table(t(fuzzy_result), paste(opts$directory, i, "/", i, "_metadata.tsv", sep = ""),
                row.names = TRUE, col.names = TRUE, na = "", sep = "\t", quote = FALSE)
    
  }
  
  cat("\n")
  return()
  
}

# Function to call the corresponding db downloader
downloadFile = function(i, opts) {
  
  # If GEO accession number is found, use GEOquery    
  if (grepl("^GSM|GDS|GPL|GSE", i)) {
    
    geoDownload(i, opts)
    
  }
  # If ENCODE accession number is found, use ENCODExplorer
  else if (grepl("^ENC", i)) {
    
    encodeDownload(i, opts)
    
  }
  # If SRA accession number is found, use sratoolkit
  else if (grepl("^SR", i)) {
    
    sraDownload(i, opts)
    
  }
  else {
    
    cat(paste("\nNo file with accession number", i, "was found.\n"))
    
  }
  
}

# If accession numbers where given and not query,
# try to download them
if (!is.null(opts$i)) {
  
  for (i in opts$i) {
    
    # If input corresponds to a valid file, try to read accession numbers from there, one per line
    if(file_test("-f", i)) {
      
      con = file(i, "r")
      while(TRUE) {
        
        line = readLines(con, n = 1)
        if (length(line) == 0) {
          break
        }
        
        downloadFile(line, opts)
        
      }
      
      close(con)
      
    }
    else
      downloadFile(i, opts)
    
  }
  
} else if (!is.null(opts$query)) {
  
  #TODO: Implement the query case
  
}

closeAllConnections()