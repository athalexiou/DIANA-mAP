#!/usr/bin/R
## Check, and install if needed, the packages required throughout the pipeline
packages <- c("ggplot2", "tictoc","tcR","optparse")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos='http://cran.us.r-project.org')  
}

## Import the default configuration file for the default values
source("diana_map_config.r")

## Argument Input
library("optparse", quietly = TRUE, warn.conflicts = FALSE)

option_list <- list(
  make_option(c("-a", "--adapter"),
              help="The adapter sequence used for the dataset preparation. \n\t\tLeave empty if the sequence is unknown. \n\t\tChange to 'Raw_Input' if there is no adapter sequence present \n\t\tin the dataset. (default = '')"),
  make_option(c("-e", "--exp_name"),
              help="A name for the experiment of the specific sample. \n\t\tText. (default = '')"),
  make_option(c("-p", "--preprocess"),
              help="Apply the pre-process analysis step. \n\t\tBoolean. (default = 'TRUE')",
              metavar="Boolean"),
  make_option(c("-y", "--adapter_type"),
              help="Adapter orientation information. \n\t\t'a' for 3' adapters or 'g' for 5' adapters (default = 'a')"),
  make_option(c("-q", "--trim_quality"), type="integer",
              help="The minimum base quality allowed. \n\t\tInteger values from 0. (default = 10)",
              metavar="number"),
  make_option(c("-r","--trim_max_err_rate"), type = "double",
              help = "The percentage of the error rate defining the allowed mismatches \n\t\tbetween the adapter sequence and the read during the trimming process. \n\t\tFloat values from 0 to 1. (default = 0.1)",
              metavar="number"),
  make_option(c("-v", "--trim_min_adapter_overlap"), type="numeric",
              help="The minimum overlap with adapter sequence required (in bases) \n\t\tto trim a sequence. \n\t\tInteger values from 0. (default = 3)",
              metavar="number"),
  make_option(c("-l", "--trim_min_len"), type="numeric",
              help="The minimum read length accepted after the trimming process. \n\t\tInteger values from 0. (default = 18)",
              metavar="number"),
  make_option(c("-b", "--adapter_library"), type="character",
              help="The position of the adapter library fasta file. \n\t\tFile path. (default = '')",
              metavar="Path"),
  make_option(c("-t", "--adapter_identity_threshold"), type="numeric",
              help="The minimum percentage of identity required for a match \n\t\tof the inferred adapter to an adapter sequence \n\t\tin the adapter library. \n\t\tInteger values from 0 to 100. (default = 90)",
              metavar="number"),
  make_option(c("-k", "--adapter_kmer_size"), type="integer",
              help="The length of the k-mers exracted from the adapter sequence \n\t\tfor the k-mer loop pre-processing stage. \n\t\tInteger values from 0. (default = 10)",
              metavar="number"),
  make_option(c("-m", "--max_multimaps"), type="integer",
              help="The maximum number of multimaps allowed for each read \n\t\tduring the mapping process, reads with more maps \n\t\twill be discarded. \n\t\tInteger values from 0. (default = 5)",
              metavar="number"),
  make_option(c("-j", "--threads"), type="integer",
              help="The number of threads to be used. (default = 1)",
              metavar="number"),
  make_option(c("-g", "--path_ref_genome_prefix"),
              help="The indexed reference genome directory plus the prefix. \n\t\tText. (default = '')",
              metavar="Path"),
  make_option(c("-i", "--large_index"),
              help="Set TRUE if the genome index ends in '.ebwtl',\n\t\tset to FALSE if it ends in '.ebwt'\n\t\tBoolean. (default = TRUE)",
              metavar="Boolean"),
  make_option(c("-s", "--species"),
              help="The three-letter species notation for the \n\t\tspecies used in the analysis. \n\t\tText. (default = 'hsa')",
              metavar="Species"),
  make_option(c("-n", "--path_to_hairpin"),
              help="The miRBase hairpin miRNAs file location. \n\t\tText. (default = '')",
              metavar="Path"),
  make_option(c("-u", "--path_to_mature"),
              help="The miRBase mature miRNAs file location. \n\t\tText. (default = '')",
              metavar="Path"),
  make_option(c("-c", "--user_config"),
              help="A configuration file provided by the user. \n\t\tPath. (default = '')",
              metavar="Path")
)

argv = parse_args(OptionParser(usage = "\n%prog [options] input_file output_directory", option_list = option_list, description = "\nAn automated miRNA analysis pipeline from pre-processing to quantification.",epilogue = "\nThe variable loading order is as follows: Command_Line_Arguments > User_Configuration_File_Arguments > Default_Configuration_File_Arguments\n"), positional_arguments = 2 )

file_in = argv$args[[1]]
out_dir = argv$args[[2]]

if (substr(out_dir,nchar(out_dir),nchar(out_dir)) != "/") {
  out_dir = paste0(out_dir,"/")
}

## If provided, import the user configuration file to overwrite the default values
if ("user_config" %in% names(argv$options)) {
  source(argv$options[["user_config"]])
}

## Finally, import any command line arguments provided to overwrite the above values
for (arg in 1:length(argv$options)) {
  if (names(argv$options[arg]) != "help") {
    assign(names(argv$options)[arg],argv$options[[arg]])
  }
}

input_file_name = strsplit(basename(file_in),"\\.")[[1]][1]

if (large_index == TRUE){
  large_index_val = " --large-index"
}else {
  large_index_val = ""
}

## Start timer
library(tictoc, quietly = TRUE, warn.conflicts = FALSE)
tictoc::tic(input_file_name)

## Create new experiment directory 
if (exp_name != "") {
  exp_name = paste0("_",exp_name)
}
cat(paste0(">>>> CREATING DIRECTORY... ",input_file_name,exp_name,"\n"))
main_dir = paste0(out_dir,input_file_name,exp_name,format(Sys.time(), "_%Y%m%d_%H%M%S"),"/")
dir.create(main_dir)

setwd(main_dir)
final_report_file = paste0(input_file_name,exp_name,"_Final_Report.txt")
writeLines(c(paste0("\t\tFINAL REPORT FOR ",input_file_name,exp_name),"\n"),final_report_file)

if (preprocess) {
  ## Quality Check
  cat(paste0(">>>> QUALITY CHECK... ",input_file_name,exp_name,"\n"))
  system(paste(path_fastqc, "fastqc ", file_in, " -t ", threads, " --outdir ", getwd(), sep=""))
  
  unzip(paste0(input_file_name,"_fastqc.zip"),paste0(input_file_name,"_fastqc/fastqc_data.txt"))
  file.copy(paste0(input_file_name,"_fastqc.zip"),paste0(input_file_name,"_fastqc/"))
  file.remove(paste0(input_file_name,"_fastqc.zip"))
  
  adapter_exists = TRUE
  
  ## If an adapter is not provided, use DNApi to infer it
  if (adapter == "") {
    cat(paste0(">>>> ADAPTER FINDING... ",input_file_name,exp_name,"\n"))
    
    ## DNApi Exaustive search
    adapter_dnapi = tryCatch({
      system(paste0("python ",path_dnapi, "dnapi.py --output-dir ",getwd()," --map-command '",path_aligner_executable," ",path_ref_genome_prefix,large_index_val," -p",threads," -v0 -k1 -S -f @in > @out' ",file_in),intern = TRUE)
    }, error = function(e) {
      stop("ADAPTER FINDING ERROR: The adapter finding software failed to launch.")
    })
    adapter_dnapi = strsplit(strsplit(adapter_dnapi[1],"=")[[1]][2],"/")
    if (is.na(adapter_dnapi[[1]][2]) == FALSE) {
      cat(paste0(">> ADAPTER WARNING: Due to low mapping rates, your data is considered ",adapter_dnapi[[1]][2],", please consider checking and/or replacing your dataset.\n"))
    }
    adapter = adapter_dnapi[[1]][1]
    
    ## Check FastQC adapter results for possible adapters
    if (grepl("RAW_INPUT",adapter)) {
      fastqc_adapter_content = read.delim(pipe(paste0("awk '/>>Adapter/{f=1;next} /END_MODULE/{f=0} f' ",paste0(input_file_name,"_fastqc/fastqc_data.txt"))),header = TRUE, comment.char = "",stringsAsFactors = FALSE)
      fastqc_adapter_max_values = lapply(fastqc_adapter_content[,-1], function(x) max(x, na.rm = T))
      
      potential_adapters = which(sapply(fastqc_adapter_max_values, function(z) isTRUE(z > 7)))
      if (is.na(adapter_dnapi[[1]][2]) == FALSE) {
        potential_adapters = which(sapply(fastqc_adapter_max_values, function(z) isTRUE(z > 2)))
      }
      
      if (length(potential_adapters) > 0) {
        cat(">> NOTICE: Unable to infer significant adapter, trying FastQC adapter content input instead.\n")
        if (length(potential_adapters) == 1) {
          potential_adapter = names(potential_adapters)
        }else {
          potential_adapter = names(which.max(fastqc_adapter_max_values[potential_adapters]))
        }
        
        # Compare the adapter presence on the first and second half of the reads to determine the orientation/type of the possible adapter
        first_half_sum = sum(fastqc_adapter_content[1:(nrow(fastqc_adapter_content)%/%2),potential_adapter])
        second_half_sum = sum(fastqc_adapter_content[((round(nrow(fastqc_adapter_content)/2))+1):nrow(fastqc_adapter_content),potential_adapter])
  
        if (first_half_sum >= second_half_sum) {
          adapter_type = "g"
        }else {
          adapter_type = "a"
        }
        
        if (potential_adapter == "Illumina.Universal.Adapter") {
          adapter = "AGATCGGAAGAG"
        }else if (potential_adapter == "Illumina.Small.RNA.3..Adapter") {
          adapter = "ATGGAATTCTCG"
        }else if (potential_adapter == "Illumina.Small.RNA.5..Adapter") {
          adapter = "ATGGAATTCTCG"
        }else if (potential_adapter == "Nextera.Transposase.Sequence") {
          adapter = "CTGTCTCTTATA"
        }else if (potential_adapter == "SOLID.Small.RNA.Adapter") {
          adapter = "CGCCTTGGCCGT"
        }
      }
    }
    
    ## If DNApi outputs "RAW_INPUT" there is no adapter present in the dataset
    if (grepl("RAW_INPUT",adapter)) {
      adapter_exists = FALSE
      adapter_label = "Not Present"
    }else if (grepl("NA",adapter)) {
      stop("ADAPTER FINDING ERROR: The  adapter finding software did not run properly, please resolve the error above and try again.")
    }else {
      ## Cross reference against Adapter Library
      lib_res = strsplit(system(paste0("swan --key-value -id ",adapter_identity_threshold," -r ",adapter_library," -qs ",adapter),intern = TRUE)," ")
      
      if (length(lib_res) >= 1){
        if (length(lib_res) > 1){
          adp_identities = list()
          for (adp_hit in 1:length(lib_res)) {
            adp_identities[[adp_hit]] = as.integer(strsplit(lib_res[[adp_hit]][8],"=")[[1]][2])
          }
          max_iden = which.max(adp_identities)
          if (length(which(adp_identities==adp_identities[[max_iden]])) == 1) {
            top_adp_hit = max_iden
          }else {
            top_adp_hit = 0
            cat(paste0(">> NOTICE: Multiple top matches found in Adapter Library. As there is no way of knowing the exact one used for the dataset, proceeding with inferred adapter: ",adapter,"\n"))
            adapter_label = "Inferred"
          }
        }else {
          top_adp_hit = 1
        }
        if (top_adp_hit != 0) {
          adapter_iden = as.integer(strsplit(lib_res[[top_adp_hit]][8],"=")[[1]][2])
          adapter_label = gsub(".*\\{(.*)\\}.*", "\\1", lib_res[[top_adp_hit]][9])
          adapter_seq = gsub(".*\\{(.*)\\}.*", "\\1", lib_res[[top_adp_hit]][4])
          cat(paste0(">> NOTICE: Adapter match FOUND in Adapter Library! Label: ",adapter_label," - Seq: ",adapter_seq," - Identity: ",adapter_iden,"\n"))
          adapter = adapter_seq
        }
      }else {
        cat(paste0(">> NOTICE: Adapter match NOT FOUND in Adapter Library. Proceeding with inferred adapter: ",adapter,"\n"))
        adapter_label = "Inferred"
      }
    }
  }else {
    adapter_label = "Provided"
  }
  
  
  ## Quality Trimming and adapter removal
  cat(paste0(">>>> QUALITY TRIMMING AND ADAPTER REMOVAL... ",input_file_name,exp_name,"\n"))
  
  if (!(identical(adapter_type,"a") || identical(adapter_type,"g"))) {
    stop("ADAPTER TYPE ERROR: Adapter_type not 'a' or 'g'. Please provide 'a' for 3' adapters or 'g' for 5' adapters. (Default: 'a')")
  }
  
  
  if (adapter_exists) {
    untrimmed_output = paste0(" --untrimmed-output reads_without_adapter.fq")
    trim_adapter_input = paste0(" -",adapter_type," ",adapter," -O ",trim_min_adapter_overlap," -e ",trim_max_err_rate,untrimmed_output)
  }else {
    cat(">> NOTICE: Your data seem to be already processed, adapter removal will be omitted.\n")
    trim_adapter_input = ""
  }
  system(paste0(path_cutadapt, "cutadapt -q ",trim_quality,trim_adapter_input," -m ",trim_min_len," -M ",trim_max_len," -o ",input_file_name,"_trimmed.fq ",file_in," > ",input_file_name,"_trimming_report.txt"))
  
  
  ## Adapter Cleansing Loops
  reads_cleansed = 0
  total_untrimmed_reads = 0
  trim_report_file = readLines(paste0(input_file_name,"_trimming_report.txt"))
  
  if (adapter_exists) {
    ## Check "reads with adapter" for low percentage, warn and attempt the complementary adapter_type removal
    adapter_content = strsplit(grep("Reads with adapters:",trim_report_file,value = TRUE)," ")[[1]]
    adapter_content_percent = as.double(gsub("[()%]","",adapter_content[length(adapter_content)]))
    if (adapter_content_percent < 70) {
      
      original_adapter_type = adapter_type
      if (identical(adapter_type,"a")) {
        adapter_type = "g"
        adp_side = "5prime"
      }else {
        adapter_type = "a"
        adp_side = "3prime"
      }
      cat(paste(">> ADAPTER WARNING: The percentage of reads containing the given/inferred adapter is too low!",paste0("Adapter: ",adapter),grep("Reads with adapters:",trim_report_file,value = TRUE),paste0("Adapter removal from ",adp_side," side of the reads will be attempted..."),"\n", sep = "\n"))
    
      system(paste0(path_cutadapt, "cutadapt -q ",trim_quality," -",adapter_type," ",adapter," -O ",trim_min_adapter_overlap," -e ",trim_max_err_rate," --untrimmed-output reads_without_",adp_side,"_adapter.fq -m ",trim_min_len," -M ",trim_max_len," -o ",input_file_name,"_",adp_side,"_trimmed.fq ",file_in," > ",input_file_name,"_",adp_side,"_trimming_report.txt"))
      
      new_trim_report_file = readLines(paste0(input_file_name,"_",adp_side,"_trimming_report.txt"))
      new_adapter_content = strsplit(grep("Reads with adapters:",new_trim_report_file,value = TRUE)," ")[[1]]
      new_adapter_content_percent = as.double(gsub("[()%]","",new_adapter_content[length(new_adapter_content)]))
      
      if (new_adapter_content_percent <= adapter_content_percent) {
        # continue with original adapter type
        adapter_type = original_adapter_type
        
        cat(paste0(">> WARNING: Adapter removal from ",adp_side," side provided even lower percentage of reads containing the given/inferred adapter. Please consider checking the adapter given or the dataset.\n"))
        
      }else if (new_adapter_content_percent > adapter_content_percent) {
        # replace all cutadapt default output files with the new ones created above and continue with this adapter type
        file.rename(paste0("reads_without_",adp_side,"_adapter.fq"),"reads_without_adapter.fq")
        file.rename(paste0(input_file_name,"_",adp_side,"_trimmed.fq"),paste0(input_file_name,"_trimmed.fq"))
        file.rename(paste0(input_file_name,"_",adp_side,"_trimming_report.txt"),paste0(input_file_name,"_trimming_report.txt"))
        trim_report_file = new_trim_report_file
        
        cat(paste(paste0(">> INFO: Adapter removal from ",adp_side," side provided better percentage of reads containing the given/inferred adapter."),paste0("Adapter: ",adapter),grep("Reads with adapters:",trim_report_file,value = TRUE),paste0("The analysis will continue using the ",adp_side," side for adapter detection/removal."),"\n", sep = "\n"))
      }
    
      if (adapter_label != "Inferred" && new_adapter_content_percent < 70 && exists("adapter_dnapi")) {
        cat(paste0(">> WARNING: The percentage of reads containing an adapter from the known adapter library on both 3' and 5' is lower than 70%, the pipeline will fall back to the inferred adapter to continue the analysis as it is more likely to provide more accurate results. Please consider checking the adapter given or the dataset.\n"))
        
        # continue with original adapter type
        adapter_type = original_adapter_type
        
        adapter_label = "Inferred"
        adapter = adapter_dnapi[[1]][1]
        trim_adapter_input = paste0(" -",adapter_type," ",adapter," -O ",trim_min_adapter_overlap," -e ",trim_max_err_rate,untrimmed_output)
        file.remove("reads_without_adapter.fq")
        file.remove(paste0(input_file_name,"_trimmed.fq"))
        file.remove(paste0(input_file_name,"_trimming_report.txt"))
        system(paste0(path_cutadapt, "cutadapt -q ",trim_quality,trim_adapter_input," -m ",trim_min_len," -M ",trim_max_len," -o ",input_file_name,"_trimmed.fq ",file_in," > ",input_file_name,"_trimming_report.txt"))
        trim_report_file = readLines(paste0(input_file_name,"_trimming_report.txt"))
      }
      
    }
    
    ## Kmer Loop on reads without adapter
    suppressPackageStartupMessages(library(tcR, quietly = TRUE, warn.conflicts = FALSE))
    
    loop_dir = paste0(main_dir,input_file_name,"_Adapter_Loops/")
    dir.create(loop_dir)
    setwd(loop_dir)
    
    adapter_kmers = get.kmers(adapter,.k = adapter_kmer_size)[[1]]
    loop_count = 1
    
    ## Loop report
    writeLines(c(paste0("+++++++++++++++++++ K-mer: ",adapter_kmer_size," +++++++++++++++")," "),"Adapter_Loops_Report.txt") 
    
    for (i in grep("Total reads processed:",trim_report_file):grep("Reads written \\(passing filters\\):",trim_report_file)) { 
      loop_results = readLines("Adapter_Loops_Report.txt") 
      writeLines(c(loop_results,trim_report_file[i]),"Adapter_Loops_Report.txt") 
    } 
    ## /Loop report
    
    for (kmer in adapter_kmers) {
      print(paste0("********* K-mer = ",kmer))
      loop_id = paste0("loop",loop_count)
      if (loop_count <=1) {
        loop_file_in = paste0(main_dir,"reads_without_adapter.fq")
      }else {
        loop_file_in = paste0("reads_without_adapter_loop",(loop_count-1),".fq")
      }
      untrimmed_output = paste0(" --untrimmed-output reads_without_adapter_",loop_id,".fq")
      trim_adapter_input = paste0(" -",adapter_type," ",kmer," -O ",trim_min_adapter_overlap," -e ",trim_max_err_rate,untrimmed_output)
      system(paste0(path_cutadapt, "cutadapt -q ",trim_quality,trim_adapter_input," -m ",trim_min_len," -M ",trim_max_len," -o ",input_file_name,"_",loop_id,"_trimmed",".fq ",loop_file_in," > ",input_file_name,"_",loop_id,"_trimming_report.txt"))
      
      ## Loop report
      loop_trim_report_file = readLines(paste0(input_file_name,"_",loop_id,"_trimming_report.txt"))
      
      loop_results = readLines("Adapter_Loops_Report.txt")
      writeLines(c(loop_results," ",paste0("=========================> ",loop_id, " - K-mer: ",kmer)),"Adapter_Loops_Report.txt")
      
      if (length(grep("No reads processed!",loop_trim_report_file)) >= 1) {
        loop_results = readLines("Adapter_Loops_Report.txt")
        writeLines(c(loop_results,"No reads left to process."),"Adapter_Loops_Report.txt")
        break
      }
      
      for (i in grep("Total reads processed:",loop_trim_report_file):grep("Reads written \\(passing filters\\):",loop_trim_report_file)) {
        loop_results = readLines("Adapter_Loops_Report.txt")
        writeLines(c(loop_results,loop_trim_report_file[i]),"Adapter_Loops_Report.txt")
      }
      
      a = strsplit(grep("Reads with adapters:",loop_trim_report_file,value = TRUE)," ")[[1]]
      a = as.integer(gsub(",","",a[length(a)-1]))
      b = strsplit(grep("Reads that were too short:",loop_trim_report_file,value = TRUE)," ")[[1]]
      b = as.integer(gsub(",","",b[length(b)-1]))
      
      reads_cleansed = reads_cleansed + (a - b)
      
      if (loop_count <=1) {
        tmp = strsplit(grep("Total reads processed:",loop_trim_report_file,value = TRUE)," ")[[1]]
        tmp = as.integer(gsub(",","",tmp[length(tmp)]))
        total_untrimmed_reads = tmp
      }
      ## /Loop report
      
      loop_count = loop_count + 1
    }
    
    ## Loop report
    loop_results = readLines("Adapter_Loops_Report.txt")
    writeLines(c(loop_results," ",paste0("Total Reads Cleansed: ",reads_cleansed," - ",round((reads_cleansed/total_untrimmed_reads)*100,digits = 2),"%")),"Adapter_Loops_Report.txt")
    ## /Loop report
    
    ## Merge cleansed reads from kmer loops to initially cleansed output
    system(paste0("cat ",input_file_name,"_loop*_trimmed.fq >> ",main_dir,input_file_name,"_trimmed.fq"))
    
  }

}

setwd(main_dir)

if (preprocess) {
  input_file_trimmed = paste0(input_file_name,"_trimmed.fq")
  input_file_name_trimmed = paste0(input_file_name,"_trimmed")
}else {
  input_file_trimmed = file_in
  input_file_name_trimmed = input_file_name
}

# Quality Check after Preprocess
cat(paste0(">>>> QUALITY CHECK AFTER PRE-PROCESS... ",input_file_name,exp_name,"\n"))
system(paste0(path_fastqc, "fastqc ", input_file_trimmed, " -t ", threads," --outdir ", getwd()))

unzip(paste0(input_file_name_trimmed,"_fastqc.zip"),paste0(input_file_name_trimmed,"_fastqc/fastqc_data.txt"))
file.copy(paste0(input_file_name_trimmed,"_fastqc.zip"),paste0(input_file_name_trimmed,"_fastqc/"))
file.remove(paste0(input_file_name_trimmed,"_fastqc.zip"))

## Genome mapping
cat(paste0(">>>> MAPPING... ",input_file_name,exp_name,"\n"))
system(paste0(path_mirdeep, "mapper.pl ",input_file_trimmed," -e -i -j -l ",trim_min_len," -o ",threads," -m -r ",max_multimaps," -p '",path_ref_genome_prefix,large_index_val,"' -q -s experiment_processed_reads.fa -t experiment_aligned.arf -v -h > mapping_log.txt 2>&1"))

if (length(grep("Mapping statistics",readLines("mapping_log.txt"))) < 1) {
  stop("MAPPING ERROR: Please refer to mapping_log.txt in the results directory for more information.")
}
mapping_table = tryCatch({
  read.table(pipe(paste0("awk '/Mapping statistics/{f=1;next} /seq:/{f=0} f' ","mapping_log.txt")),header = FALSE, comment.char = "",stringsAsFactors = FALSE)
}, error = function(e) {
  stop("MAPPING ERROR: Please refer to mapping_log.txt in the results directory for more information.")
})
# Check mapped reads percentage
if (mapping_table[2,5] < 0.5) {
  cat(">> MAPPING WARNING: Mapped percentage to genome is less than 50% , consider checking your dataset or the reference genome provided.\n")
}

## Quantification
cat(paste0(">>>> QUANTIFICATION... ",input_file_name,exp_name,"\n"))
system(paste0(path_mirdeep, "quantifier.pl -r experiment_processed_reads.fa -p ",path_to_hairpin," -m ",path_to_mature," -W -d -t ",species," -g ",max_mismatch_to_precursors," > quantification_log.txt 2>&1"))

if (length(grep("Mapping statistics",readLines("quantification_log.txt"))) < 1) {
  stop("QUANTIFICATION ERROR: Please refer to quantification_log.txt in the results directory for more information.")
}
quantification_table = tryCatch({
  read.table(pipe(paste0("awk '/Mapping statistics/{f=1;next} /seq:/{f=0} f' ","quantification_log.txt")),header = FALSE, comment.char = "",stringsAsFactors = FALSE)
}, error = function(e) {
  stop("QUANTIFICATION ERROR: Please refer to quantification_log.txt in the results directory for more information.")
})
# Check quantified reads percentage
if (quantification_table[2,5] < 0.5) {
  cat(">> QUANTIFICATION WARNING: Quantified percentage is less than 50% , consider checking your dataset or the reference data provided.\n")
}

counts_data = read.csv(normalizePath(dir(getwd(),"miRNAs_expressed_all_samples*")), sep="\t" )[ ,c(1,2,6)]
colnames(counts_data) = c("miRNA_ID","read_count","read_count_norm")

if (transform_miRNAs_to_IDs){
  ## Transform the miRNA names of the quantifier result into mirBase IDs
  cat(paste0(">>>> TRANSFORMING MIRNAS TO MIRBASE IDS... ",input_file_name,exp_name,"\n"))
  
  ## Create a matching table of miRNA names with mirBase IDs from the mature miRNA file.
  sys_command = paste("grep '>", species, "' ", path_to_mature, " | cut -d ' ' -f 1,2 | cut -d '>' -f 2",sep = "")
  con_to_file = pipe(sys_command)
  matching_table = read.csv(con_to_file,header = FALSE,sep = " ")
  
  ## Replace miRNA names with mirBase IDs using the matching table.
  for (row in 1:nrow(matching_table)) {
    counts_data$miRNA_ID = replace(as.character(counts_data$miRNA_ID), counts_data$miRNA_ID == matching_table[row,1], as.character(matching_table[row,2]))
  }
}
## Sort and aggregate counts on IDs
final_counts = aggregate(.~miRNA_ID,data = counts_data,FUN = sum)

## Add log2 of normalized counts
final_counts["log2_read_count_norm"] = log2(final_counts["read_count_norm"]+1)

write.table(final_counts,file = paste(getwd(),"/",input_file_name,exp_name,"_Counts.txt",sep=""),sep = "\t",row.names = FALSE,quote = FALSE)


## Prepare Final Report file
final_report = readLines(final_report_file)

if (preprocess) {
  if (total_untrimmed_reads > 0) {
    reads_without_adapters = paste0("Reads without adapter:\t\t\t",format(total_untrimmed_reads,big.mark = ",",scientific = FALSE))
    report_cleansed_reads = paste0("Reads cleansed with Kmer loops:\t\t",format(reads_cleansed,big.mark = ",",scientific = FALSE)," (",round((reads_cleansed/total_untrimmed_reads)*100,digits = 1),"%)")
  }else {
    reads_without_adapters = "No untrimmed reads detected."
    report_cleansed_reads = "No cleansed reads produced."
  }
  adapter_used = paste0(" - Adapter: ",adapter," - ",adapter_label)
  mappable_reads = as.integer(strsplit(system(paste0("wc -l ",input_file_name,"_trimmed.fq"),intern = TRUE)," ")[[1]][1])/4
  con = pipe(paste0(" sed -n -e '/WARNING:/,/beginning of the adapter sequence/p' ",paste0(input_file_name,"_trimming_report.txt")))
  adapter_warning = readLines(con)[1:4]
  close(con)
  if (!is.na(adapter_warning[1])) {
    warnings = adapter_warning
  }else {
    warnings = ""
  }
  reads_with_Ns = paste0("Reads removed because they contained 'N's: ",system(paste0("grep -v '@|+|#' ",input_file_name,"_trimmed.fq | grep 'N' | wc -l"),intern = TRUE))

  preprocess_report = paste(grep("Total reads processed:",trim_report_file,value = TRUE),paste0(grep("Reads with adapters:",trim_report_file,value = TRUE),adapter_used),paste0(grep("Reads that were too short:",trim_report_file,value = TRUE),paste0(" - Length cutoff: ",trim_min_len)),paste0(grep("Reads that were too long:",trim_report_file,value = TRUE),paste0(" - Length cutoff: ",trim_max_len)),grep("Reads written \\(passing filters\\):",trim_report_file,value = TRUE)," ",reads_without_adapters,report_cleansed_reads," ",paste0(grep("Quality-trimmed:",trim_report_file,value = TRUE)," - Quality cutoff: ",trim_quality),"\n",paste0("Mappable Reads:\t\t\t",format(mappable_reads,big.mark = ",",scientific = FALSE),"\n"),warnings, sep = "\n")  
}else {
  preprocess_report = "The Pre-processing step has been omitted by user choice."
  reads_with_Ns = ""
}

## Stop timer
tictoc::toc(log=TRUE)
## Write the Final Report file
writeLines(c(final_report,"\t\t\tGENERAL STATISTICS:\n","PREPROCESSING:\n",preprocess_report,"\nMAPPING:\n",toString(mapping_table[1,]),toString(mapping_table[2,]),"\nCollapsed Reads:",gsub("-m","surpassing max limit of multimaps",readLines("bowtie.log")),"\n\nQUANTIFICATION:\n",reads_with_Ns," ",toString(quantification_table[1,]),toString(quantification_table[2,]),"\nCollapsed Reads:",readLines(Sys.glob(file.path(getwd(),"expression_analyses/*/bowtie_reads.out"))),"\n\n",paste0("Time to complete ",tictoc::tic.log()[[1]][1])),final_report_file)
tictoc::tic.clearlog()

if (preprocess) {
  ## Plot Read Length Distributions Before and After the analysis (parse the read lengths even if they are intervals instead of integers)
  library("ggplot2", quietly = TRUE, warn.conflicts = FALSE)
  seq_length_distr_before = read.table(pipe(paste0("awk '/Sequence Length Distribution/{f=1;next} /END_MODULE/{f=0} f' ",paste0(input_file_name,"_fastqc/fastqc_data.txt"))),header = TRUE, comment.char = "",stringsAsFactors = FALSE)
  if (grepl("-",seq_length_distr_before[length(seq_length_distr_before$X.Length),1]) || grepl("-",seq_length_distr_before[1,1])) {
    
    if (grepl("-",seq_length_distr_before[length(seq_length_distr_before$X.Length),1]) && grepl("-",seq_length_distr_before[1,1])) {
      seq_length_distr_before_maxLength = as.numeric(strsplit(as.character(seq_length_distr_before[length(seq_length_distr_before$X.Length),1]),"-")[[1]][2])
      seq_length_distr_before_minLength = as.numeric(strsplit(as.character(seq_length_distr_before[1,1]),"-")[[1]][2])
    }else if (grepl("-",seq_length_distr_before[length(seq_length_distr_before$X.Length),1])) {
      seq_length_distr_before_maxLength = as.numeric(strsplit(as.character(seq_length_distr_before[length(seq_length_distr_before$X.Length),1]),"-")[[1]][2])
      seq_length_distr_before_minLength = seq_length_distr_before[1,1]
    }else {
      seq_length_distr_before_maxLength = seq_length_distr_before[length(seq_length_distr_before$X.Length),1]
      seq_length_distr_before_minLength = as.numeric(strsplit(as.character(seq_length_distr_before[1,1]),"-")[[1]][2])
    }
    
    seq_length_distr_before = rbind(seq_length_distr_before, c(seq_length_distr_before_maxLength+1, 0))
    seq_length_distr_before = rbind(c(seq_length_distr_before_minLength-1, 0), seq_length_distr_before)
    
    ggplot(seq_length_distr_before, aes(x=factor(seq_length_distr_before$X.Length,levels=seq_length_distr_before$X.Length),group = 1)) + geom_line(aes(y=seq_length_distr_before$Count)) + scale_y_continuous(labels = scales::comma) + scale_x_discrete(breaks = seq_length_distr_before$X.Length) + labs(title = "Read Length Distribution Before Analysis",x="Read Length(bp)",y="Read Count") + theme(text = element_text(size=16),axis.text.x = element_text(angle = 75, hjust = 1))
    ggsave(paste0(input_file_name,exp_name,"_Seq_Length_Distribution_BeforeAnalysis.png"))
    
  }else {
    seq_length_distr_before_maxLength = as.numeric(seq_length_distr_before[length(seq_length_distr_before$X.Length),1])
    seq_length_distr_before_minLength = as.numeric(seq_length_distr_before[1,1])
    
    seq_length_distr_before = rbind(seq_length_distr_before, c(seq_length_distr_before_maxLength+1, 0))
    seq_length_distr_before = rbind(c(seq_length_distr_before_minLength-1, 0), seq_length_distr_before)
    
    ggplot(seq_length_distr_before, aes(X.Length)) + geom_line(aes(y=Count)) + scale_y_continuous(labels = scales::comma) + scale_x_continuous(breaks = seq_length_distr_before$X.Length[1]:seq_length_distr_before$X.Length[length(seq_length_distr_before$X.Length)]) + labs(title = "Read Length Distribution Before Analysis",x="Read Length(bp)",y="Read Count") + theme(text = element_text(size=16),axis.text.x = element_text(angle = 75, hjust = 1))
    ggsave(paste0(input_file_name,exp_name,"_Seq_Length_Distribution_BeforeAnalysis.png"))
  }
  
  seq_length_distr_after = read.table(pipe(paste0("awk '/Sequence Length Distribution/{f=1;next} /END_MODULE/{f=0} f' ",paste0(input_file_name,"_trimmed_fastqc/fastqc_data.txt"))),header = TRUE, comment.char = "",stringsAsFactors = FALSE)
  if (grepl("-",seq_length_distr_after[length(seq_length_distr_after$X.Length),1])) {
    seq_length_distr_after_maxLength = as.numeric(strsplit(as.character(seq_length_distr_after[length(seq_length_distr_after$X.Length),1]),"-")[[1]][2])
    seq_length_distr_after = rbind(seq_length_distr_after, c(seq_length_distr_after_maxLength+1, 0))
    
    ggplot(seq_length_distr_after, aes(x=factor(seq_length_distr_after$X.Length,levels=seq_length_distr_after$X.Length),group = 1)) + geom_line(aes(y=seq_length_distr_after$Count)) + scale_y_continuous(labels = scales::comma) + scale_x_discrete(breaks = seq_length_distr_after$X.Length) + labs(title = "Read Length Distribution After Analysis",x="Read Length(bp)",y="Read Count") + theme(text = element_text(size=16),axis.text.x = element_text(angle = 75, hjust = 1))
    ggsave(paste0(input_file_name,exp_name,"_Seq_Length_Distribution_AfterAnalysis.png"))
    
  }else {
    seq_length_distr_after_maxLength = seq_length_distr_after[length(seq_length_distr_after$X.Length),1]
    seq_length_distr_after = rbind(seq_length_distr_after, c(seq_length_distr_after_maxLength+1, 0))
    
    ggplot(seq_length_distr_after, aes(X.Length)) + geom_line(aes(y=Count)) + scale_y_continuous(labels = scales::comma) + scale_x_continuous(breaks = seq_length_distr_after$X.Length[1]:seq_length_distr_after$X.Length[length(seq_length_distr_after$X.Length)]) + labs(title = "Read Length Distribution After Analysis",x="Read Length(bp)",y="Read Count") + theme(text = element_text(size=16),axis.text.x = element_text(angle = 75, hjust = 1))
    ggsave(paste0(input_file_name,exp_name,"_Seq_Length_Distribution_AfterAnalysis.png"))
  }
  
  ## Plot Pie Chart
  total_reads = strsplit(grep("Total reads processed:",trim_report_file,value = TRUE)," ")[[1]]
  total_reads = as.integer(gsub(",","",total_reads[length(total_reads)]))
  
  reads_too_short = strsplit(grep("Reads that were too short:",trim_report_file,value = TRUE)," ")[[1]]
  reads_too_short = as.integer(gsub(",","",reads_too_short[length(reads_too_short)-1]))
  
  reads_too_long = strsplit(grep("Reads that were too long:",trim_report_file,value = TRUE)," ")[[1]]
  reads_too_long = as.integer(gsub(",","",reads_too_long[length(reads_too_long)-1]))
  
  slices = c((total_untrimmed_reads-reads_cleansed), reads_cleansed, (mappable_reads-reads_cleansed), reads_too_short, reads_too_long)
  lbls = c("Without_Adapter", "Mappable (Cleansed)", "Mappable", "Too_Short", "Too_Long")
  colrs = c("Red3","Blue4","Royalblue1","Yellow2","Green2")
  pct = paste0(format(slices,big.mark = ",",scientific = FALSE), " (",round(slices/total_reads*100,digits = 2),"%)")
  
  ## Remove zero count categories if any
  zero_read_pos = which(slices == 0)
  if (length(zero_read_pos) > 0) {
    slices = slices[-zero_read_pos]
    lbls = lbls[-zero_read_pos]
    colrs = colrs[-zero_read_pos]
    pct = pct[-zero_read_pos]
  }
  
  pct[1] = paste0(pct[1],"\n\n")
  png(filename=paste0(input_file_name,exp_name,"_PieChart.png"), width = 6, height = 6, units = 'in', res = 600)
  par(mar=c(3,3,6,3))
  pie(slices,labels = pct, col=colrs,main=paste0("Total Reads Processed: ",format(total_reads,big.mark = ",",scientific = FALSE),"\nMappable Reads: ",paste0(format(mappable_reads,big.mark = ",",scientific = FALSE), " (",round(mappable_reads/total_reads*100,digits = 2),"%)"),"\n\n"),init.angle = 90,radius = 1, cex = 1.3, cex.main = 1.5) 
  legend("center", lbls, cex = 1.3, fill = colrs)
  dev.off()
  
  if (file.exists("Rplots.pdf")) {
    system("rm Rplots.pdf")
  }
}

## Folder Cleanup
intermed_results = paste0(main_dir,input_file_name,"_Intermediate_Results/")
dir.create(intermed_results)
qc_results = paste0(main_dir,input_file_name,"_QC/")
dir.create(qc_results)

system(paste0("mv mapper_logs/ expression_analyses/ *.f* *.log *.csv *.arf ",intermed_results))
if (preprocess) {
  system(paste0("mv *_report.txt ",intermed_results))
}
system(paste0("mv *_fastqc/ ",qc_results))

setwd(out_dir)
