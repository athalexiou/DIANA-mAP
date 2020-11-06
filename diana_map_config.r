#!/usr/bin/R

## General paths
path_fastqc <- "/path/to/FastQC/folder/" # path to fastqc executable folder
path_dnapi <- "/path/to/DNApi/folder/" # path to DNApi executable folder
path_cutadapt <- "/usr/bin/" # path to cutadapt executable folder or '/usr/bin/'
path_mirdeep <- "/path/to/mirdeep2/folder/bin/" # path to mirdeep2 executable folder
path_aligner_executable <- "/path/to/mirdeep2/folder/essentials/bowtie-1.1.1/bowtie" # path to the aligner executable in the essentials subfolder of mirdeep2 (including the executable), used in the DNApi adapter inference through mapping.
# the path to mirdeep-2/essentials folder needs to be on the path variable!

#######################################################################

## General info
exp_name = '' # Optional experiment name to be used in the created directory.
threads = 1 # The number of threads to be used by all the parallelizable processes of the analysis.

#######################################################################

## Specify analysis-specific variables and paths.

# Pre-processing (Quality Trimming and Adapter Removal)
preprocess = TRUE # Flag to perform the Pre-process step or not.
adapter = "" # Please provide the adapter used for your data within the quotes. Leave empty quotes if adapter is unknown or "RAW_INPUT" if your data are already processed. 
adapter_type = "a" # Please provide "a" for 3' adapters or "g" for 5' adapters.
trim_quality = 10 # The minimum quality of a base allowed, lower quality bases will be trimmed.
trim_max_err_rate = 0.1 # The allowed mismatch error rate between adapter and read sequence during a match.
trim_min_adapter_overlap = 3 # Minimum overlap between adapter and read sequence required to trim a sequence.
trim_min_len = 18 # The minimum allowed length of a read after trimming, reads with less bases will be discarded.
trim_max_len = 50 # The maximum allowed length of a read after trimming, reads with more bases will be discarded.
adapter_library = "/path/to/known_adapter_library.fa" # The path to the known adapter library fasta file.
adapter_identity_threshold = 90 # The minimum identity percetange threshold of an adapter with a known adapter in the above library.
adapter_kmer_size = 10 # The kmer size an adapter is split during the pre-processing loop process, please read the publication for more info.

# Mapping
max_multimaps = 5 # The maximum number of allowed multimaps for reads during the alignment process.
path_ref_genome_prefix = "/path/to/indexed/reference/genome/folder/ref_genome_prefix" # Path to the reference genome folder (including the genome prefix).
large_index = TRUE # Bowtie large-index flag option, for more info please read the bowtie manual online.

# Quantification
max_mismatch_to_precursors = 1 # The maximum allowed mismatched bases during the alignment of reads to known miRNA precursors.
species = "hsa" # The species studied in the experiment as required by miRBase.
path_to_hairpin = "/path/to/miRBase/hairpin/file/hairpin_v22.fa" # Path to the miRBase hairpin fasta file.
path_to_mature = "/path/to/miRBase/mature/file/mature_v22.fa" # Path to the miRBase mature fasta file.
transform_miRNAs_to_IDs = FALSE # Transform the miRNA names to their equivalent miRBase IDs for the quantification results.
