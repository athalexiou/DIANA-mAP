## Check, and install if needed, the packages required throughout the pipeline
packages <- c("parallel","optparse","RColorBrewer")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos='http://cran.us.r-project.org')  
}

## Argument Input
library("optparse", quietly = TRUE, warn.conflicts = FALSE)

option_list <- list(
  make_option(c("-d", "--difexpression"),
              help="Add this option to apply Differential Expression \n\t\tanalysis to the results of the samples run.\n\t\tProvide the path to the condition table for the \n\t\tdifferential expression analysis.\n\t\tIt should contain two columns with their headers,\n\t\tthe first being the sample column containing all \n\t\tthe sample names and the second should be the \n\t\tcondition column, containing a condition for each sample. \n\t\tNOTICE: The different conditions must be exactly two in number.",
              metavar="Path")
)
pargv = parse_args(OptionParser(usage = "\n%prog [options] input_directory output_directory number_of_cores pipeline_arguments", option_list = option_list, description = "\nRun multiple pipelines in parallel using forking with 'mclapply' and possibly apply Differential Expression analysis to their results. The 'number_of_cores' specifies how many parallel pipelines will be executed simultaneously at each time.\n\nThe pipeline arguments will be passed to every pipeline run with this script and can be in the form of: %prog \"input_directory\" \"output_directory\" \"number_of_cores\" \"-c 'value' -m 12 -d 'value'\" etc. If no parameters are to be passed to the pipeline, please use empty quotes to nominate that (eg. \"\"). For more information on the parameters available please use the help flag with the pipeline script.\n\nWARNING: If the 'thread' parameter is provided for each pipeline and its value is above 1, the total number or threads that will be required for this script to run is: (number_of_cores * threads)! Please make sure your system has enough available resources. Requesting for more resources than available can lead to multiple fatal errors and is not recommended."), positional_arguments = 4)

input_directory = pargv$args[[1]]
output_directory = pargv$args[[2]]
core_num = pargv$args[[3]]
pipeline_args = pargv$args[[4]]

if (substr(input_directory,nchar(input_directory),nchar(input_directory)) == "/") {
  input_directory = substr(input_directory, 1, nchar(input_directory)-1)
}

for (arg in 1:length(pargv$options)) {
  if (names(pargv$options[arg]) != "help") {
    assign(names(pargv$options)[arg],pargv$options[[arg]])
  }
}

library("parallel", quietly = TRUE, warn.conflicts = FALSE)

cat("ANALYSIS START TIME: ")
system("date")

## Gather the fastq files from the input folder
files = list.files(path = input_directory, pattern="*.f*q*", full.names = T, recursive = FALSE)

## Create the "_Failed_Runs" directory
system(paste0("mkdir ",output_directory,"_Failed_Runs/"))

## Run multiple instances of the miRNA_pipeline in parallel
cat(paste0(">>>> INITIALIZING PARALLEL PIPELINES FOR ",length(files)," SAMPLES, PLEASE WAIT... \n"))

pipeline_results = mclapply(files, mc.cores = core_num, function(infile) {
  input_base_name = strsplit(basename(infile),"\\.")[[1]][1]

  pipeline = system(paste0("Rscript diana_map.r ",pipeline_args," ",infile," ",output_directory, " > ", output_directory, input_base_name,"_output_log.txt 2>&1"))
  if (pipeline == 0) {
    cat(paste0(input_base_name,"> Finished successfully!\n\n"))
    system(paste0("mv ",output_directory, input_base_name,"_output_log.txt ",output_directory,input_base_name,"_*/"))
  }else {
    cat(paste0(input_base_name,"> Something went wrong; Please check the output log file.\n\n"))

    ## Move the output_log file in the sample folder and move the folder to the "_Failed_Runs" directory
    if (length(file.exists(Sys.glob(file.path(paste0(output_directory,input_base_name,"_*[0-9]/"))))) < 1) {
      system(paste0("mkdir ",output_directory, input_base_name,"_/"))
    }
    system(paste0("mv ",output_directory, input_base_name,"_output_log.txt ",output_directory,input_base_name,"_*/"))
    system(paste0("mv ",output_directory, input_base_name,"_*/ ",output_directory,"_Failed_Runs/"))
  }
})
cat("ANALYSIS END TIME: ")
system("date")

## Differential Expression Analysis
if (exists("difexpression")) {
  cat("DE ANALYSIS START TIME: ")
  system("date")
  cat(">>>> INITIATING DIFFERENTIAL EXPRESSION ANALYSIS... \n")
  bioconductor_packages <- c("DESeq2","pcaExplorer")
  if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(bioconductor_packages)
  }
  suppressMessages(library("DESeq2", quietly = TRUE, warn.conflicts = FALSE))
  suppressMessages(library("pcaExplorer", quietly = TRUE, warn.conflicts = FALSE))
  
  de_dir = paste0(output_directory,"DE_Analysis",format(Sys.time(), "_%Y%m%d_%H%M%S"),"/")
  dir.create(de_dir, showWarnings = FALSE)
  setwd(de_dir)
  
  ## Read condition table
  colData = read.csv(file = difexpression, stringsAsFactors = F)
  colData[[colnames(colData)[2]]] = factor(colData[[colnames(colData)[2]]])
  if (nlevels(colData[[colnames(colData)[2]]]) != 2) {
    stop(paste0("The conditions provided for the Differential Expression of the samples are not 2.\nPlease change the input table accordingly in order to provide exactly 2 different conditions.\n\nConditions given: ",toString(levels(factor(colData[[colnames(colData)[2]]]))),"\n"))
  }
  
  ## Populate raw_counts_table for DE Analysis
  for(sample in colData[,1]) {
    sample_counts_file = Sys.glob(file.path(output_directory,paste0(sample,"*"),paste0(sample,"*_Counts.txt")))
    if (length(sample_counts_file) > 0 && file.exists(sample_counts_file)) {
      sample_counts = read.table(sample_counts_file, sep="\t", col.names = c("miRNA_ID",paste0(sample,"_raw_counts"),as.character(sample),paste0(sample,"_log2")), header = FALSE, skip = 1 )
      
      if (!exists("master_table")){
        master_table = sample_counts[,c(1,2)]
      }else {
        master_table = merge(master_table,sample_counts[,c(1,2)],by = "miRNA_ID")
      }
    }
  }
  
  write.csv(master_table, file = "raw_counts_table.csv",row.names=FALSE)
  
  ## Read in the counts table
  countFilePath = "raw_counts_table.csv"
  countData = read.table(file = countFilePath, header = TRUE, sep = ",", row.names = 1)
  countData = round(countData,0)
  colnames(countData) = sub("_raw_counts","",colnames(countData))
  
  group_name = colnames(colData)[2]
  group1 = levels(factor(colData[[colnames(colData)[2]]]))[1]
  group2 = levels(factor(colData[[colnames(colData)[2]]]))[2]
  
  ## Import the counts table and run DESeq2 for the DE Analysis
  dataset = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = as.formula(paste0("~",group_name)))
  
  dds = DESeq(dataset)
  
  result = results(dds, contrast=c(group_name,group1,group2))
  result = result[complete.cases(result),] # Remove any rows with NA
  
  ## Write the full results
  write.csv(as.data.frame(result),file = "DE_full_results.csv")
  writeLines(c("\t\t\t\tDIFFERENTIAL EXPRESSION REPORT\n\n","\t\t\tSUMMARY:"),"DE_report.txt")
  capture.output(summary(result),file = "DE_report.txt", append = TRUE)
  
  ## DE Plots
  # PCA
  cat(">> PLOTTING PCA... \n")
  rld = rlogTransformation(dds, blind = TRUE)
  pcaobj = prcomp(t(assay(rld)))
  pdf("DE_PCA.pdf")
  print(pcaplot(rld, intgroup = group_name,title = paste0("Principal Component Analysis on ",group_name)))
  print(pcascree(pcaobj,type = "pev",title = "Proportion of Explained Variance per Principal Component"))
  dev.off()
  
  # MA
  cat(">> PLOTTING MA-PLOT... \n")
  png(filename="DE_MAplot.png")
  plotMA(result, main=paste0("MA-plot on ",group_name,": ",group1," vs. ",group2), ylim=c(-5,5))
  dev.off()
  
  ## Extract results for the top n up-regulated and the top n down-regulated miRNAs, sorted by adjusted p-value
  cat(">> EXTRACTING TOP RESULTS... \n")
  n = 100
  resOrdered = result[order(result$padj),]
  topResults = rbind(resOrdered[resOrdered[,"log2FoldChange"] > 0, ][1:n,], resOrdered[resOrdered[,"log2FoldChange"] < 0, ][n:1,])
  write.csv(as.data.frame(topResults),file = paste0("DE_top",n,"_by_padj_results.csv"))

  ## Report top and bottom 5 miRNAs
  de_report = readLines("DE_report.txt")
  writeLines(c(de_report," ","\t\t\tTOP AND BOTTOM 5 miRNAs BASED ON ADJUSTED P-VALUE:"),"DE_report.txt")
  capture.output(topResults[c(1:5,(2*n-4):(2*n)), c("baseMean","log2FoldChange","padj")],file = "DE_report.txt", append = TRUE)

  ## Plot count comparisons per miRNA from the top results
  cat(">> PLOTTING COUNT COMPARISON PLOTS PER miRNA... \n")
  pdf(paste0("top",n,"_miRNA_count_comparison_plots.pdf"))
  for (gn in row.names(topResults)) {
    plotCounts(dds, gene = gn, intgroup = group_name, pch = 19)
  }
  dev.off()
  
  ## DE Plot Heatmaps
  # Plot Top 100 miRNAs heatmap and cluster by similarity
  cat(">> PLOTTING HEATMAPS... \n")
  suppressMessages(library("RColorBrewer", quietly = TRUE, warn.conflicts = FALSE))
  hmcol = brewer.pal(11,"RdBu")
  nCounts = counts(dds, normalized = TRUE)
  png(filename=paste0("DE_top",n,"_by_padj_heatmap.png"), width = 8, height = 8, units = 'in', res = 600)
  heatmap(as.matrix(nCounts[row.names(topResults),]), Rowv = NA, col = hmcol, mar = c(10,2))
  dev.off()
  
  # Plot Top and Bottom 25 miRNAs heatmap and cluster by similarity
  m = 25
  png(filename=paste0("DE_top",m,"_by_padj_heatmap.png"), width = 6, height = 6, units = 'in', res = 600)
  heatmap(as.matrix(nCounts[row.names(topResults)[c(1:m,(n-m+1):n)], ]), Rowv = NA, col = hmcol, mar = c(10,2))
  dev.off()
  
  cat(">>>> DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED SUCCESSFULLY! \n")
  cat("DE ANALYSIS END TIME: ")
  system("date")
}