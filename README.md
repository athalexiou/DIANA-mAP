# DIANA-mAP
DIANA-mAP is an automated microRNA analysis pipeline with an emphasis on pre-processing.

The tool is free to use under MIT licence, please cite:
>Alexiou, A.; Zisis, D.; Kavakiotis, I.; Miliotis, M.; Koussounadis, A.; Karagkouni, D.; Hatzigeorgiou, A.G. DIANA-mAP: Analyzing miRNA from Raw NGS Data to Quantification. Genes 2021, 12, 46. (DOI: [10.3390/genes12010046](https://doi.org/10.3390/genes12010046))

## INSTALLATION

### DOCKER

DIANA-mAP is available as a Docker image in the public repository of Docker Hub under the name of _**athalexiou/diana_map**_. Inside the image all the scripts are available, along with the human GRCh38 reference genome, the v21 and v22 of the miRBase hairpin and mature known miRNAs and finally 4 test sample datasets allowing for a test run upon loading the image on a container. This is the recommended way of running this tool, as there is no need for any complex dependency installations.

Through the volume option of Docker (-v flag explained below), it is easy to provide the data to be analyzed by DIANA-mAP and the results to be accessible and/or stored locally outside of the container. Please do note that any results not transfered outside of the container will be lost upon exit.


1. **Follow the installation instructions on the docker website to install docker on your machine.**
	- If you are a user in an HPC or cluster and do not have root access, ask the IT administrator to install Docker for you and to follow [this guide](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user) in order to allow your user to run Docker commands without the need for root access or sudo.

2. **Start a container using the DIANA_mAP image from Docker Hub.**
    1. Provide a directory in your local machine to connect with the container in the flag -v of the command below. By doing that, you will be able to transfer data in or out of the container through the directory "/home/my-data/" inside the container.
         - `docker run -it --name diana_map_container -v /local/directory/to/connect/to/container/:/home/my-data/ -w="/home/DIANA_mAP/" --entrypoint "/home/init.sh" athalexiou/diana_map:v1.0 /bin/bash`
    2. The above command with start the container and position your terminal into the "DIANA_mAP" folder, where all the DIANA_mAP scripts are.
    3. To test the pipeline, the following command will perform an analysis on the 4 test samples provided in the image. Please provide the desired output directory and the number of cores to use in the following command:
        - `Rscript diana_map_multi.r -d /home/DIANA_mAP/example_DE_conditions_table.csv /home/test_data/ /path/to/output/directory/ number_of_cores ""`


### STANDALONE

DIANA_mAP can be downloaded from this repository and run standalone on unix-based systems. The following dependencies must be installed and their paths adjusted into the "diana_map_config.r" file. All the downloaded files of this repository must remain in the same directory. Due to the amount of dependencies and the required R libraries, it is highly recommended to run this tool through Docker as described above.

For Red Hat, CentOS, and Fedora systems substitute all the "apt-get" with "yum" on the following commands.

1. **Install R** (tested with v3.6)

2. **Install JAVA jre** (for FastQC)
    - `apt-get install openjdk-8-jre`
    - test: `java -version`

3. **Install FastQC** (tested with v.0.11.7)
    - Download the fastqc.zip from the fastqc repository
    - Unzip into folder of choice, get into created folder
    - `chmod 755 fastqc`	(provide execution rights to the fastqc executable)
    - `sudo ln -s /path/to/fastqc/fastqc /usr/local/bin/fastqc`		(to be able to run fastqc from anywhere)

4. **Install Python 2** (for DNApi)
    - `apt-get install python2.7`
    - `sudo ln -s /usr/bin/python2.7 /usr/bin/python`
    - test: `python --version`

5. **Install DNApi** (tested with v1.1)
    - Download the DNApi-master.zip from the DNApi github
    - Unzip, move into the created directory
    - test: `python dnapi.py -h`

6. **Install CURL** (for PIP)
    - `apt-get install curl`

7. **Install PIP** (for Cutadapt)
    - `curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py`
    - `python get-pip.py`

8. **Install Cutadapt** (tested with v1.18)
    - `pip install cutadapt`
    - test: `cutadapt -h`

9. **Install Perl** (for miRDeep2)
    - `apt-get install perl`
    - test: `perl --version`

10. **Install miRDeep2** (tested with v0.1.0)
    - Download the zip file from github
    - Unpack and get in the folder
    - `perl install.pl` (follow instructions)
    - test: get into the tutorial_dir/ and: `bash run_tut.sh`

11. **Install swan**
    - Download the precompiled binary file from the EMBL-EBI Kraken directory (Reaper-Tally-Minion package), and save as "swan"
    - `sudo ln -s /path/to/swan /usr/local/bin/swan`  (to be able to be run from anywhere)

12. **Install SRAtoolkit** (for downloader script, tested with v2.9.6)
    - Download and unzip into the directory of choice

13. **Acquire a reference Genome**		(tested predominatly with human genome GRCh37 and GRCh38)
    - Download the reference genome to the directory of choice
    - Use the bowtie tool from miRDeep2/essentials to create an index for your genome
    - `/path/to/miRDeep2/essentials/bowtie-1.1.1/bowtie-build reference_genome outfile_base`	(for more info please visit the bowtie 1.x manual)

14. **Download the known miRNAs from miRBase**		(tested with miRBase v20,21,22)
    - visit the miRBase website and download the hairpin and mature sequences

15. **Adjust the path variables in the Configuration file**
    - Inside the "diana_map_config.r" file adjust the path variables with the proper locations on your machine for: FastQC, DNApi, Cutadapt, mirDeep2 (and its aligner executable), the known_adapter_library.fa file, your indexed reference genome, the miRBase hairpin and mature fasta files.
    - Make sure that the path to mirDeep2/essentials folder is in your path variable.
    - If needed for downloading samples from SRA, adjust the path to your sratoolkit folder inside the "diana_map_downloader.r" file.



## DIANA_mAP Downloader script:

A DIANA-mAP R-script to automatically download biological datasets from various online databases using their Accession Numbers (Tested using R version 3.6). Additionally, instead of standalone accession numbers in the terminal, a text file can be provided containing one accession number per line.
Current databases supported: NCBI GEO, NCBI SRA, ENCODE.

Note: Absolute path to sratoolkit needs to be set manually through line 2 in the script.
```
Usage: 
Rscript diana_map_downloader.r [-d <output_path>] [--geofilter <filter_regex>] [--encodefilter <filter_type>] [--srafilter <filter_type>] (-i <input>...)

Options:

	-d <output_path>, --directory <output_path>
		Output saving path [default: ./]
		
	--geofilter <filter_regex>
		Filter downloaded GEO files by regular expression ex. BED [default: NULL] (only works for R version > 3.5.0)
		
	--encodefilter <filter_type>
		Filter downloaded ENCODE files ex. fastq [default: all]
		
	--srafilter <filter_type>
		File type to be download ex.fastq or sam to call the corresponding sratoolkit dumper. Only fastq and sam files supported. [default: fastq]
	
	-h --help
		Show the help screen```

Arguments:
	-i <input>...
		One or more file accession numbers to download. Also, a file path may be provided including one accession number per line. Τhe -i flag is needed before each argument.
```
**Example:**
`Rscript diana_map_downloader.r -d /home/ -i SRR033731 -i GSM508011 -i ENCFF522BMH`



## DIANA_mAP script:

The main DIANA_mAP script used for the analysis of a single sample. The file input can be either in .fq/.fastq format or in compressed .gz format.

The diana_map_config.r file provides the parameters for the analysis including the required software paths. However the user can provide a configuration file of their own that matches the structure of the diana_map_config.r file using the appropriate option (-c). Please use one of the two ways to provide the appropriate paths required.

Finally almost all of the parameters, except the software paths, can be provided using the option flags when calling the script as explained bellow.

```
Usage: 
diana_map.r [options] input_file output_directory


Options:
	-a ADAPTER, --adapter=ADAPTER
		The adapter sequence used for the dataset preparation. 
		Leave empty if the sequence is unknown. 
		Change to 'Raw_Input' if there is no adapter sequence present 
		in the dataset. (default = '')

	-e EXP_NAME, --exp_name=EXP_NAME
		A name for the experiment of the specific sample. 
		Text. (default = '')

	-p BOOLEAN, --preprocess=BOOLEAN
		Apply the pre-process analysis step. 
		Boolean. (default = 'TRUE')

	-y ADAPTER_TYPE, --adapter_type=ADAPTER_TYPE
		Adapter orientation information. 
		'a' for 3' adapters or 'g' for 5' adapters (default = 'a')

	-q NUMBER, --trim_quality=NUMBER
		The minimum base quality allowed. 
		Integer values from 0. (default = 10)

	-r NUMBER, --trim_max_err_rate=NUMBER
		The percentage of the error rate defining the allowed mismatches 
		between the adapter sequence and the read during the trimming process. 
		Float values from 0 to 1. (default = 0.1)

	-v NUMBER, --trim_min_adapter_overlap=NUMBER
		The minimum overlap with adapter sequence required (in bases) 
		to trim a sequence. 
		Integer values from 0. (default = 3)

	-l NUMBER, --trim_min_len=NUMBER
		The minimum read length accepted after the trimming process. 
		Integer values from 0. (default = 18)

	-b PATH, --adapter_library=PATH
		The position of the adapter library fasta file. 
		File path. (default = '')

	-t NUMBER, --adapter_identity_threshold=NUMBER
		The minimum percentage of identity required for a match 
		of the inferred adapter to an adapter sequence 
		in the adapter library. 
		Integer values from 0 to 100. (default = 90)

	-k NUMBER, --adapter_kmer_size=NUMBER
		The length of the k-mers exracted from the adapter sequence 
		for the k-mer loop pre-processing stage. 
		Integer values from 0. (default = 10)

	-m NUMBER, --max_multimaps=NUMBER
		The maximum number of multimaps allowed for each read 
		during the mapping process, reads with more maps 
		will be discarded. 
		Integer values from 0. (default = 5)

	-j NUMBER, --threads=NUMBER
		The number of threads to be used. (default = 1)

	-g PATH, --path_ref_genome_prefix=PATH
		The indexed reference genome directory plus the prefix. 
		Text. (default = '')

	-i BOOLEAN, --large_index=BOOLEAN
		Set TRUE if the genome index ends in '.ebwtl',
		set to FALSE if it ends in '.ebwt'
		Boolean. (default = TRUE)

	-s SPECIES, --species=SPECIES
		The three-letter species notation for the 
		species used in the analysis. 
		Text. (default = 'hsa')

	-n PATH, --path_to_hairpin=PATH
		The miRBase hairpin miRNAs file location. 
		Text. (default = '')

	-u PATH, --path_to_mature=PATH
		The miRBase mature miRNAs file location. 
		Text. (default = '')

	-c PATH, --user_config=PATH
		A configuration file provided by the user. 
		Path. (default = '')

	-h, --help
		Show this help message

The variable loading order is as follows: Command_Line_Arguments > User_Configuration_File_Arguments > Default_Configuration_File_Arguments
```
**Example:**
`Rscript diana_map.r -e 'experiment_1' -j 2 -l 16 /home/data/SRR033731.fq /home/results/`


**Warning:** If the pre-process step is omitted using the -p flag, the input has to be pre-processed (i.e. without contaminants) and cannot be in compressed .gz format, it has to be provided in .fq/.fastq format.




## DIANA_mAP Multiple script:

The DIANA_mAP multi script is used to run multiple pipelines in parallel using forking and possibly apply Differential Expression analysis to their results if the -d flag is given along with the conditions table file. The script invoces the diana_map.r script multiple times.

The "number_of_cores" specifies how many parallel pipelines will be executed simultaneously at each time.

The "pipeline_arguments" will be passed to every pipeline run with this script and can be in the form of: diana_map_multi.r "input_directory" "output_directory" "number_of_cores" "-c 'value' -m 12 -d 'value'" etc. If no parameters are to be passed to the pipeline, please use empty quotes to nominate that (eg. ""). For more information on the parameters available please use the help flag on the diana_map.r script.

**Warning:** If the "thread" parameter is provided for each pipeline and its value is above 1, the total number or threads that will be required for this script to run is: (number_of_cores * threads)! Please make sure your system has enough available resources. Requesting for more resources than available can lead to multiple fatal errors and is not recommended.
```
Usage: 
diana_map_multi.r [options] input_directory output_directory number_of_cores pipeline_arguments

Options:
	-d PATH, --difexpression=PATH
		Add this option to apply Differential Expression 
		analysis to the results of the samples run.
		Provide the path to the condition table for the 
		differential expression analysis.
		It should contain two columns with their headers,
		the first being the sample column containing all 
		the sample names and the second should be the 
		condition column, containing a condition for each sample. 
		NOTICE: The different conditions must be exactly two in number.

	-h, --help
		Show this help message and exit
```
**Example:**
`Rscript diana_map_multi.r -d /home/DIANA_mAP/example_DE_conditions_table.csv /home/test_data/ /home/test_results/ 4 "-l 16 -e 'experiment_1'"`

## Results
The results are organized into folders for each sample using the input file name, any value provided in the '-e' parameter (DIANA-mAP script) and the timestamp. Upon the successful completion of the analysis each sample folder contains:
  - A **Counts.txt result file** containing per known miRNA the raw quantified reads and the normalized quantified reads (RPM)
  - A **Final_Report.txt file** containing a summary of the analysis statistics
  - A **Pie chart** showing the distribution of reads after the pre-processing step of the analysis
  - Two **Sequence Length Distribution graphs**, one before and one after the pre-processing step
  - A **QC folder** containing the FastQC output files depicting the sample quality before and after pre-processing
  - An **Intermediate Results folder** produced during the various analysis steps

## KNOWN ISSUES

- `dnapi.py: error: [Errno 28] No space left on device`
>When analyzing multiple samples in parallel, sometimes python gets overloaded and crashes providing the error above. Please set the parameter 'number_of_cores' of diana_map_multi.r to 5 or less and run the analysis on failed samples again.
This issue has been occasionally observed when multiple large samples (5 million reads or more) are analyzed simultaneously.
