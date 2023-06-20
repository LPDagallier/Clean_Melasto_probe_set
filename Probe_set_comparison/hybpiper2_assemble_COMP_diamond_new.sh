#!/bin/bash

############      SLURM CONFIGURATION      ###################
#SBATCH --job-name=hybpiper2_assemble_Probe_set_comparison_diamond_new
#SBATCH --account=<insert_your_account_ID_here>
#SBATCH --qos=<insert_your_QOS_here>
#SBATCH --cpus-per-task=8 
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64gb 
#SBATCH --time=4-00:00:00
#SBATCH --mail-user=<insert_your_email_here>
#SBATCH --mail-type=ALL
#SBATCH --output=slurm-%x-%j.out
############################################################

echo "JOB CONFIGURATION"
echo "Job ID: " $SLURM_JOB_ID
echo "Name of the job: " $SLURM_JOB_NAME
echo "Node allocated to the job: " $SLURM_JOB_NODELIST
echo "Number of nodes allocated to the job: " $SLURM_JOB_NUM_NODES
echo "Number of CPU tasks in this job: " $SLURM_NTASKS
echo "Directory from which sbatch was invoked: " $SLURM_SUBMIT_DIR
echo "Temporary folder in which the job runs: " $SLURM_TMPDIR

###################################################
#### 0 Preparation of files and transfer to cluster
###################################################

#### Set up PATHS and VARIABLES
# change depending on what data is being analysed
# The analyse identifier that will be used as a preffix or suffix in scripts, file names and/or folder names
# must be unique to THIS run of analysis (that is this combination of samples, reference file, HYbPiper parameters, etc.)
analysis_ID="Probe_set_comparison"
step_ID="hybpiper2_assemble_no_stiched_diamond_new"

# - Samples: Tibouchina and Memecylon samples from Jantzen et al 2020 (NCBI BioProject numbers PRJNA573947 and PRJNA576018)
# - Target file used: PROBE_SET_CLEAN_v4_prot.FAA
# - HybPiper version 2.1.1
# - BWA

# REFERENCE FILE:
path_to_ref=$path_to_wd/CLEAN_PROBE_SET/;
reference_fasta_file="PROBE_SET_CLEAN_v4_prot.FAA"

# DATA-RELATED FILES:
# Prepare (locally) the following files:
#namelist.txt: contains the sample names that will be analysed
#input_fastq: contains the copying comands for batch-copy the fastq files
#files_renaming.txt: contains the renaming comands for batch-rename the fastq files
#These files have to be in a folder that has the same name as analysis_ID situated in path_to_dir_in:
path_to_dir_in=$path_to_wd/Probe_set_comparison;
path_to_fastq=$path_to_wd/Probe_set_comparison/fastq;

# OUTPUT FOLDER
# where the ouptputs will be stored
path_to_dir_out=$path_to_wd/Probe_set_comparison/$analysis_ID"_"$step_ID"_"$SLURM_JOB_ID/;

#temporary folder (intermediate files)
#make temporary directory to store files/run analyses in
# this will be automatically deleted at the end
path_to_tmp=$SLURM_TMPDIR
# for checking at the begining:
# path_to_tmp=$path_to_dir_out
# mkdir $path_to_tmp

cd $path_to_tmp


# COPY FILES
echo "COPYING FILES";

echo "copying fasta reference"
scp $path_to_ref/$reference_fasta_file $path_to_tmp
echo "The reference file used is "$reference_fasta_file
echo "done copying fasta reference"

echo "copying data-related files"
scp $path_to_dir_in/$analysis_ID/namelist.txt $path_to_tmp

# format text files properly (sometimes windows does weird things on txt files, making them unix-incompatible)
dos2unix *.txt
echo "done copying data-related files"

#list the file in the log file to allow checks
echo "LIST OF THE FILES:"
ls

# load module HybPiper
module load hybpiper/2.1.1

###################################################
#### 1 HybPiper assemble 
###################################################
cd $path_to_tmp

echo "Starting HybPiper assemble";

while read name;
do hybpiper assemble -t_aa $reference_fasta_file -r $path_to_fastq/$name"_"*.fastq --prefix $name --cpu 8 --diamond --unpaired $path_to_fastq/$name.fastq --no_stitched_contig;
done < namelist.txt

echo "Done HybPiper assemble";

###################################################
#### 2 Summary statistics & visualizing results
###################################################

echo "Starting computing summary statistics for GENES";
hybpiper stats -t_dna $reference_fasta_file --seq_lengths_filename genes_sequences_lengths --stats_filename hybpiper_genes_statistics gene namelist.txt
echo "Done computing summary statistics";


echo "Starting visualizing results with HybPiper script - might lead to bad figure";
hybpiper recovery_heatmap --heatmap_dpi 300 --heatmap_filetype pdf --heatmap_filename recovery_heatmap_exons genes_sequences_lengths.tsv
echo "Done visualizing results";

##################################################
#### 5 Transfer
##################################################

echo "Transfert data node -> storage";

# make output folder in home directory
mkdir $path_to_dir_out

# copy everything, for testing
scp -rp $path_to_tmp/* $path_to_dir_out/

echo "done moving, FINISHED";

# All the data are automatically deleted on the node by SLURM