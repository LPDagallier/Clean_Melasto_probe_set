#!/bin/bash

############      SLURM CONFIGURATION      ###################
#SBATCH --job-name=prepare_1KP
#SBATCH --account=soltis
#SBATCH --qos=soltis
#SBATCH --cpus-per-task=8 
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16gb 
#SBATCH --time=24:00:00
#SBATCH --mail-user=<insert email>
#SBATCH --mail-type=ALL
############################################################

echo "JOB CONFIGURATION"
echo "Job ID: " $SLURM_JOB_ID
echo "Name of the job: " $SLURM_JOB_NAME
echo "Node(s) allocated to the job: " $SLURM_JOB_NODELIST
echo "Number of nodes allocated: " $SLURM_JOB_NUM_NODES
echo "Memory allocated per node: " $SBATCH_MEM_PER_NODE
echo "Number of CPU tasks: " $SLURM_NTASKS
echo "Amount of CPU per tasks: " $SLURM_CPUS_PER_TASK
echo "Amount of RAM per tasks: " $SBATCH_MEM_PER_CPU
echo "Directory from which sbatch was invoked: " $SLURM_SUBMIT_DIR
echo "Temporary folder in which the job runs: " $SLURM_TMPDIR

###################################################
#### 0 Preparation of files and transfer to cluster
###################################################

#### Set up PATHS and VARIABLES
# change depending on what data is being analysed

path_to_dir_in="/blue/soltis/dagallierl/DATASETS/PHYLOGENOMICS/onekp/v2";

# OUTPUT FOLDER
# where the ouptputs will be stored
path_to_dir_out="/blue/soltis/dagallierl/DATA_ANALYSES/PROBE_SET_CLEANING/prepare1KP_$SLURM_JOB_ID/";
# make output folder in home directory
mkdir $path_to_dir_out

#temporary folder (intermediate files)
#make temporary directory to store files/run analyses in
# this will be automatically deleted at the end
path_to_tmp=$SLURM_TMPDIR

cd $path_to_tmp

#module load parallel

# COPY FILES
echo "COPYING FILES";

scp $path_to_dir_in/SWGX-translated-reference-names.tsv $path_to_tmp
scp $path_to_dir_in/WWQZ-translated-reference-names.tsv $path_to_tmp
scp $path_to_dir_in/SWGX-translated-nucleotides.fa $path_to_tmp
scp $path_to_dir_in/SWGX-translated-protein.fa $path_to_tmp
scp $path_to_dir_in/WWQZ-translated-nucleotides.fa $path_to_tmp
scp $path_to_dir_in/WWQZ-translated-protein.fa $path_to_tmp

echo "DONE COPYING ALL FILES";

# LOAD MODULES
module load seqkit/2.0.0
module load ncbi_blast/2.10.1
module load csvtk/0.23.0

###################################################
#### 1 Group the sequences according to their ortholog group
###################################################
cd $path_to_tmp

echo "Start grouping the sequences according to their ortholog group...";
mkdir orthogroups
scp $path_to_dir_in/orthogroups/* ./orthogroups
echo "Done grouping the sequences according to their ortholog group...";

###################################################
#### 1bis Filter out sequences less than 30 AA (90 nt)
###################################################

cd orthogroups
seqkit seq -w0 -m 90 SWGX-WWQZ-translated-nucleotides-grouped.FNA > SWGX-WWQZ-translated-nucleotides-grouped.FNA.temp
seqkit seq -w0 -m 30 SWGX-WWQZ-translated-prot-grouped.FAA > SWGX-WWQZ-translated-prot-grouped.FAA.temp
rm -r SWGX-WWQZ-translated-prot-grouped.FAA SWGX-WWQZ-translated-nucleotides-grouped.FNA
rename -v '.FAA.temp' '.FAA' ./*
rename -v '.FNA.temp' '.FNA' ./*
cd $path_to_tmp

###################################################
#### 2 Build BLAST databases
###################################################

echo "Starting building BLAST databases...";
mkdir BLAST_custom_DB
makeblastdb -in orthogroups/SWGX-WWQZ-translated-nucleotides-grouped.FNA -out BLAST_custom_DB/SWGX-WWQZ_nucleotides_grouped_BLAST_DB -dbtype nucl -parse_seqids
makeblastdb -in orthogroups/SWGX-WWQZ-translated-prot-grouped.FAA -out BLAST_custom_DB/SWGX-WWQZ_protein_grouped_BLAST_DB -dbtype prot -parse_seqids
echo "Done building BLAST databases.";

####################################################
##### 3 Self-Blast1
####################################################

echo "Self-Blasting1..."
blastp -task blastp-fast -db BLAST_custom_DB/SWGX-WWQZ_protein_grouped_BLAST_DB -query orthogroups/SWGX-WWQZ-translated-prot-grouped.FAA -evalue 1e-06 -outfmt "10 std qframe qcovs qseq" -num_threads $SLURM_CPUS_PER_TASK > selfblast_SWGX-WWQZ_protein_grouped_blastp_out_SB1.csv
cat selfblast_SWGX-WWQZ_protein_grouped_blastp_out_SB1.csv | awk -F ',' 'OFS = "," {
temp1 = $1
temp2 =$2
sub(".+-","", temp1) 
sub(".+-","", temp2)
if (temp1 != temp2)
print $0
}' > different_query_and_subject_loci_SB1.csv
# cat different_query_and_subject_loci_SB1.csv | csvtk filter -H -f'3>60' | csvtk filter -H -f'4>50' | csvtk filter -H -f'14>50' | csvtk filter2 -H -f '$12 > 100' > filtered_different_query_and_subject_loci_SB1.csv
cat different_query_and_subject_loci_SB1.csv | csvtk filter2 -H -f '$3 > 60 && $4 > 50 && $14 > 50 && $12 > 100 || $3 > 80 && $4 > 50 && $12 > 50' > filtered_different_query_and_subject_loci_SB1.csv
echo "Done Self-Blast1."


####################################################
##### 4 Associate the loci according to the self blast
####################################################

echo "Starting associating the loci according to the self blast..."
echo "Starting associating loci names..."

cat filtered_different_query_and_subject_loci_SB1.csv | awk -F ',' '{
query = $1
subject = $2
sub(".+-","", query) 
sub(".+-","", subject)
print query"," subject
}' > pairs_SB1.csv

query_list=$(cat pairs_SB1.csv | awk -F ',' '{print $1}' | sort | uniq)
:> selfblast_table.txt
#echo "locusCustomID, correspondingLoci" > selfblast_table.txt
#query='LOC00745'
#query='LOC01143'
#query='LOC01796'
for query in $query_list
do
  echo $query
  subject=$(grep -e $query, pairs_SB1.csv | awk -F ',' 'ORS = "\n" {print $2}' | sort -n | uniq | tr '\n' ',' | sed 's/,$/\n/')
  echo $subject
  to_search=$(echo $query,$subject | sed 's/,/|/g')
  #egrep -qs -e $to_search selfblast_table.txt
  if ! egrep -qs -e $to_search selfblast_table.txt; then
    echo Query or Subject not in the table
    echo $query,$subject >> selfblast_table.txt
  else
    echo Query or Subject YET in the table
    replacement=$(echo $(egrep -e $to_search selfblast_table.txt | sed -z 's/\n/,/g' | sed 's/.$/\n/'),$query,$subject | sed -e $'s/,/\\\n/g' | sort -n | uniq | tr '\n' ',' | sed 's/.$/\n/')
    line_numbers=$(egrep -n -e $to_search selfblast_table.txt | sed 's/:.*//g')
    #replacement=$(echo $query','$subject | sed -e $'s/,/\\\n/g' | sort -n | uniq | tr '\n' ',' | sed 's/.$/\n/')
    for N in $line_numbers;
    do
      sed -i $N's/.*/'$replacement'/g' selfblast_table.txt
    done
    cat selfblast_table.txt | sort | uniq > selfblast_table.txt.temp
    rm selfblast_table.txt 
    mv selfblast_table.txt.temp selfblast_table.txt 
  fi
done
printf '\a'
cat selfblast_table.txt > selfblast_table_fromSB1.txt #intermediate backup
echo "Done associating loci names..."

#check for no duplicates:
echo Total number of IDs: $(cat selfblast_table.txt | sed 's/,/\n/g' | wc -l)
echo "Number of IDs (dups removed):" $(cat selfblast_table.txt | sed 's/,/\n/g' | sort | uniq | wc -l)

echo "Starting associating loci in .FNA and .FAA files..."
:> selfblast_table_paired.txt
for line in $(cat selfblast_table.txt)
do
#  echo $line
  locus=$(echo $line | cut -d ',' -f 1)
  echo $locus
  others=$(echo $line | sed 's/'$locus',//g' | sed -e $'s/,/\\\n/g')
  for other in $others
  do
    echo $locus,$other >> selfblast_table_paired.txt
  done
done
cat selfblast_table_paired.txt > selfblast_table_paired_fromSB1.txt #intermediate backup
cat selfblast_table_paired.txt | awk -F ',' '{print $2"\t"$1}' > replacement.txt
seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key orthogroups/SWGX-WWQZ-translated-prot-grouped.FAA > orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA
seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key orthogroups/SWGX-WWQZ-translated-nucleotides-grouped.FNA > orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA
#printf '\a'

# sort by locusCustomID
sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA.temp
sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA.temp
rm orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA
rename -v '.FAA.temp' '.FAA' ./orthogroups/*
rename -v '.FNA.temp' '.FNA' ./orthogroups/*
echo "Done associating loci in .FNA and .FAA files..."

echo "Done associating the loci according to the self blast."

###################################################
#### 5 Build BLAST databases over previous databases
###################################################

echo "Starting building SB1 BLAST databases...";
makeblastdb -in orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA -out BLAST_custom_DB/SWGX-WWQZ_nucleotides_grouped_SB1_BLAST_DB -dbtype nucl -parse_seqids
makeblastdb -in orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA -out BLAST_custom_DB/SWGX-WWQZ_protein_grouped_SB1_BLAST_DB -dbtype prot -parse_seqids
echo "Done building SB1 BLAST databases.";

######################################################################################################
#### 5bis intermediate transfer
# copy everything to output folder
scp -rp $path_to_tmp/* $path_to_dir_out/
######################################################################################################

# 
# ####################################################
# ##### 6 Self-Blast2
# ####################################################
# 
# echo "Self-Blasting2..."
# blastp -task blastp-fast -db BLAST_custom_DB/SWGX-WWQZ_protein_grouped_SB1_BLAST_DB -query orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA -evalue 1e-06 -outfmt "10 std qframe qcovs qseq" -num_threads $SLURM_CPUS_PER_TASK > selfblast_SWGX-WWQZ_protein_grouped_blastp_out_SB2.csv
# 
# cat selfblast_SWGX-WWQZ_protein_grouped_blastp_out_SB2.csv | awk -F ',' 'OFS = "," {
# temp1 = $1
# temp2 =$2
# sub(".+-","", temp1) 
# sub(".+-","", temp2)
# if (temp1 != temp2)
# print $0
# }' > different_query_and_subject_loci_SB2.csv
# cat different_query_and_subject_loci_SB2.csv | csvtk filter -H -f'3>60' | csvtk filter -H -f'4>50' | csvtk filter -H -f'14>50' | csvtk filter2 -H -f '$12 > 100' > filtered_different_query_and_subject_loci_SB2.csv
# echo "Done Self-Blast2."
# 
# ####################################################
# ##### 7 Associate the loci according to the self blast
# ####################################################
# 
# echo "Starting associating the loci according to the self blast..."
# echo "Starting associating loci names..."
# 
# cat filtered_different_query_and_subject_loci_SB2.csv | awk -F ',' '{
# query = $1
# subject = $2
# sub(".+-","", query) 
# sub(".+-","", subject)
# print query"," subject
# }' > pairs_SB2.csv
# 
# query_list=$(cat pairs_SB2.csv | awk -F ',' '{print $1}' | sort | uniq)
# #query="LOC00013"
# for query in $query_list
# do
#   echo $query
#   subject=$(grep -e $query, pairs_SB2.csv | awk -F ',' 'ORS = "\n" {print $2}' | sort -n | uniq | tr '\n' ',' | sed 's/,$/\n/')
#   echo $subject
#   to_search=$(echo $query,$subject | sed 's/,/|/g')
#   #egrep -qs -e $to_search selfblast_table.txt
#   if ! egrep -qs -e $to_search selfblast_table.txt; then
#     echo Query or Subject not in the table
#     echo $query,$subject >> selfblast_table.txt
#   else
#     echo Query or Subject YET in the table
#     replacement=$(echo $(egrep -e $to_search selfblast_table.txt | sed -z 's/\n/,/g' | sed 's/.$/\n/'),$query,$subject | sed -e $'s/,/\\\n/g' | sort -n | uniq | tr '\n' ',' | sed 's/.$/\n/')
#     line_numbers=$(egrep -n -e $to_search selfblast_table.txt | sed 's/:.*//g')
#     for N in $line_numbers;
#     do
#       sed -i $N's/.*/'$replacement'/g' selfblast_table.txt
#     done
#     cat selfblast_table.txt | sort | uniq > selfblast_table.txt.temp
#     rm selfblast_table.txt 
#     mv selfblast_table.txt.temp selfblast_table.txt 
#   fi
# done
# printf '\a'
# cat selfblast_table.txt > selfblast_table_fromSB2.txt #intermediate backup
# 
# echo "Done associating loci names..."
# 
# #check for no duplicates:
# echo Total number of IDs: $(cat selfblast_table.txt | sed 's/,/\n/g' | wc -l)
# echo "Number of IDs (dups removed):" $(cat selfblast_table.txt | sed 's/,/\n/g' | sort | uniq | wc -l)
# 
# echo "Starting associating loci in .FNA and .FAA files..."
# for line in $(cat selfblast_table.txt)
# do
# #  echo $line
#   locus=$(echo $line | cut -d ',' -f 1)
#   echo $locus
#   others=$(echo $line | sed 's/'$locus',//g' | sed -e $'s/,/\\\n/g')
#   for other in $others
#   do
#     echo $locus,$other >> selfblast_table_paired.txt
#   done
# done
# cat selfblast_table_paired.txt > selfblast_table_paired_fromSB2.txt #intermediate backup
# 
# cat selfblast_table_paired.txt | awk -F ',' '{print $2"\t"$1}' > replacement.txt
# seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA > orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB2.FNA
# seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA > orthogroups/SWGX-WWQZ-translated-prot-grouped_SB2.FAA
# #printf '\a'
# 
# # sort by locusCustomID
# sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' orthogroups/SWGX-WWQZ-translated-prot-grouped_SB2.FAA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > orthogroups/SWGX-WWQZ-translated-prot-grouped_SB2.FAA.temp
# sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB2.FNA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB2.FNA.temp
# rm orthogroups/SWGX-WWQZ-translated-prot-grouped_SB2.FAA orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB2.FNA
# # rm orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB2.FNA
# rename -v '.FAA.temp' '.FAA' ./orthogroups/*
# rename -v '.FNA.temp' '.FNA' ./orthogroups/*
# echo "Done associating loci in .FNA and .FAA files..."
# 
# echo "Done associating the loci according to the self blast."
# 
# ###################################################
# #### 8 Build BLAST databases over previous databases
# ###################################################
# 
# echo "Starting building SB2 BLAST databases...";
# # seqkit translate -w0 orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB2.FNA > orthogroups/SWGX-WWQZ-translated-prot-grouped_SB2.FAA
# makeblastdb -in orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB2.FNA -out BLAST_custom_DB/SWGX-WWQZ_nucleotides_grouped_SB2_FINAL_BLAST_DB -dbtype nucl -parse_seqids
# makeblastdb -in orthogroups/SWGX-WWQZ-translated-prot-grouped_SB2.FAA -out BLAST_custom_DB/SWGX-WWQZ_protein_grouped_SB2_FINAL_BLAST_DB -dbtype prot -parse_seqids
# echo "Done building SB2 BLAST databases.";

####################################################
##### 9 Check Self-Blast (SB2)
####################################################

echo "Starting check self-blast (SB2)..."
blastp -task blastp-fast -db BLAST_custom_DB/SWGX-WWQZ_protein_grouped_SB1_BLAST_DB -query orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA -evalue 1e-06 -outfmt "10 std qframe qcovs qseq" -num_threads $SLURM_CPUS_PER_TASK > check_selfblast_SWGX-WWQZ_protein_grouped_blastp_out.csv
cat check_selfblast_SWGX-WWQZ_protein_grouped_blastp_out.csv | awk -F ',' 'OFS = "," {
temp1 = $1
temp2 =$2
sub(".+-","", temp1)
sub(".+-","", temp2)
if (temp1 != temp2)
print $0
}' > check_different_query_and_subject_loci.csv
cat check_different_query_and_subject_loci.csv | csvtk filter -H -f'3>60' | csvtk filter -H -f'4>50' | csvtk filter -H -f'14>50' | csvtk filter2 -H -f '$12 > 100' > check_filtered_different_query_and_subject_loci.csv
echo "Done check self-blast."

##################################################
#### 10 Transfer & clean
##################################################

echo "Transfert data node -> blue storage";

# make output folder in home directory
#mkdir $path_to_dir_out

# copy everything to output folder
scp -rp $path_to_tmp/* $path_to_dir_out/

echo "done moving"
echo "FINISHED"

# All the data are automatically deleted on the node by SLURM