#!/bin/bash

############      SLURM CONFIGURATION      ###################
#SBATCH --job-name=prepare_full_reference_DB
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
# SLURM_CPUS_PER_TASK=8
path_to_onekp="/blue/soltis/dagallierl/DATA_ANALYSES/PROBE_SET_CLEANING/prepare1KP_51843049";
#Angio353_FNA="/blue/soltis/dagallierl/DATASETS/PHYLOGENOMICS/target_references/Angio353_oneline.fa";
Angio353_FNA="/blue/soltis/dagallierl/DATASETS/PHYLOGENOMICS/target_references/mega353_Myrtales.fa";

prot_DB=$path_to_onekp/BLAST_custom_DB/SWGX-WWQZ_protein_grouped_SB1_BLAST_DB

# OUTPUT FOLDER
# where the ouptputs will be stored
path_to_dir_out="/blue/soltis/dagallierl/DATA_ANALYSES/PROBE_SET_CLEANING/full_reference_DB_$SLURM_JOB_ID/";
#path_to_dir_out="/blue/soltis/dagallierl/DATA_ANALYSES/PROBE_SET_CLEANING/full_reference_DB_test2";

#temporary folder (intermediate files)
#make temporary directory to store files/run analyses in
# this will be automatically deleted at the end
path_to_tmp=$SLURM_TMPDIR
# mkdir $path_to_dir_out
# path_to_tmp=$path_to_dir_out

cd $path_to_tmp

#module load parallel

# COPY FILES
echo "COPYING FILES";

cat $path_to_onekp/orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA > FULL_REFERENCE_nucl.FNA
cat $path_to_onekp/orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA > FULL_REFERENCE_prot.FAA
nuc_FNA="FULL_REFERENCE_nucl.FNA"
prot_FAA="FULL_REFERENCE_prot.FAA"
dos2unix $Angio353_FNA
# cat $Angio353_FNA | sed 's/-v/_v/g' | sed 's/\./_/g' | sed 's/-ASGPB/_ASGPB/g' | sed 's/-iTAG/_iTAG/g' | sed 's/-Genoscope/_Genoscope/g' > Angio353.FNA
cat $Angio353_FNA | sed 's/ .*//g' | sed 's/-[[:digit:]]\{4\}_exonerate_grafted/_EG/g' | sed 's/grafted_with/GW/g' | sed 's/-v/_v/g' | sed 's/\./_/g' | sed 's/-ASGPB/_ASGPB/g' | sed 's/-iTAG/_iTAG/g' | sed 's/-Genoscope/_Genoscope/g' > Angio353.FNA
Angio353_FNA=$path_to_tmp/Angio353.FNA
echo "DONE COPYING ALL FILES";

# LOAD MODULES
module load seqkit/2.0.0
module load csvtk/0.23.0
module load ncbi_blast/2.10.1

###################################################
#### 1 Associate Angio353 loci with SWGX and WWQZ transcriptomes
###################################################
cd $path_to_tmp

echo "Start associating Angio353 loci with SWGX and WWQZ transcriptomes...";

seqkit seq -n $Angio353_FNA | sed 's/.*-//g' | sort | uniq > Angio353_list.txt
:> Angio353_to_custom_loci_table.txt
# locus="5596"
# locus="4932"
for locus in $(cat Angio353_list.txt)
do
  echo STARTING FOR $locus
  mkdir $locus
  cd $locus
  seqkit grep -w0 -nrp $locus $Angio353_FNA > $locus.FNA
  query=$locus".FNA"
  cat $query | seqkit translate -w0 > $locus".FAA"
  blastp -task blastp-fast -matrix BLOSUM62 -db $prot_DB -query $locus".FAA" -evalue 1e-06 -outfmt "10 std qframe qcovs qseq sframe sseq" -num_threads $SLURM_CPUS_PER_TASK > $locus"_blastp_out_unfiltered.csv"
  #blastx -task blastx-fast -matrix BLOSUM62 -db $prot_DB -query $query -culling_limit 1 -evalue 1e-06 -outfmt "10 std qframe qcovs qseq sframe sseq" -num_threads $SLURM_CPUS_PER_TASK > $locus"_blastx_out_unfiltered.csv"
  if [ -s $locus"_blastp_out_unfiltered.csv" ]; then
    #cat $locus"_blastx_out_unfiltered.csv" | csvtk filter -H -f"3>60" > $locus"_blastx_out.csv"
    cat $locus"_blastp_out_unfiltered.csv" | csvtk filter2 -H -f '$3 > 60 && $4 > 50 && $14 > 50 && $12 > 100 && $13 > 0' > $locus"_blastp_out.csv"
  fi
  if [ -s $locus"_blastp_out.csv" ]; then
    cat $locus"_blastp_out.csv" | awk -F ',' '{print $2}' | sed 's/.*-LOC/LOC/g' | sort | uniq > DB_locus
    nb_DB_loci=$(cat DB_locus | wc -l)
    if [ $nb_DB_loci -gt 1 ]; then
      echo WARNING: several DB loci found for locus $locus. No transcriptome appended. You might want to check manually for $locus.
      cat $locus".FNA" >> ../$nuc_FNA
      cat $locus".FAA" >> ../$prot_FAA
      echo $locus multiplematches >> ../Angio353_to_custom_loci_table.txt
    else
      DB_locus=$(cat DB_locus)
      cat $locus".FNA" | sed '/^>/s/-/__/g' | sed '/^>/s/$/__-'$DB_locus'/g' >> ../$nuc_FNA
      cat $locus".FAA" | sed '/^>/s/-/__/g' | sed '/^>/s/$/__-'$DB_locus'/g' >> ../$prot_FAA
      echo $locus $(cat DB_locus) >> ../Angio353_to_custom_loci_table.txt
    fi
  else
    echo No match found for $locus.
    cat $locus".FNA" >> ../$nuc_FNA
    cat $locus".FAA" >> ../$prot_FAA
    echo $locus nomacth >> ../Angio353_to_custom_loci_table.txt
  fi
cd ..
done
# printf '\a'

echo "Done associating Angio353 loci with SWGX and WWQZ transcriptomes...";

###################################################
#### 2 Build BLAST databases
###################################################
echo "Starting building BLAST databases...";
mkdir FULL_REFERENCE_BLAST_DB
makeblastdb -in FULL_REFERENCE_nucl.FNA -out FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_nucl_BLAST_DB -dbtype nucl -parse_seqids
makeblastdb -in FULL_REFERENCE_prot.FAA -out FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_prot_BLAST_DB -dbtype prot -parse_seqids
echo "Done building BLAST databases.";

####################################################
##### 3 Self-Blast
####################################################

echo "Self-Blasting..."
blastp -task blastp-fast -db FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_prot_BLAST_DB -query FULL_REFERENCE_prot.FAA -evalue 1e-06 -outfmt "10 std qframe qcovs qseq" -num_threads $SLURM_CPUS_PER_TASK > full_ref_blastp_out.csv
cat full_ref_blastp_out.csv | awk -F ',' 'OFS = "," {
temp1 = $1
temp2 =$2
sub(".+-","", temp1) 
sub(".+-","", temp2)
if (temp1 != temp2)
print $0
}' > different_query_and_subject_loci.csv
cat different_query_and_subject_loci.csv | csvtk filter2 -H -f '$3 > 60 && $4 > 50 && $14 > 50 && $12 > 100 && $13 > 0' > filtered_different_query_and_subject_loci.csv
echo "Done Self-Blast."

####################################################
##### 4 Associate the loci according to the self blast
####################################################

echo "Starting associating the loci according to the self blast..."
echo "Starting associating loci names..."

cat filtered_different_query_and_subject_loci.csv | awk -F ',' '{
query = $1
subject = $2
sub(".+-","", query) 
sub(".+-","", subject)
print query"," subject
}' | sort | uniq > pairs.csv

query_list=$(cat pairs.csv | awk -F ',' '{print $1}' | sort | uniq)
:> selfblast_table.txt
#query="5899"
for query in $query_list
do
  echo $query
  subject=$(grep -e $query, pairs.csv | awk -F ',' 'ORS = "\n" {print $2}' | sort -n | uniq | tr '\n' ',' | sed 's/,$/\n/')
  echo $subject
  to_search=$(echo $query,$subject | sed 's/,/|/g')
  if ! egrep -qs -e $to_search selfblast_table.txt; then
    echo Query or Subject not in the table
    echo $query,$subject >> selfblast_table.txt
  else
    echo Query or Subject YET in the table
    replacement=$(echo $(egrep -e $to_search selfblast_table.txt | sed -z 's/\n/,/g' | sed 's/.$/\n/'),$query,$subject | sed -e $'s/,/\\\n/g' | sort -n | uniq | tr '\n' ',' | sed 's/.$/\n/')
    line_numbers=$(egrep -n -e $to_search selfblast_table.txt | sed 's/:.*//g')
    for N in $line_numbers;
    do
      sed -i $N's/.*/'$replacement'/g' selfblast_table.txt
    done
    cat selfblast_table.txt | sort | uniq > selfblast_table.txt.temp
    rm selfblast_table.txt 
    mv selfblast_table.txt.temp selfblast_table.txt 
  fi
done

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

:> replacement.txt
# cat selfblast_table_paired.txt | awk -F ',' '{print "-"$2"\t-"$1}'
cat selfblast_table_paired.txt | grep ',[[:digit:]]\{4\}$' | awk -F ',' '{print $2"\t__"$2"__-"$1}' >> replacement.txt
cat selfblast_table_paired.txt | grep -v ',[[:digit:]]\{4\}$' | awk -F ',' '{print $2"\t"$1}' >> replacement.txt
seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key FULL_REFERENCE_prot.FAA | seqkit replace -w0 -p '-_' -r '_' > FULL_REFERENCE_prot_SB1.FAA
seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key FULL_REFERENCE_nucl.FNA | seqkit replace -w0 -p '-_' -r '_' > FULL_REFERENCE_nucl_SB1.FNA

# sort by locusCustomID
sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' FULL_REFERENCE_prot_SB1.FAA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > FULL_REFERENCE_prot_SB1.FAA.temp
sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' FULL_REFERENCE_nucl_SB1.FNA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > FULL_REFERENCE_nucl_SB1.FNA.temp
rm FULL_REFERENCE_prot_SB1.FAA FULL_REFERENCE_nucl_SB1.FNA
rename -v '.FAA.temp' '.FAA' ./*
rename -v '.FNA.temp' '.FNA' ./*

echo "Done associating loci in .FNA and .FAA files..."
echo "Done associating the loci according to the self blast."

###################################################
#### 5 Build BLAST databases over previous databases
###################################################

echo "Starting building SB1 BLAST databases...";
makeblastdb -in FULL_REFERENCE_nucl_SB1.FNA -out FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_nucl_BLAST_DB_FINAL -dbtype nucl -parse_seqids
makeblastdb -in FULL_REFERENCE_prot_SB1.FAA -out FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_prot_BLAST_DB_FINAL -dbtype prot -parse_seqids
echo "Done building SB1 BLAST databases.";

####################################################
##### 6 Check Self-Blast (SB2)
####################################################

echo "Starting check self-blast (SB2)..."
blastp -task blastp-fast -db FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_prot_BLAST_DB_FINAL -query FULL_REFERENCE_prot_SB1.FAA -evalue 1e-06 -outfmt "10 std qframe qcovs qseq" -num_threads $SLURM_CPUS_PER_TASK > check_selfblast_full_ref_blastp_out.csv

cat check_selfblast_full_ref_blastp_out.csv | awk -F ',' 'OFS = "," {
temp1 = $1
temp2 =$2
sub(".+-","", temp1)
sub(".+-","", temp2)
if (temp1 != temp2)
print $0
}' > check_different_query_and_subject_loci.csv
cat check_different_query_and_subject_loci.csv | csvtk filter2 -H -f '$3 > 60 && $4 > 50 && $14 > 50 && $12 > 100 && $13 > 0' > check_filtered_different_query_and_subject_loci.csv
echo "Done check self-blast."
#printf '\a'

##################################################
#### 4 Transfer & clean
##################################################

echo "Transfert data node -> blue storage";

# make output folder in home directory
mkdir $path_to_dir_out

# copy everything to output folder
scp -rp $path_to_tmp/* $path_to_dir_out/

echo "done moving"
echo "FINISHED"

# All the data are automatically deleted on the node by SLURM