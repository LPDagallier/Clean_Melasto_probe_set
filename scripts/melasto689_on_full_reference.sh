#!/bin/bash

############      SLURM CONFIGURATION      ###################
#SBATCH --job-name=melasto689_on_full_ref
#SBATCH --account=soltis
#SBATCH --qos=soltis
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16gb 
#SBATCH --time=24:00:00
#SBATCH --mail-user=ldagallier@nybg.org
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
path_to_Melasto689="/blue/soltis/dagallierl/DATA_ANALYSES/PROBE_SET_CLEANING/Melastomataceae_689_original.FNA"
path_to_scripts="/blue/soltis/dagallierl/DATA_ANALYSES/PROBE_SET_CLEANING"

# to change:
path_to_full_ref="/blue/soltis/dagallierl/DATA_ANALYSES/PROBE_SET_CLEANING/full_reference_DB_51849446";
# sequences

# blast DB
prot_DB=$path_to_full_ref/FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_prot_BLAST_DB_FINAL
nucl_DB=$path_to_full_ref/FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_nucl_BLAST_DB_FINAL

# OUTPUT FOLDER
# where the ouptputs will be stored
path_to_dir_out="/blue/soltis/dagallierl/DATA_ANALYSES/PROBE_SET_CLEANING/melasto689_on_full_ref_$SLURM_JOB_ID/";
#path_to_dir_out="/blue/soltis/dagallierl/DATA_ANALYSES/PROBE_SET_CLEANING/melasto689_on_full_ref_51852846";

#temporary folder (intermediate files)
#make temporary directory to store files/run analyses in
# this will be automatically deleted at the end
path_to_tmp=$SLURM_TMPDIR
mkdir $path_to_dir_out
#path_to_tmp=$path_to_dir_out

cd $path_to_tmp
#module load parallel

# COPY FILES
echo "COPYING FILES";
scp -r $path_to_Melasto689 $path_to_tmp
scp -r $path_to_full_ref/FULL_REFERENCE_nucl.FNA $path_to_tmp
echo "DONE COPYING ALL FILES";

# LOAD MODULES
module load seqkit/2.0.0
module load csvtk/0.23.0
module load ncbi_blast/2.10.1
module load bedtools/2.30.0
module load R/4.2

###################################################
#### 1 Remove Angio353 template sequences from Melasto689
###################################################
cd $path_to_tmp
echo "Start removing Angio353 template sequences from Melasto689...";

seqkit grep -w 0 -n -r -p "[[:upper:]]{4}-" Melastomataceae_689_original.FNA > Angio353_from_Melasto689.FNA;
seqkit grep -w 0 -n -r -v -p "[[:upper:]]{4}-" Melastomataceae_689_original.FNA > Melasto689_trimed.FNA;
Melasto689="Melasto689_trimed.FNA"
seqkit seq -n Melasto689_trimed.FNA > template_sequences_list.txt
echo "Done removing Angio353 template sequences from Melasto689...";

###################################################
#### 2 Blastx each TS from Melasto689_trimmed and extract the hits
###################################################
echo "Start blastx each TS from Melasto689_trimmed and extract the hits..."
# extracted hits will be in extracted_hits
mkdir extracted_hits
mkdir nomatch_TS
mkdir NMR_hits
# TS="Brachyotum-5699"
# TS="Tibouchina1-5699"
#TS="Tibouchina2-5699"
#TS="Tibouchina-5554"
#TS="Tibouchina-NM_001198514.1"
#grep 'NM_124699.3' template_sequences_list.txt
#grep 'Tibouchina-5554' template_sequences_list.txt
for TS in $(cat template_sequences_list.txt);
do
  echo $TS
  mkdir $TS
  cd $TS
  seqkit grep -w0 -nrp $TS ../Melasto689_trimed.FNA > $TS.FNA
  blastx -task blastx-fast -matrix BLOSUM62 -db $prot_DB -query $TS.FNA -culling_limit 1 -evalue 1e-06 -ungapped -comp_based_stats F -outfmt "10 std qframe qcovs sframe qseq" > $TS"_blastx_out_unfiltered.csv"
  if [ -s $TS"_blastx_out_unfiltered.csv" ]; then
    cat $TS"_blastx_out_unfiltered.csv" | csvtk filter2 -H -f '$12 > 50 && $14 > 50 || $3 > 60 && $14 > 10 && $12 > 50' > $TS"_blastx_out.csv"
  fi
  if [ -s $TS"_blastx_out.csv" ]; then
    number_of_seqs_BLASTED=$(cat $TS"_blastx_out.csv" | wc -l)
  else
    number_of_seqs_BLASTED=(0)
  fi
  if [ $number_of_seqs_BLASTED -gt 0 ]; then
    for line in $(cat $TS"_blastx_out.csv");
    do
    # Retrieve the hit
      refLocus=$(echo $line | awk -F ',' '{print $2}' | sed 's/.*-//g')
      frame=$(echo $line | awk -F ',' '{print $13}')
      matchingSubject=$(echo $line | awk -F ',' '{print $2}')
      echo $line | awk -F ',' '{print $2}' >> matching_subject-$refLocus.txt
      #echo $line 
      #echo $refLocus
      #echo $frame
      if [ $frame -gt 0 ]; then # HIT IS REGULAR
        hit=$(echo $line | awk -F ',' '{print $1":"$7-1"-"$8}' | sed 's/-/__/g' | sed 's/$/-'$refLocus'/g')
        echo Hit is: $hit
        echo "...(regular hit)"
        echo $line | awk -F ',' '{print($1"\t"$7-1"\t"$8)}' | bedtools sort > $hit".bed"
        bedtools getfasta -bed $hit".bed" -fi $TS.FNA | sed 's/-/__/g' | sed '/^>/s/$/-'$refLocus'/g' > $hit"_hit.FNA"
      fi
      if [ $frame -lt 0 ]; then # HIT IS RC
        hit=$(echo $line | awk -F ',' '{print $1":"$8-1"-"$7}' | sed 's/-/__/g' | sed 's/$/-'$refLocus'/g')
        echo Hit is: $hit
        echo "... (reverse complemented hit)"
        echo $line | awk -F ',' '{print($1"\t"$8-1"\t"$7)}' | bedtools sort > $hit".bed"
        bedtools getfasta -bed $hit".bed" -fi $TS.FNA | seqkit seq -rp -t dna -v -w 0 | sed 's/-/__/g' | sed '/^>/s/$/-'$refLocus'/g' > $hit"_hit.FNA"
      fi
    done
    cp *_hit.FNA ../extracted_hits
    
    # Retrieve the NMR
    qcovs=$(cat $TS"_blastx_out.csv" | csvtk summary -H -n0 -f 14:sum)
    if [ $qcovs -lt 90 ]; then
      echo ...NMR found for $TS, extract NMR...
      cat *.bed | bedtools sort | bedtools merge > $TS"_merged_ALL.bed"
      seqkit fx2tab -ln  $TS".FNA" | awk '{print($1"\t"$2"\t"$2)}' | bedtools sort > $TS".bed"
      bedtools complement -L -i $TS"_merged_ALL.bed" -g $TS".bed" > $TS"_nonmatching_regions.bed"
      cat $TS"_nonmatching_regions.bed" >> ../TS_nonmatching_regions.bed
      bedtools getfasta -bed $TS"_nonmatching_regions.bed" -fi $TS".FNA" | seqkit seq -g -m 30 -w 0 | sed 's/-/_/g' | sed '/^>/s/$/-NMR/g' > $TS"NMR.FNA"
      if [ ! -s $TS"NMR.FNA" ]; then
        echo ... NMR less than 30 nucleotides, discard it...
        rm $TS"NMR.FNA"
      else
        cp *NMR.FNA ../NMR_hits
      fi
    fi
  else
    echo No matching hit for $TS
    cat $TS".FNA" > $TS"_nomatch.FNA"
    cp *_nomatch.FNA ../nomatch_TS
  fi
  cd ../
done
echo "Done blastx each TS from Melasto689_trimmed and extract the hits..."
printf '\a'

###################################################
#### 3 Merge the hits belogning to a same TS and a same locus
###################################################
echo "Start merging the hits belogning to a same TS and a same locus..."
# Cleaned hit will be stored in cleaned_hits
cd extracted_hits
mkdir cleaned_hits
echo "locusCustomID, correspondingLoci" > $path_to_tmp/loci_matching_table.csv

loci_list=$(ls -1 *.FNA | sed 's/.*-//g' | sed 's/_hit.FNA//g' | sort | uniq)

for locus in $loci_list
do
  # echo $locus
  # rm -r $locus
  mkdir $locus
  cp *-$locus"_hit.FNA" $locus
done

#locus="5699"
#locus="LOC08266"
for locus in $loci_list
do
  echo $locus
  cd $locus
  TS_list=$(ls -1 *hit.FNA | sed 's/:.*//g' | sort | uniq)
  TS_loci=$(ls -1 *hit.FNA | sed 's/:.*//g' | sed 's/.*__//g' | sort | uniq)
  echo $locus, $TS_loci >> $path_to_tmp/loci_matching_table.csv
  seqkit grep -w0 -nrp $locus $path_to_tmp/FULL_REFERENCE_nucl.FNA > $locus"_original.FNA"
  seqkit seq -n $locus"_original.FNA" > IDs_to_remove_after_alignment.txt
  
  #TS="Tetrazygia__NM_113708.2"
  for TS in $TS_list
  do
    echo $TS
    TS_locus=$(echo $TS | sed 's/.*__//g')
    nb_hits=$(ls -1 $TS*_hit.FNA | wc -l)
    if [ $nb_hits -gt 1 ]; then
      echo "More than 1 hit found for "$TS". Draw consensus sequence."
      cat $locus"_original.FNA" $TS*_hit.FNA > $TS"_hits_to_align.FNA"
    # Align the hits to the reference (as translated amino acids) and build consensus upon aligned hits
      Rscript $path_to_scripts/align_translated_to_ref_and_draw_consensus.R $TS"_hits_to_align.FNA" $TS"_hits_aligned_to_ref.FNA" IDs_to_remove_after_alignment.txt $TS"_hits_consensus" $TS"_hits_consensus.FNA"
    # Remove the gaps to clean the consensus sequence
      seqkit seq -g -w 0 $TS"_hits_consensus.FNA" | sed '/^>/s/_'$TS_locus'_/_'$TS_locus'__/g' | sed '/^>/s/$/'-$locus'/g' > $TS"_hits_consensus_clean-"$locus".FNA"
    elif [ $nb_hits -eq 1 ]; then
      echo "Only 1 hit found for "$TS"."
      seqkit replace -w 0 -p ':.*' -r '__hit' $TS*_hit.FNA | sed '/^>/s/$/'-$locus'/g' > $TS"_hit_clean-"$locus".FNA"
    fi
  done
  cp *clean-$locus.FNA ../cleaned_hits
  cd ..
done
printf '\a'
echo "Done merging the hits belogning to a same TS and a same locus..."

###################################################
#### 4 Create a sequence set with the loci to which the cleaned hits matched
###################################################
echo "Start creating a sequence set with the loci to which the cleaned hits matched..."

cd cleaned_hits
mkdir $path_to_tmp/cleaned_sequences_set
cat *.FNA > $path_to_tmp/cleaned_sequences_set/cleaned_hits_all.FNA

cd $path_to_tmp
mkdir cleaned_sequences_set/matching_subjects
find . -type f -name "matching_subject*.txt" -exec cp {} ./cleaned_sequences_set/matching_subjects \; # error messages are normal as the destination folder is embeded within the searching directory
cd cleaned_sequences_set/matching_subjects
for file in $(ls -1 *.txt)
do
  cat $file | sort | uniq > $file.temp
done
rm *.txt
rename -v '.txt.temp' '.txt' *

cd $path_to_tmp/cleaned_sequences_set
seqkit seq -n cleaned_hits_all.FNA | sed 's/.*-//g' | sort | uniq > loci_list.txt
cat cleaned_hits_all.FNA > final_sequences_set.FNA
#locus='LOC00263'
for locus in $(cat loci_list.txt);
do
  #echo $locus
  nb_TTS=$(seqkit grep -w0 -nrp -$locus ../FULL_REFERENCE_nucl.FNA | seqkit seq -n | grep 'scaffold' | wc -l) # nb transcriptomes TS
  if [ $nb_TTS -gt 5 ]; then # if number of transcript (from 1KP) template sequences is above 4 just select a few (ie the matching ones)
    # retrieve the matching TS(s) for the locus (SWGX WWQZ sequences transcripts and/or Angio353 seqs)
    #cat matching_subjects/*-$locus.txt
    seqkit grep -w0 -nrf matching_subjects/*-$locus.txt ../FULL_REFERENCE_nucl.FNA >> final_sequences_set.FNA
    # retrieve the corresponding Angio353 sequences for the locus
    # seqkit grep -w0 -nrp $locus ../FULL_REFERENCE_nucl.FNA | seqkit seq -n | grep  'scaffold'
    seqkit grep -w0 -nrp -$locus ../FULL_REFERENCE_nucl.FNA | seqkit seq -n | grep -v 'scaffold' > $locus"_TS_to_get.txt"
    if [ -s $locus"_TS_to_get.txt" ]; then
      seqkit grep -w0 -nrf $locus"_TS_to_get.txt" ../FULL_REFERENCE_nucl.FNA >> final_sequences_set.FNA
    fi
    rm $locus"_TS_to_get.txt"
  else
    seqkit grep -w0 -nrp -$locus ../FULL_REFERENCE_nucl.FNA >> final_sequences_set.FNA
  fi
done
echo "... done retrieving the sequences..."
echo "Done creating a sequence set with the loci to which the cleaned hits matched..."

###################################################
#### 5 Append the Angio353 and remove duplicates
###################################################
cd $path_to_tmp/cleaned_sequences_set
echo 'Start recovering the Angio353 and remove duplicates...'
seqkit seq -n ../Angio353_from_Melasto689.FNA | sed 's/.*-//g' | sort | uniq > Angio353_from_Melasto689_list.txt
#seqkit seq -n final_sequences_set.FNA | grep -o -e '[[:punct:]][[:digit:]]\{4\}\($\|[[:punct:]]\)' | grep -o '[[:digit:]]\{4\}'| sort | uniq > Angio353_in_final_set.txt
#grep -v -f Angio353_in_final_set.txt Angio353_from_Melasto689_list.txt > Angio353_to_add.txt

locus="5162"
for locus in $(cat Angio353_from_Melasto689_list.txt);
do
  echo $locus
  seqkit grep -w0 -nrp -$locus ../FULL_REFERENCE_nucl.FNA >> final_sequences_set.FNA
  seqkit grep -w0 -nrp "__"$locus"__" ../FULL_REFERENCE_nucl.FNA >> final_sequences_set.FNA
done

seqkit rmdup -w0 -D name_duplicates.txt -n final_sequences_set.FNA > final_sequences_set_rmdups.FNA
seqkit rmdup -w0 -D sequence_duplicates.txt -s final_sequences_set_rmdups.FNA > final_sequences_set.FNA

echo 'Done recovering the Angio353 and remove duplicates...'

###################################################
#### 6 Finalisation
###################################################
cd $path_to_tmp/cleaned_sequences_set
echo "Start finalization"
# sort by locusCustomID
cat final_sequences_set.FNA | seqkit rmdup -n | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' | seqkit sort -w0 -n | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > final_sequences_set.FNA.temp
#cat final_sequences_set.FNA | seqkit rmdup -n | sed 's/-ASGPB/_ASGPB/g' | sed 's/-Genoscope/_Genoscope/g' | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' | seqkit sort -w0 -n | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > final_sequences_set.FNA.temp
#cat final_sequences_set.FNA | seqkit seq -n | sort | uniq -d
rm final_sequences_set.FNA
#rename -v '.FNA.temp' '.FNA' ./*

echo "... renaming the template sequences..."
seqkit replace -w0 -p 'hits_consensus' -r 'hc' final_sequences_set.FNA.temp > final_sequences_set.FNA.temp2
seqkit replace -w0 -p 'transcriptome_combo' -r 'tr_combo' final_sequences_set.FNA.temp2 > final_sequences_set.FNA
#sed -i 's/__iTAG/_iTAG/g' final_sequences_set.FNA

cat final_sequences_set.FNA | seqkit seq -n > TS_full_list.txt

#seqkit translate -w0 final_sequences_set.FNA > final_sequences_set.FAA
rm -r *temp*
echo "Done finalization"


######################################################################################################
#### 6bis intermediate transfer
# copy everything to output folder
scp -rp $path_to_tmp/* $path_to_dir_out/
######################################################################################################

####################################################
##### 7 Self-Blast
####################################################

echo "Self-Blasting..."
echo "Starting building BLAST databases...";
cd $path_to_tmp
mkdir BLAST_DB
makeblastdb -in cleaned_sequences_set/final_sequences_set.FNA -out BLAST_DB/final_sequences_set_nucl -dbtype nucl -parse_seqids
makeblastdb -in cleaned_sequences_set/final_sequences_set.FAA -out BLAST_DB/final_sequences_set_prot -dbtype prot -parse_seqids
echo "Done building BLAST databases.";

blastp -task blastp-fast -db BLAST_DB/final_sequences_set_prot -query cleaned_sequences_set/final_sequences_set.FAA -evalue 1e-06 -outfmt "10 std qframe qcovs qseq" -num_threads $SLURM_CPUS_PER_TASK > selfblast_final_sequences_set_blastp_out_SB.csv
cat selfblast_final_sequences_set_blastp_out_SB.csv | awk -F ',' 'OFS = "," {
temp1 = $1
temp2 =$2
sub(".+-","", temp1) 
sub(".+-","", temp2)
if (temp1 != temp2)
print $0
}' > different_query_and_subject_loci_SB.csv
cat different_query_and_subject_loci_SB.csv | csvtk filter -H -f'3>60' > filtered_different_query_and_subject_loci_SB.csv
echo "Done Self-Blast."

####################################################
##### 8 Remove transcriptome sequences that are supersets of mega353 sequences
####################################################
cat selfblast_final_sequences_set_blastp_out_SB.csv | csvtk filter -H -f'14=100' | csvtk filter -H -f'3=100' | awk -F ',' 'OFS = "," {
if ($1 != $2)
print $0
}' | awk -F, '$2 ~ /scaffold/ {print}' | awk -F, '$1 ~ /^WWQZ|^SWGX/ {print}' > different_query_and_subject_SB.csv
cat different_query_and_subject_SB.csv

cat different_query_and_subject_SB.csv | awk -F ',' '{print $2}' > TS_to_remove.txt
cd $path_to_tmp/cleaned_sequences_set
seqkit grep -w0 -nrvf ../TS_to_remove.txt final_sequences_set.FNA > final_sequences_set.FNA.temp
rm final_sequences_set.FNA.
mv final_sequences_set.FNA.temp final_sequences_set.FNA

##################################################
#### 9 Transfer & clean
##################################################

echo "Transfert data node -> blue storage";

# make output folder in home directory
#mkdir $path_to_dir_out

# copy everything to output folder
scp -rp $path_to_tmp/* $path_to_dir_out/

echo "done moving"
echo "FINISHED"

# All the data are automatically deleted on the node by SLURM
