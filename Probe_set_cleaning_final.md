Cleaning of the Melastomataceae probe set for target enrichment
================

The basic idea of the cleaning process is to build a set of reference
sequences and to align every template sequence from Melasto689 to this
set of references, correcting the frame, reverse-complementing if
needed, and removing stop codons and poorly aligned regions. The cleaned
template sequences are then associated to the reference sequences they
matched into a final sequence set. Finally, the final sequence set go
through different checks and fine-tuned cleaning.

The process will be divided in 3 main steps:

-   **Step 1**: build the set of reference sequences  
-   **Step 2**: align and associate the Melasto689 template sequences to
    the set of reference into a final sequence set  
-   **Step 3**: final checks and fine-tuned cleaning: extra steps to
    make sure the sequences are clean

**Directory organization.** For simplicity, all the following assumes
analyses are run in a single directory `$path_to_wd`. This directory is
divided in 6 sub-directory:

-   [`prepare1KP`](prepare1KP) (Step 1: 1. - Prepare the 1KP
    transcriptomes reference sequences set)
-   [`prepare_mega353`](prepare_mega353) (Step 1: 2. - Prepare the
    mega353 reference sequences set)
-   [`full_reference_DB`](full_reference_DB) (Step 1: 3. - Merge the two
    previous reference sequences set)
-   [`original_template_sequences`](original_template_sequences) (Step
    2: 1. - Prepare the original Melast689 template sequence set)
-   [`melasto689_on_full_ref`](melasto689_on_full_ref) (Step 2: 2. -
    Align and associate the Melasto689 template sequences to the set of
    reference)
-   [`CLEAN_PROBE_SET`](CLEAN_PROBE_SET) (Step 3 - Contains the final
    and cleaned probe set)

``` bash
cd $path_to_wd
mkdir prepare1KP
mkdir prepare_mega353
mkdir full_reference_DB
mkdir original_template_sequences
mkdir melasto689_on_full_ref
mkdir CLEAN_PROBE_SET
```

**Naming convention.** Sequences names follow the convention used in
many probe sets for targeted sequencing and associated programs
(e.g.??Angio353, HybPiper, Captus, SECAPR). The name of a sequence is
composed of the template sequence ID (generally identifying the taxon
from which the sequence initially come from, sometimes with additional
info) and the locus name separated by an hyphen,
e.g.??`>template_sequence_ID-locusID`

**Scripts.** Every step of the cleaning process is presented and
explained in detail below. The command are presented as if run in
interactive mode. For some steps, the running time is long and was
actually not achieved in interactive mode, but using `.sh` scripts
(scripts were run on a computing cluster that uses the SLURM manager,
`sbatch` command). The `.sh` scripts are made available in this
repository (see the [`script`](scripts) folder) and a statement is
present at the beginning of each step specifying which script to refer
to.

**Files.** Note that not all the intermediate files produced by the
cleaning process are present in this repository. The missing files are
some intermediate .FNA files or heavy (\>100 Mo) blast output tables.
The most important files are still included.

**Dependencies**

-   [SeqKit 2.0.0](https://bioinf.shenwei.me/seqkit/)
-   [NCBI BLAST 2.10.1](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
-   [csvtk 0.23.0](https://bioinf.shenwei.me/csvtk/)
-   [bedtools 2.30.0](https://bedtools.readthedocs.io/en/latest/)
-   [R 4.2](https://www.r-project.org/) (needs packages
    [DECIPHER](http://www2.decipher.codes/),
    [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html))
-   [Python 3.10](https://www.python.org/downloads/release/python-3100/)
-   [HybPiper 2.0.1](https://github.com/mossmatters/HybPiper)
-   [AliView](https://ormbunkar.se/aliview/)
-   [Muscle v3.8.425](https://drive5.com/muscle/) (as wrapped in
    AliView)

**Glossary**:

-   **Angio353** (sometimes referred as **full Angio353**): the original
    Angiosperms 353 probe set
-   **locusCustomID**: temporary locus identifier given to the loci
    during the cleaning process
-   **mega353**: the Angiosperms 353 probe set extended with
    transcriptome sequences (see:
    <https://github.com/chrisjackson-pellicle/NewTargets>)
-   **Melasto689**: the full Melastomataceae probe set. It is composed
    of 384 loci, each of them being composed of 1 to 4 template
    sequences, for a total of 689 template sequences. It is divided into
    4 subsets of loci:
    -   **angio353** (with small **a**): a set of 266 loci (439 template
        sequences) that come from the Angiosperms 353 probe set
    -   **markerminer**: a set of 104 loci (215 template sequences) that
        come from putative single-copy loci identified with MarkerMiner
    -   **functional_extra**: a set of 6 loci (12 template sequences)
        that come from functional genes
    -   **RM2016**: a set of 8 loci (23 template sequences) that come
        from [Reginato & Michelangeli
        (2016)](https://doi.org/10.3732/apps.1500092)
-   **reference sequence** (or **RS**): a sequence from the reference
    set; in the sequence set, a same locus can have several template
    sequences
-   **template sequence** (or **TS**): a sequence from the Melasto689
    set or final sequence set; in the sequence set, a same locus can
    have several template sequences

# Step 1: build the set of reference sequences

The set of reference is composed of
[mega353](https://doi.org/10.1002/aps3.11420) (that is an extended
version of [Angio353](https://doi.org/10.1093/sysbio/syy086)) and the
transcriptomes of two Melastomataceae species from the [1KP
project](https://sites.google.com/a/ualberta.ca/onekp/?pli=1).

The set of reference sequence will contain different sequences. The
sequences will be grouped under different ???loci??? based on their
similarity. During the process of building the set of references, a
temporary locus name (sometimes referred as `*locusCustomID*`) will be
assigned to the loci in the form `LOC#####` (where `#` represents a
digit, e.g.??`LOC12345`). This `*locusCustomID*` will be used as locus
name for most of the sequences during the whole cleaning process, but
the loci names will be reverted to their ???old??? names (i.e.??the names
used in Melasto689 and Angio353) in Step 3: 3.

## 1. Prepare the 1KP transcriptomes reference database

> The `sbatch` script corresponding to this section is:
> [`prepare_1KP_DB.sh`](scripts/prepare_1KP_DB.sh)

We will use the transcriptomes of the two Melastomataceae samples
present in the 1KP project dataset: SWGX and WWQZ (*Tetrazygia bicolor*
and *Medinilla magnifica*, respectively).

### 1.1 Retrieve the sequences

Download SWGX and WWQZ data from the [1KP project (updated
release)](https://dx.doi.org/10.5524/100910), and decompress (unzip) the
files into `$path_to_wd/prepare1KP`. We will use of the
`XXXX-translated-protein.fa` and `XXXX-translated-nucleotides.fa` fasta
files.

### 1.2 Group the sequences according to their ortholog group

Sequences in the SWGX and WWQZ data have a unique identifier (referred
as `transcriptSequencesID`), in the form
`>scaffold-XXXX-2#######-Genus_species`, where XXXX is either SWGX or
WWQZ, \# represents a digit, and Genus_species is either
Tetrazygia_bicolor (SWGX) or Medinilla_magnifica (WWQZ),
e.g.??`>scaffold-WWQZ-2000004-Medinilla_magnifica`.

The sequences need to be grouped according to their ortholog group. To
do so, we???ll use the 2<sup>nd</sup> column in the
`XXXX-translated-reference-names.tsv` files, containing the gene
identifiers according to OrthoFinder (see [One Thousand Plant
Transcriptomes Initiative
(2019)](https://doi.org/10.1038/s41586-019-1693-2), referred as
`locusOrthologID`. The identifier is in the form
???gi\|255582030\|ref\|XP_002531812.1\|???.

-   get the list of the unique orthogroup identifiers from the
    `SWGX-translated-reference-names.tsv`
    `WWQZ-translated-reference-names.tsv` files into the file
    `orthoIDs_list.txt`

``` bash
cd $path_to_wd/prepare1KP
cat SWGX-translated-reference-names.tsv WWQZ-translated-reference-names.tsv | awk '{print $2}' | sort | uniq > orthoIDs_list.txt # get the list of the unique gene identifiers according to OrthoFinder
```

-   append a unique custom locus name (*`locusCustomID`*) in the form
    LOC##### to every orthogroup (*`locusOrthologID`*) and store them in
    the array `$loc_list`. The array contains locus names in the form
    `locusOrthologID__locusCustomID`.

``` bash
loc_list=$(cat orthoIDs_list.txt | awk '{ printf $1"__LOC" "%05d\n", ++a }') # append a unique custom locus name in the form LOC##### (locusCustomID) to every orthogroup (locusOrthologID) and store them in the array loc_list
```

-   create a folder called `orthogroups` and initialize the files:
    -   `ortho_loci_table.txt`: will contain the correspondance between
        *`locusCustomID`*, *`locusOrthologID`* and
        *`transcriptSequencesID`*
    -   `SWGX-WWQZ-translated-nucleotides-grouped.FNA`: will contain the
        transcriptome sequences, grouped by custom locus ID
    -   `SWGX-WWQZ-translated-prot-grouped.FAA`: will contain the
        transcriptome sequences, grouped by custom locus ID

``` bash
mkdir orthogroups
echo locusCustomID locusOrthologID transcriptSequencesID > orthogroups/ortho_loci_table.txt # initialize file
:> orthogroups/SWGX-WWQZ-translated-nucleotides-grouped.FNA # initialize file
:> orthogroups/SWGX-WWQZ-translated-prot-grouped.FAA # initialize file
```

-   for every locus name in the list of loci `$loc_list`:
    -   define \$locID as the custom locus ID (*`locusCustomID`*)
    -   define \$orthoID as the orthogroup ID (*`locusOrthologID`*)
    -   retrieve the `transcriptSequencesID`(s) corresponding to the
        *`locusOrthologID`* in SWGX and WWQZ databases
    -   if `transcriptSequencesID`(s) were retrieved:
        -   extract the sequence(s) from the databases, append the
            custom locus ID (*`locusCustomID`*) to its name, and append
            the renamed sequence to
            `SWGX-WWQZ-translated-nucleotides-grouped.FNA` and
            `SWGX-WWQZ-translated-prot-grouped.FAA`
    -   populate the `ortho_loci_table.txt` with the corresponding
        *`locusCustomID`*, *`locusOrthologID`* and
        *`transcriptSequencesID`*

``` bash
for name in $loc_list;
do
  locID=$(echo $name | sed 's/.*__//g')
  orthoID=$(echo $name | sed 's/__.*//g')

  # retrieve the sequence ID(s) corresponding to the locus in SWGX and WWQZ databases
  grep -e $orthoID SWGX-translated-reference-names.tsv | awk '{print $1}' > SWGX_seqs
  grep -e $orthoID WWQZ-translated-reference-names.tsv | awk '{print $1}' > WWQZ_seqs
  
  if [ -s SWGX_seqs ]; then
    seqkit grep -w0 -n -f SWGX_seqs SWGX-translated-nucleotides.fa | sed 's/-/_/g' | sed '/^>/s/$/-'$locID'/g' >> orthogroups/SWGX-WWQZ-translated-nucleotides-grouped.FNA
    seqkit grep -w0 -n -f SWGX_seqs SWGX-translated-protein.fa | sed 's/-/_/g' | sed '/^>/s/$/-'$locID'/g' >> orthogroups/SWGX-WWQZ-translated-prot-grouped.FAA
  fi
  if [ -s WWQZ_seqs ]; then
    seqkit grep -w0 -n -f WWQZ_seqs WWQZ-translated-nucleotides.fa | sed 's/-/_/g' | sed '/^>/s/$/-'$locID'/g' >> orthogroups/SWGX-WWQZ-translated-nucleotides-grouped.FNA
    seqkit grep -w0 -n -f WWQZ_seqs WWQZ-translated-protein.fa | sed 's/-/_/g' | sed '/^>/s/$/-'$locID'/g' >> orthogroups/SWGX-WWQZ-translated-prot-grouped.FAA
  fi
  
  echo $locID $orthoID $(cat SWGX_seqs) $(cat WWQZ_seqs) >> orthogroups/ortho_loci_table.txt
  rm SWGX_seqs WWQZ_seqs
done
```

-   remove sequences less than 30 amino-acids (AA) and less than 90
    nucleotides (NT) long

``` bash
cd orthogroups
seqkit seq -w0 -m 90 SWGX-WWQZ-translated-nucleotides-grouped.FNA > SWGX-WWQZ-translated-nucleotides-grouped.FNA.temp
seqkit seq -w0 -m 30 SWGX-WWQZ-translated-prot-grouped.FAA > SWGX-WWQZ-translated-prot-grouped.FAA.temp
rm -r SWGX-WWQZ-translated-prot-grouped.FAA SWGX-WWQZ-translated-nucleotides-grouped.FNA
rename -v '.FAA.temp' '.FAA' ./*
rename -v '.FNA.temp' '.FNA' ./*
cd $path_to_wd/prepare1KP
```

### 1.2 Check for sequences with high similarity

Sequences with high similarity need to be grouped into a same ???locus???,
otherwise the template sequences from Melasto689 will align to several
loci which is not desirable since we want to distribute these template
sequences into a single loci.  
To group sequences with high similarity, we will run a ???self-Blast???
using `blastp-fast` (protein to protein alignment).

#### 1.2.1 Build BLAST databases

Make a custom BLAST database with the 2 transcriptomes (SWGX and WWQZ)
concatenated and grouped by orthogroup. Build a database for both
amino-acids and nucleotides sequences.

``` bash
cd $path_to_wd/prepare1KP
mkdir BLAST_custom_DB
makeblastdb -in orthogroups/SWGX-WWQZ-translated-nucleotides-grouped.FNA -out BLAST_custom_DB/SWGX-WWQZ_nucleotides_grouped_BLAST_DB -dbtype nucl -parse_seqids # actually not necessary
makeblastdb -in orthogroups/SWGX-WWQZ-translated-prot-grouped.FAA -out BLAST_custom_DB/SWGX-WWQZ_protein_grouped_BLAST_DB -dbtype prot -parse_seqids
```

#### 1.2.2 Run the self-BLAST

``` bash
cd $path_to_wd/prepare1KP
blastp -task blastp-fast -db BLAST_custom_DB/SWGX-WWQZ_protein_grouped_BLAST_DB -query orthogroups/SWGX-WWQZ-translated-prot-grouped.FAA -evalue 1e-06 -outfmt "10 std qframe qcovs qseq" -num_threads $SLURM_CPUS_PER_TASK > selfblast_SWGX-WWQZ_protein_grouped_blastp_out_SB1.csv
```

The raw result of a self-BLAST needs to be filtered since every sequence
will inevitably have a 100% match with itself. So, we filter out all the
matches between two sequences belonging to the same locus. Also, we only
keep the matches that have either 1) a percentage of identity greater
than 60%, an alignment length greater than 50 AA, a bitscore greater
than 100 and a query coverage greater than 50% (meaning at least 50% of
the query covers the subject sequence), or 2) a percentage of identity
greater than 80%, an alignment length greater than 50 AA and a bitscore
greater than 50.

``` bash
cat selfblast_SWGX-WWQZ_protein_grouped_blastp_out_SB1.csv | awk -F ',' 'OFS = "," {
temp1 = $1
temp2 =$2
sub(".+-","", temp1) 
sub(".+-","", temp2)
if (temp1 != temp2)
print $0
}' > different_query_and_subject_loci_SB1.csv

cat different_query_and_subject_loci_SB1.csv | csvtk filter2 -H -f '$3 > 60 && $4 > 50 && $14 > 50 && $12 > 100 || $3 > 80 && $4 > 50 && $12 > 50' > filtered_different_query_and_subject_loci_SB1.csv
```

Some of the sequences with high similarity were assigned to different
orthogroups. We will thus group these sequences, which mean grouping
some supposedly ortholog loci.

### 1.3 Group loci with high sequence similarity

#### 1.3.1 Build a similarity table

Retrieve the pairs of loci that need to be grouped.

``` bash
cat filtered_different_query_and_subject_loci_SB1.csv | awk -F ',' '{
query = $1
subject = $2
sub(".+-","", query) 
sub(".+-","", subject)
print query"," subject
}' > pairs_SB1.csv
```

A locus can have high sequence similarity with more than one other
locus, we need to retrieve all the loci associations, and print them in
a similarity table (`selfblast_table.txt`).

For every locus in the list of loci from the pairs of loci:

-   retrieve the other locus (loci) it matches to, and arrange them into
    a list
-   search if any of the locus in the list is already in the table
    `selfblast_table.txt`
-   if none of the locus is in the similarity table, then append the
    list at the end of the table
-   if any of the locus is found in the similarity table, then extract
    the matching list of loci in the similarity table, append the list
    to the extracted list, remove duplicates loci and replace the
    original list with the newly generated list

``` bash
query_list=$(cat pairs_SB1.csv | awk -F ',' '{print $1}' | sort | uniq)
:> selfblast_table.txt
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
cat selfblast_table.txt > selfblast_table_fromSB1.txt #intermediate backup

#check for no duplicates:
echo Total number of IDs: $(cat selfblast_table.txt | sed 's/,/\n/g' | wc -l)
echo "Number of IDs (dups removed):" $(cat selfblast_table.txt | sed 's/,/\n/g' | sort | uniq | wc -l)
```

Each line of the table `selfblast_table.txt` contains a list of
comma-separated matching loci in the form:
*`locusCustomID, correspondingLocus1, correspondingLocus2, ...`*. The
name retained as *`locusCustomID`* is the one having the lowest suffix
number among the list of corresponding loci on a same line.

#### 1.3.2 Associate the loci in the fasta files

Now the list of corresponding loci based on similarity has been
established, we need to associate them in the `.FNA` and `.FAA` files.

Transform the table into a list of pairs of
*`locusCustomID,correspondingLocus`* (`selfblast_table_paired.txt`).

``` bash
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
```

Transform the comma-separated list into a tab-separated list
(`replacement.txt`), that can be read by `seqkit replace` command.

``` bash
cat selfblast_table_paired.txt | awk -F ',' '{print $2"\t"$1}' > replacement.txt
```

Use `seqkit replace` to replace every correspondingLocus to its matching
*`locusCustomID`*.

``` bash
seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key orthogroups/SWGX-WWQZ-translated-prot-grouped.FAA > orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA
seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key orthogroups/SWGX-WWQZ-translated-nucleotides-grouped.FNA > orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA
```

Finally sort the fasta files by *`locusCustomID`*.

``` bash
sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA.temp
sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA.temp
rm orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA
rename -v '.FAA.temp' '.FAA' ./orthogroups/*
rename -v '.FNA.temp' '.FNA' ./orthogroups/*
```

### 1.4 Final self-BLAST check

#### 1.4.1 Build BLAST database

``` bash
cd $path_to_wd/prepare1KP
makeblastdb -in orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA -out BLAST_custom_DB/SWGX-WWQZ_nucleotides_grouped_SB1_BLAST_DB -dbtype nucl -parse_seqids
makeblastdb -in orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA -out BLAST_custom_DB/SWGX-WWQZ_protein_grouped_SB1_BLAST_DB -dbtype prot -parse_seqids
```

#### 1.4.2 Run the self-BLAST check

``` bash
blastp -task blastp-fast -db BLAST_custom_DB/SWGX-WWQZ_protein_grouped_SB1_BLAST_DB -query orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA -evalue 1e-06 -outfmt "10 std qframe qcovs qseq" > check_selfblast_SWGX-WWQZ_protein_grouped_blastp_out.csv
cat check_selfblast_SWGX-WWQZ_protein_grouped_blastp_out.csv | awk -F ',' 'OFS = "," {
temp1 = $1
temp2 =$2
sub(".+-","", temp1)
sub(".+-","", temp2)
if (temp1 != temp2)
print $0
}' > check_different_query_and_subject_loci.csv
cat check_different_query_and_subject_loci.csv | csvtk filter -H -f'3>60' | csvtk filter -H -f'4>50' | csvtk filter -H -f'14>50' | csvtk filter2 -H -f '$12 > 100' > check_filtered_different_query_and_subject_loci.csv
```

The file `check_filtered_different_query_and_subject_loci.csv` resulting
from the check self-BLAST is empty, meaning no sequence from a locus has
similarity with a sequence of any other locus.

## 2. Prepare the mega353 reference database

The mega353 is an expansion of the Angio353 probe set, incorporating
sequences from more than 500 transcriptomes from the the 1KP project.
Along with the dataset, a script is provided in order to filter the
sequences according to the group of interest. Here we will filter the
mega353 dataset to keep the sequences from the samples in the Myrtales
order.

Download the NewTargets directory from
<https://github.com/chrisjackson-pellicle/NewTargets> and unzip the
files. Store the `NewTargets` folder into `$path_to_wd/prepare_mega353`.

A filtering text file is created, in order to filter on the Myrtales
order
([`select_file_Myrtales.txt`](prepare_mega353/select_file_Myrtales.txt)).

The `filter_mega353.py` script is then run.

``` bash
cd $path_to_wd/prepare_mega353/NewTargets
python filter_mega353.py mega353.fasta select_file_Myrtales.txt -filtered_target_file ../mega353_Myrtales.fa
```

The mega353 reference sequences set used is the
[`mega353_Myrtales.fa`](prepare_mega353/mega353_Myrtales.fa) file.

## 3. Associate the transcriptome and mega353 references databases into a single set

> The `sbatch` script corresponding to this section is:
> [`prepare_full_reference_DB.sh`](scripts/prepare_full_reference_DB.sh)

We will now associate the the transcriptome and mega353 references set
into a single set of references. The association will be done based on
sequence similarity.

### 3.1 Prepare the paths and data

Define the paths to previously prepared databases.

``` bash
cd $path_to_wd/full_reference_DB
path_to_onekp=$path_to_wd/prepare1KP;
Angio353_FNA=$path_to_wd/prepare_mega353/mega353_Myrtales.fa;
```

Concatenate the transcriptome reference databases in the current
directory into the `FULL_REFERENCE_nucl.FNA` and
`FULL_REFERENCE_prot.FAA` files. Sequences from the mega353 database
will be appended to these files.

``` bash
cat $path_to_wd/prepare1KP/orthogroups/SWGX-WWQZ-translated-nucleotides-grouped_SB1.FNA > FULL_REFERENCE_nucl.FNA
cat $path_to_wd/prepare1KP/orthogroups/SWGX-WWQZ-translated-prot-grouped_SB1.FAA > FULL_REFERENCE_prot.FAA
nuc_FNA="FULL_REFERENCE_nucl.FNA"
prot_FAA="FULL_REFERENCE_prot.FAA"
```

Copy the mega353 database in the current directory and rename the
sequences to remove hyphens that do not separate sequence ID from locus
name.

``` bash
dos2unix $Angio353_FNA
cat $Angio353_FNA | sed 's/ .*//g' | sed 's/-[[:digit:]]\{4\}_exonerate_grafted/_EG/g' | sed 's/grafted_with/GW/g' | sed 's/-v/_v/g' | sed 's/\./_/g' | sed 's/-ASGPB/_ASGPB/g' | sed 's/-iTAG/_iTAG/g' | sed 's/-Genoscope/_Genoscope/g' > Angio353.FNA
Angio353_FNA=$path_to_tmp/Angio353.FNA
```

Define the path to the transcriptome protein BLAST database from Step 1:
1.4.1.

``` bash
prot_DB=$path_to_wd/prepare1KP/BLAST_custom_DB/SWGX-WWQZ_protein_grouped_SB1_BLAST_DB
```

### 3.2 Associate mega353 loci with SWGX and WWQZ transcriptomes

Extract the list of mega353 loci in a text file. Initialize matching
table between mega353 loci and custom locus ID (*`locusCustomID`*).

``` bash
cd $path_to_wd/full_reference_DB
seqkit seq -n $Angio353_FNA | sed 's/.*-//g' | sort | uniq > Angio353_list.txt
:> Angio353_to_custom_loci_table.txt
```

For each locus in the list of mega353 loci:

-   create a sub-directory and set it as current directory
-   extract the locus sequence from the mega353 database into a .FNA
    file, translate it (.FAA)
-   search for similar protein sequence into the transcriptome protein
    database using `blastp` and store the results in a .csv file
-   if the .csv file exists and is not empty (i.e.??blast found matching
    sequences), then:
    -   filter the .csv to keep only query sequences with percentage
        identity greater than 60% to the matching subject, alignment
        length greater than 50 amino-acids, bitscore greater than 100,
        frame greater than 0 (i.e.??sequences not in reverse complement)
        and a query coverage greater than 50% (meaning at least 50% of
        the query covers the subject sequence)
-   if the filtered .csv file exists and is not empty (i.e.??matching
    sequences remain after filtering), then:
    -   extract the *`customLocusID`*(s) from the matching sequence(s)
    -   if more than 1 *`customLocusID`* were found, then:
        -   print a warning
        -   append the locus sequences to the full reference database as
            they are (no renaming)
        -   print the mega353 locus name and the mention
            ???multiplematches??? to the matching table
    -   else (i.e.??only 1 *`customLocusID`* was found), then:
        -   rename the locus sequence so that the mega353 locus name
            (4-digit sequence) is enclosed in 2 \_\_ (double underscore)
            and moved to the sequence ID, and so that the
            *`customLocusID`* is appended at the end of the sequence
        -   append the renamed locus sequence to the full reference
            database
        -   print the mega353 locus name and its *`customLocusID`* to
            the matching table
-   else (i.e.??no matching sequences remain after filtering), then:
    -   append the locus sequences to the full reference database as
        they are (no renaming)
    -   print the mega353 locus name and the mention ???nomatch??? to the
        matching table
-   set the upper level directory as current directory

``` bash
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
```

### 3.4 Check for sequence similarity

This step checks sequence similarity for sequences from different loci.
If highly similar sequences from different loci are found, the loci have
to be grouped under the same name. To check for sequence similarity, we
run a ???self-Blast??? using `blastp-fast` (protein to protein alignment).

#### 3.4.1 Build BLAST database

``` bash
mkdir FULL_REFERENCE_BLAST_DB
makeblastdb -in FULL_REFERENCE_nucl.FNA -out FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_nucl_BLAST_DB -dbtype nucl -parse_seqids
makeblastdb -in FULL_REFERENCE_prot.FAA -out FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_prot_BLAST_DB -dbtype prot -parse_seqids
```

#### 3.4.2 Run the self-BLAST

``` bash
blastp -task blastp-fast -db FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_prot_BLAST_DB -query FULL_REFERENCE_prot.FAA -evalue 1e-06 -outfmt "10 std qframe qcovs qseq" -num_threads $SLURM_CPUS_PER_TASK > full_ref_blastp_out.csv
```

The raw result of a self-BLAST needs to be filtered since every sequence
will inevitably have a 100% match with itself. So, we filter out all the
matches between two sequences belonging to the same locus. Also, we only
keep the query sequences with percentage identity greater than 60% to
the matching subject, alignment length greater than 50 amino acids,
bitscore greater than 100, frame greater than 0 (i.e.??sequences not in
reverse complement) and a query coverage greater than 50% (meaning at
least 50% of the query covers the subject sequence).

``` bash
cat full_ref_blastp_out.csv | awk -F ',' 'OFS = "," {
temp1 = $1
temp2 =$2
sub(".+-","", temp1) 
sub(".+-","", temp2)
if (temp1 != temp2)
print $0
}' > different_query_and_subject_loci.csv
cat different_query_and_subject_loci.csv | csvtk filter2 -H -f '$3 > 60 && $4 > 50 && $14 > 50 && $12 > 100 && $13 > 0' > filtered_different_query_and_subject_loci.csv
```

### 3.5 Associate the loci having sequences with high similarity

Highly similar sequences from different loci were found, these loci have
to be grouped under the same name.

#### 3.5.1 Build a similarity table

Retrieve the pairs of loci that need to be grouped. The file contains a
list of pairs of query-subject.

``` bash
cat filtered_different_query_and_subject_loci.csv | awk -F ',' '{
query = $1
subject = $2
sub(".+-","", query) 
sub(".+-","", subject)
print query"," subject
}' | sort | uniq > pairs.csv
```

A same locus may have to be associated with more than one other locus.
We thus need to retrieve all the loci associations, and print them in a
similarity table (`selfblast_table.txt`).

For every query list of the pairs of loci: - retrieve the other locus
(loci) it matches to (subject), and arrange them into a list - search if
any of the locus in the list is already in the table
`selfblast_table.txt` - if none of the locus is in the similarity table,
then append the list at the end of the table - if any of the locus is
found in the similarity table, then extract the matching list of loci in
the similarity table, append the list to the extracted list, remove
duplicates loci and replace the original list with the newly generated
list

``` bash
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
```

#### 3.5.2 Associate the loci in the fasta files

Now the list of corresponding loci based on similarity has been
established (`selfblast_table.txt`), we need to associate them in the
`.FNA` and `.FAA` files.

Transform the table into a list of pairs of
*`locusCustomID,correspondingLocus`* (`selfblast_table_paired.txt`).

``` bash
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
```

Transform the comma-separated list into a tab-separated list
(`replacement.txt`), that can be read by `seqkit replace` command.

``` bash
:> replacement.txt
cat selfblast_table_paired.txt | grep ',[[:digit:]]\{4\}$' | awk -F ',' '{print $2"\t__"$2"__-"$1}' >> replacement.txt
cat selfblast_table_paired.txt | grep -v ',[[:digit:]]\{4\}$' | awk -F ',' '{print $2"\t"$1}' >> replacement.txt
```

Use `seqkit replace` to replace every *`correspondingLocus`* to its
matching *`locusCustomID`*.

``` bash
seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key FULL_REFERENCE_prot.FAA | seqkit replace -w0 -p '-_' -r '_' > FULL_REFERENCE_prot_SB1.FAA
seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key FULL_REFERENCE_nucl.FNA | seqkit replace -w0 -p '-_' -r '_' > FULL_REFERENCE_nucl_SB1.FNA
```

Finally sort the fasta files by *`locusCustomID`*.

``` bash
sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' FULL_REFERENCE_prot_SB1.FAA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > FULL_REFERENCE_prot_SB1.FAA.temp
sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' FULL_REFERENCE_nucl_SB1.FNA | seqkit sort -w0 | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > FULL_REFERENCE_nucl_SB1.FNA.temp
rm FULL_REFERENCE_prot_SB1.FAA FULL_REFERENCE_nucl_SB1.FNA
rename -v '.FAA.temp' '.FAA' ./*
rename -v '.FNA.temp' '.FNA' ./*
```

### 3.6 Final self-BLAST check

#### 3.6.1 Build BLAST database

``` bash
cd $path_to_wd/full_reference_DB
makeblastdb -in FULL_REFERENCE_nucl_SB1.FNA -out FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_nucl_BLAST_DB_FINAL -dbtype nucl -parse_seqids
makeblastdb -in FULL_REFERENCE_prot_SB1.FAA -out FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_prot_BLAST_DB_FINAL -dbtype prot -parse_seqids
```

#### 3.6.2 Run the self-BLAST check

``` bash
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
```

The file `check_filtered_different_query_and_subject_loci.csv` resulting
from the check self-BLAST is empty, meaning no sequence from a locus has
similarity with a sequence of any other locus.

**The full set of reference sequences is ready**
in`$path_to_wd/full_reference_DB`:
[`FULL_REFERENCE_nucl_SB1.FNA`](full_reference_DB/FULL_REFERENCE_nucl_SB1.FNA.gz)
and
[`FULL_REFERENCE_prot_SB1.FAA`](full_reference_DB/FULL_REFERENCE_prot_SB1.FAA.gz)

# Step 2: align the Melasto689 template sequences to the set of reference

## 1. Prepare the raw target sequences set (Melasto689)

Download or copy the [original Melastomataceae 689 target sequence
set](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8931zcrm2),
and decompress the `Template_sequences.zip` folder.

-   De-align the sequences
-   Remove all the ? characters

``` bash
cd $path_to_wd/original_template_sequences
cd Template_sequences
sed -i '/^>/!s/-//g' *.fasta
sed -i '/^>/!s/\?//g' *.fasta
```

### 1.1 Clean the KT######.1 loci

#### 1.1.1 KT377070.1 and KT377086.1

KT377070.1 and KT377086.1 are each split into 2 different files
(`KT377070.1_alignment.fasta` and `KT377070_alignment.fasta`, and
`KT377086.1_alignment.fasta` and `KT377086-Affzeli.fasta`,
respectively), supposedly representing 2 different ???versions??? of the
genes with high divergence (Jantzen, pers. comm.). Actually, the
alternative ???version??? is not highly divergent but in reverse complement
of what it should be. Here we just merge the fasta files in order to
have only 1 file per gene. We arbitrarily keep the ???.1??? suffix in the
names of the genes.

**KT377070.1**  
To avoid duplicates in fasta headers, we rename the sequence
???\>KT377070-Tibouchina??? (from `KT377070_alignment.fasta`) to
???\>KT377070.1-Tibouchina2???.

``` bash
cat KT377070.1_alignment.fasta KT377070_alignment.fasta | sed 's/KT377070-Tibouchina/KT377070-Tibouchina2/g' | sed 's/KT377070-/KT377070.1-/g' > KT377070.1_custom.fasta
rm KT377070.1_alignment.fasta KT377070_alignment.fasta
```

**KT377086.1**

``` bash
cat KT377086.1_alignment.fasta KT377086-Affzeli.fasta | sed 's/KT377086-/KT377086.1-/g' > KT377086.1_custom.fasta
rm KT377086.1_alignment.fasta KT377086-Affzeli.fasta
```

#### 1.1.2 Other KT######.1 loci

Modify the fasta headers of the template sequences for the KT######.1
loci so that all the headers have the ???.1??? suffix.

``` bash
cat KT377102.1_alignment.fasta | sed 's/KT377102-/KT377102.1-/g' > KT377102.1_custom.fasta
cat KT377110.1_alignment.fasta | sed 's/KT377110-/KT377110.1-/g' > KT377110.1_custom.fasta
rm KT377102.1_alignment.fasta KT377110.1_alignment.fasta
```

### 1.2 Rename the sequence names

Need to rename the sequences names to stick to the naming convention.

``` bash
ls -1 *.fasta | \
while read file
do
cat $file| sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > $file.FNA
done
rename '.fasta.FNA' '.FNA' *
rm *.fasta
```

### 1.3 Concatenate all the sequences into a single fasta file:

``` bash
cat *.FNA > ../Melastomataceae_689_original.FNA
```

## 2. Associate each template sequence (or portion of) from Melasto689 to a locus from the full reference set

> The `sbatch` script corresponding to this section is:
> [`melasto689_on_full_reference.sh`](scripts/melasto689_on_full_reference.sh)

### 2.1 Prepare the paths and data

Define the paths to previously prepared databases.

``` bash
cd $path_to_wd/melasto689_on_full_ref
path_to_Melasto689=$path_to_wd/original_template_sequences;
path_to_full_ref=$path_to_wd/full_reference_DB;
prot_DB=$path_to_full_ref/FULL_REFERENCE_BLAST_DB/FULL_REFERENCE_prot_BLAST_DB_FINAL
path_to_scripts=$path_to_scripts # define the path to the scripts
```

Copy the necessary files to the current directory.

``` bash
scp -r $path_to_Melasto689 ./
scp -r $path_to_full_ref/FULL_REFERENCE_nucl.FNA ./
```

### 2.2 Remove the Angio353 template sequences from Melasto689

The Melasto389 probe set includes template sequences from the Angio353
dataset (what we call angio353, with a small **a**, in the glossary).
They need to be removed beforehand. They are easily identifiable with
their name in the form `XXXX-####` where `X` is a letter and `#` a
digit.

``` bash
cd $path_to_wd/melasto689_on_full_ref
seqkit grep -w 0 -n -r -p "[[:upper:]]{4}-" Melastomataceae_689_original.FNA > Angio353_from_Melasto689.FNA;
seqkit grep -w 0 -n -r -v -p "[[:upper:]]{4}-" Melastomataceae_689_original.FNA > Melasto689_trimed.FNA;
Melasto689="Melasto689_trimed.FNA"
seqkit seq -n Melasto689_trimed.FNA > template_sequences_list.txt
```

### 2.3 Align every template sequence to the reference and extract the hits

Every template sequence (TS) from Melasto689 (`Melasto689_trimed.FNA`)
will be aligned to the set of references designed in Step 1. The aligner
used here is `BLAST` (`blastx-fast`). The regions of the TS that match
to a reference (the hits) are then extracted, the name of the reference
locus is appended to them, and they are place into the `extracted_hits`
folder. The non-matching regions are also separated into the `NMR_hits`
folder, and the TS with no match (i.e.??no portion of the TS matches a
reference) are placed into the `nomatch_TS` folder.

Prepare the sub-folders.

``` bash
cd $path_to_wd/melasto689_on_full_ref
mkdir extracted_hits
mkdir nomatch_TS
mkdir NMR_hits
```

For every TS in the list of template sequences:

-   create a sub-directory and set it as current directory
-   extract the TS from the whole Melasto689_trimmed file
-   align the TS to the reference with `blastx-fast` and store the
    result into a table
-   if the result table exists and is not empty, then:
    -   filter the results to keep only the matches with a bitscore
        greater than 50 and a query coverage greater than 50% (meaning
        at least 50% of the query covers the subject sequence), or
        matches with percentage identity greater than 60%, a query
        coverage greater than 10%, and a bitscore greater than 50
-   if the filtered table exists and is not empty, then:
    -   for every line in the table (i.e.??for every hit):
        -   retrieve the reference locus to which the hit aligns
            `refLocus`
        -   retrieve the frame in which the hit aligns to the reference
            `frame`
        -   print the matching subject (i.e.??the name of the reference
            TS to which the hit matched) into the file
            `matching_subject-$refLocus.txt`. This file will be used
            later during step 2.5
        -   if the `frame` value is greater than 0 (i.e.??the hit is not
            in reverse complement), then:
            -   set the hit into `.bed` format
            -   retrieve the sequence of the hit using
                `bedtools getfasta` and rename the sequence in order to
                have the old locus name separated with ???\_\_??? (double
                underscore) and the new locus name (`refLocus`) appended
                at the end of the name
            -   the hit is saved in `.FNA` (`<hitID>_hit.FNA`)
        -   if the `frame` value is lower than 0 (i.e.??the hit is in
            reverse complement), then:
            -   set the hit into `.bed` format
            -   retrieve the sequence of the hit using
                `bedtools getfasta`, reverse complement the sequence and
                rename the sequence in order to have the old locus name
                separated with ???\_\_??? (double underscore) and the new
                locus name (`refLocus`) appended at the end of the name
            -   the hit is saved in `.FNA` (`<hitID>_hit.FNA`)
        -   copy the hits (`<hitID>_hit.FNA` files) into
            `$path_to_wd/extracted_hits`
    -   retrieve the query coverage (i.e.??the percentage of query
        covered by the alignment)
    -   if the query coverage is lower than 90 (i.e.??at least 10% of the
        TS did not align to any reference), then:
        -   merge all the `.bed` files to get the region(s) of the TS
            that did match
        -   set a `.bed` file for the whole TS
        -   extract the coordinates of region(s) of the TS that did not
            match (non-matching regions, NMR) with `bedtools complement`
        -   retrieve the sequence of the NMR using `bedtools getfasta`,
            filter out the sequences less than 30 nucleotides long,
            append ???NMR??? to their locus name, and store them into
            `<TS_ID>NMR.FNA`
        -   copy the NMR (`<TS_ID>NMR.FNA` files) into
            `$path_to_wd/NMR_hits`
-   else (i.e.??the TS has no match):
    -   rename the TS to ???nomatch??? and store it into
        `$path_to_wd/nomatch_TS`
-   set the upper level directory as current directory

``` bash
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
      echo $line | awk -F ',' '{print $2}' >> matching_subject-$refLocus.txt
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
```

### 2.4 Merge the hits belonging to a same TS and a same locus

The extracted hits might or might not belong to a same TS and a same
locus. This step will merge the extracted hits that do. This is
important in order to keep a coherence in the cleaned TS. For example, a
TS that has a hit from the position 1 to 30 and another hit from the
position 60 to 90 has to be considered as a same TS and not two
separated TS. Note that for a same TS and same locus, the hits might
actually overlap. We will thus re-align them to the reference and draw a
consensus sequence of the hits alignment.

In the folder `extracted_hits`, create a new subfolder `cleaned_hits`
that will contain the merged hits.

``` bash
cd $path_to_wd/melasto689_on_full_ref/extracted_hits
mkdir cleaned_hits
```

Initialize a table containing the matches between the `locusCustomID`
and the corresponding Melasto689 locus name.

``` bash
echo "locusCustomID, correspondingLoci" > $path_to_wd/loci_matching_table.csv
```

Draw the the list of loci (containing `locusCustomID` for the TS that
matched to a reference in 2.3, and old locus name in case the TS had no
match).

``` bash
loci_list=$(ls -1 *.FNA | sed 's/.*-//g' | sed 's/_hit.FNA//g' | sort | uniq)
```

For every `locus` in the list of loci:

-   create a subfolder for the locus
-   copy the hits that belong to the locus in the subfolder
-   set the subfolder as current directory
-   retrieve the list of TS belonging to this locus `TS_list`
-   retrieve the list of old loci from the TS belonging to this locus
    `TS_loci`
-   write the locus and the corresponding old loci list into the table
    of matching loci
-   extract the sequences corresponding to the locus from the reference
    set
-   retrieve the names of these sequences and store them into a text
    file (later used by the
    [`align_translated_to_ref_and_draw_consensus.R`](scripts/align_translated_to_ref_and_draw_consensus.R)
    script)
-   for every TS in the `TS_list`:
    -   retrieve its old locus name `TS_locus`
    -   retrieve the number of hits belonging to the TS and locus
        `nb_hits`
    -   if the number of hits is greater than 1, then:
        -   concatenate these hits files together and with the reference
            sequences
        -   draw a consensus sequence for these hits using the
            [`align_translated_to_ref_and_draw_consensus.R`](scripts/align_translated_to_ref_and_draw_consensus.R)
            R script
        -   remove the gaps from the consensus sequence and rename it so
            that the old locus name (`TS_locus`) is bordered with ???\_\_???
            and so that the `-locus` name is at the end of the sequence
            name
    -   if the number of hits is equal to 1, then:
        -   rename the hit sequence so that the `-locus` name is at the
            end of the sequence name
-   copy the cleaned hits or merged hits into the `cleaned_hits`
-   set the upper level directory as current directory

``` bash
for locus in $loci_list
do
  #echo $locus
  mkdir $locus
  cp *-$locus"_hit.FNA" $locus
  cd $locus
  TS_list=$(ls -1 *hit.FNA | sed 's/:.*//g' | sort | uniq)
  TS_loci=$(ls -1 *hit.FNA | sed 's/:.*//g' | sed 's/.*__//g' | sort | uniq)
  echo $locus, $TS_loci >> $path_to_wd/loci_matching_table.csv
  seqkit grep -w0 -nrp $locus $path_to_wd/FULL_REFERENCE_nucl.FNA > $locus"_original.FNA"
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
```

At the end of this step we have the folder `cleaned_hits` full of
individual files with cleaned hits.

### 2.5 Add reference sequences to the final set

This step aims to add the relevant reference sequences (RS) to the final
probe set. We will only add the RS from loci that matched the TS in
2.3.  
Some loci have a very high number of transcriptome RS. Instead of
appending them all to the final set, we just select the transcriptome RS
that previously aligned to a TS.  
Concerning the Angio353 RS, no selection was carried out and all of them
were appended to the final set, if belonging to a locus present in the
final set.

-   create the `cleaned_sequences_set` folder
-   concatenate all the cleaned hit `.FNA` files into the file
    `cleaned_hits_all.FNA` in this folder

``` bash
cd $path_to_wd/melasto689_on_full_ref/extracted_hits/cleaned_hits
mkdir $path_to_wd/cleaned_sequences_set
cat *.FNA > $path_to_wd/cleaned_sequences_set/cleaned_hits_all.FNA
```

-   create the `matching_subjects` subfolder in the
    `cleaned_sequences_set` folder
-   fetch all the `matching_subject-$refLocus.txt` files and copy them
    into `$path_to_wd/cleaned_sequences_set/matching_subjects`
-   remove the duplicate lines in every text file

``` bash
cd $path_to_wd/melasto689_on_full_ref
mkdir cleaned_sequences_set/matching_subjects
find . -type f -name "matching_subject*.txt" -exec cp {} ./cleaned_sequences_set/matching_subjects \; # error messages are normal as the destination folder is embedded within the searching directory
cd cleaned_sequences_set/matching_subjects
for file in $(ls -1 *.txt)
do
  cat $file | sort | uniq > $file.temp
done
rm *.txt
rename -v '.txt.temp' '.txt' *
```

-   retrieve the list of loci from the concatenated clean hits file
-   initialize the `.FNA` file that will contain the final sequences and
    concatenate cleaned hits in it (`final_sequences_set.FNA`)
-   for every locus in the list of loci:
    -   retrieve the sequences that belong to this locus in the set of
        reference and get the number of those sequences that are from
        the trancriptome data (identified with ???scaffold??? in their name)
    -   if the number of these transcriptome sequences is greater than
        5, then:
        -   retrieve only the reference sequences that matched a TS in
            step 2.3 (whose ID are in the `matching_subject-locus.txt`)
            and append them to `final_sequences_set.FNA`
        -   retrieve all the reference sequences for this locus that
            come from Angio353 and if any was retrieved, append them to
            `final_sequences_set.FNA`
    -   else (i.e.??number of transcriptome sequences is lower than 5):
        -   retrieve all the reference sequences for this locus and
            append them to `final_sequences_set.FNA`

``` bash
cd $path_to_wd/melasto689_on_full_ref/cleaned_sequences_set
seqkit seq -n cleaned_hits_all.FNA | sed 's/.*-//g' | sort | uniq > loci_list.txt
cat cleaned_hits_all.FNA > final_sequences_set.FNA
for locus in $(cat loci_list.txt);
do
  #echo $locus
  nb_TTS=$(seqkit grep -w0 -nrp -$locus ../FULL_REFERENCE_nucl.FNA | seqkit seq -n | grep 'scaffold' | wc -l) # nb transcriptomes TS
  if [ $nb_TTS -gt 5 ]; then # if number of transcript (from 1KP) reference sequences is above 4 just select a few (ie the matching ones)
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
```

### 2.6 Append the Angio353 loci and remove duplicates

Get the list of names of Angio353 loci that were included in the
Melasto689 (i.e.??the list of angio353 loci) using file generated in 2.2.

``` bash
cd $path_to_wd/melasto689_on_full_ref/cleaned_sequences_set
seqkit seq -n ../Angio353_from_Melasto689.FNA | sed 's/.*-//g' | sort | uniq > Angio353_from_Melasto689_list.txt
```

For every locus in the list of angio353 loci:  
- retrieve the sequences from the reference set that belong to the locus
(i.e.??either having it as locus name or having it in their TS name) and
store them in `final_sequences_set.FNA`.

``` bash
for locus in $(cat Angio353_from_Melasto689_list.txt);
do
  #echo $locus
  seqkit grep -w0 -nrp -$locus ../FULL_REFERENCE_nucl.FNA >> final_sequences_set.FNA
  seqkit grep -w0 -nrp "__"$locus"__" ../FULL_REFERENCE_nucl.FNA >> final_sequences_set.FNA
done
```

Remove the duplicate sequences based on their name and sequence.

``` bash
seqkit rmdup -w0 -D name_duplicates.txt -n final_sequences_set.FNA > final_sequences_set_rmdups.FNA
seqkit rmdup -w0 -D sequence_duplicates.txt -s final_sequences_set_rmdups.FNA > final_sequences_set.FNA
```

### 2.7 Finalization

Sort the sequences by locus name in the final set.

``` bash
cd $path_to_wd/melasto689_on_full_ref/cleaned_sequences_set
cat final_sequences_set.FNA | seqkit rmdup -n | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' | seqkit sort -w0 -n | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > final_sequences_set.FNA.temp
rm final_sequences_set.FNA
```

Modify the names of some TS to limit the length of the names.

``` bash
seqkit replace -w0 -p 'hits_consensus' -r 'hc' final_sequences_set.FNA.temp > final_sequences_set.FNA.temp2
seqkit replace -w0 -p 'transcriptome_combo' -r 'tr_combo' final_sequences_set.FNA.temp2 > final_sequences_set.FNA
rm -r *temp*
```

Retrieve the TS names in the final set.

``` bash
cat final_sequences_set.FNA | seqkit seq -n > TS_full_list.txt
```

### 2.8 self-BLAST check

#### 2.8.1 Build BLAST database

``` bash
cd $path_to_wd/melasto689_on_full_ref
mkdir BLAST_DB
makeblastdb -in cleaned_sequences_set/final_sequences_set.FNA -out BLAST_DB/final_sequences_set_nucl -dbtype nucl -parse_seqids
makeblastdb -in cleaned_sequences_set/final_sequences_set.FAA -out BLAST_DB/final_sequences_set_prot -dbtype prot -parse_seqids
```

#### 2.8.2 Run the self-BLAST check

``` bash
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
```

The file `filtered_different_query_and_subject_loci_SB.csv` resulting
from the check self-BLAST contains a few sequences from different two
pairs of loci having high similarity.  
For the loci LOC12198 and sLOC23791 it concerns reference sequences
coming from Angio353 references (QSKP\_\_6363 and AEPI\_\_6320), but the
sequences are similar over less than 20% of their length.  
For the loci LOC20331 and LOC22857, it concerns reference sequences
coming from Angio353 and a TS from Melasto689. Their is relatively high
similarity (ca. 62%) over ca. 40% of the TS from Melasto689, but over
less than 10% for the reference sequences.  
Given the short region of similarity, here we decided to keep them
separated.

### 2.9 Remove transcriptome sequences that are supersets of mega353 sequences

The extension of the Angio353 set to the mega353 set appended some
sequences (actually part of sequences) from the SWGX and WWQZ
transcriptomes in the full reference set (Step 1: 2). These same
transcriptomes sequences were also appended in full to the full
reference set (Step 1: 1). Because of that, some sequences from the
mega353 are subsets of sequences coming from the transcriptomes. To
avoid redundancy, we remove the transcriptomes sequences that are
supersets of mega353 sequences.

Filter and select these sequences from the self-BLAST output table.

``` bash
cd $path_to_wd/melasto689_on_full_ref
cat selfblast_final_sequences_set_blastp_out_SB.csv | csvtk filter -H -f'14=100' | csvtk filter -H -f'3=100' | awk -F ',' 'OFS = "," {
if ($1 != $2)
print $0
}' | awk -F, '$2 ~ /scaffold/ {print}' | awk -F, '$1 ~ /^WWQZ|^SWGX/ {print}' > different_query_and_subject_SB.csv
cat different_query_and_subject_SB.csv
```

Retrieve the names of the transcriptomes sequences and remove them from
the final sequence set.

``` bash
cat different_query_and_subject_SB.csv | awk -F ',' '{print $2}' > TS_to_remove.txt
cd $path_to_wd/melasto689_on_full_ref/cleaned_sequences_set
seqkit grep -w0 -nrvf ../TS_to_remove.txt final_sequences_set.FNA > final_sequences_set.FNA.temp
rm final_sequences_set.FNA.
mv final_sequences_set.FNA.temp final_sequences_set.FNA
```

# Step 3. Final checks and fine-tuned cleaning

The final sequence set needs extra steps to make sure the sequences are
clean.

## 1. Treat the remaining stop codons individually

Some stop codons were found remaining in the Melasto689 TS. We decided
to assess them one by one. To do so, we aligned all the nucleotide TS
from a same locus. Unless specified, we replaced the stop codon with the
consensus sequence at the position of the stop codon.

Note that other stop codons are present in some reference sequences
coming from the mega353 dataset. They are treated later (see Step 3:
6.1).

### 1.1 Prepare the paths and data

``` bash
cd $path_to_wd
mkdir CLEAN_PROBE_SET
path_to_sequence_set=$path_to_wd/melasto689_on_full_ref
cp $path_to_sequence_set/cleaned_sequences_set/final_sequences_set.FNA ./CLEAN_PROBE_SET/
cd CLEAN_PROBE_SET
```

### 1.2 Treat the stop codons

All the following correction commands are stored in the script
[`custom_stop_codon_removal.sh`](scripts/custom_stop_codon_removal.sh).

#### Brachyotum\_\_5699\_\_hc-5699

TAG -\> NAN

``` bash
sed -i '/>Brachyotum__5699__hc-5699/,+1s/TTGCCCCTCAAATGGTAGTATCCC/TTGCCCCTCAAATGGNANTATCCC/g' final_sequences_set.FNA
```

#### Tibouchina\_\_AT3G59040\_\_hc-LOC00667

TAA -\> NAA

``` bash
sed -i '/>Tibouchina__AT3G59040__hc-LOC00667/,+1s/CGATCTGAAAKCTGGAGGAACATGAGATAA/CGATCTGAAAKCTGGAGGAACATGAGANAA/g' final_sequences_set.FNA
```

#### Affzelli\_\_5449\_\_hit-LOC00776

TGA -\> NNN

``` bash
sed -i '/>Affzelli__5449__hit-LOC00776/,+1s/ACTTGGGCTAAATTGTTGTGTTTTTGA/ACTTGGGCTAAATTGTTGTGTTTTNNN/g' final_sequences_set.FNA
```

#### Memecylon_tr_combo\_\_AT1G02020\_\_hc-LOC01351

Stop codon in a badling aligned portion at the beginning of a hit.
Remove the badly aligned portion (36 nucleotides).

``` bash
sed -i '/>Memecylon_tr_combo__AT1G02020__hc-LOC01351/,+1s/TGCAACAAAAGTAAATCCTGGTTGTATTTTATGTAG//g' final_sequences_set.FNA
```

#### Memecylon_tr_combo\_\_AT4G39520\_\_hc-LOC01520

TAG -\> NAG

``` bash
sed -i '/>Memecylon_tr_combo__AT4G39520__hc-LOC01520/,+1s/AGGCTAAACAAGTAGCCACCT/AGGCTAAACAAGNAGCCACCT/g' final_sequences_set.FNA
```

#### Tibouchina\_\_TC05G024510\_\_hc-LOC04226

2 stop codons in the middle of a badly aligning region with a lot of
ambiguous characters. Remove the whole region (27 nucleotides).

``` bash
sed -i '/>Tibouchina__TC05G024510__hc-LOC04226/,+1s/GGKTAGTTTATTCATTGAATTKKMSYT//g' final_sequences_set.FNA
```

#### Tibouchina\_\_AT4G00740\_\_hc-LOC04935

TAA -\> TAN

``` bash
sed -i '/>Tibouchina__AT4G00740__hc-LOC04935/,+1s/GCATGGTAACTTCTGATGAAACATTTTGTA/GCATGGTANCTTCTGATGAAACATTTTGTA/g' final_sequences_set.FNA
```

#### Tibouchina\_\_5639\_\_hc-LOC05244

TAG -\> TNG

``` bash
sed -i '/>Tibouchina__5639__hc-LOC05244/,+1s/TTKSTCTAGGGCAAT/TTKSTCTNGGGCAAT/g' final_sequences_set.FNA
```

#### Affzelli\_\_5426\_\_hit-LOC06188

TGA -\> TGG

``` bash
sed -i '/>Affzelli__5426__hit-LOC06188/,+1s/ATGTGAAAAAATGTTGGTAGT/ATGTGGAAAAATGTTGGTAGT/g' final_sequences_set.FNA
```

#### Tibouchina\_\_NM_179236.3\_\_hc-LOC06614

TAA -\> TNN

``` bash
sed -i '/>Tibouchina__NM_179236.3__hc-LOC06614/,+1s/TCAAGGTAACTTTTTGTGCTTTTTGCA/TCAAGGTNNCTTTTTGTGCTTTTTGCA/g' final_sequences_set.FNA
```

The other stop codon is in a badly aligning region with ambiguous
characters. Remove the region (21 nucleotides).

``` bash
sed -i '/>Tibouchina__NM_179236.3__hc-LOC06614/,+1s/YTSCTGCTGCTATGGTGACAG//g' final_sequences_set.FNA
```

#### Tibouchina\_\_5339\_\_hc-LOC07010

Remove the stop codon and the region with a lot of ambiguous characters
just before (48 nucleotides).

``` bash
sed -i '/>Tibouchina__5339__hc-LOC07010/,+1s/CAKACCAGGKATKKWTYTCMMCTTRAAAAAATCTAAACGAGGCAGTTG//g' final_sequences_set.FNA
```

#### Memecylon_tr_combo\_\_TC04G019470\_\_hc-LOC07251

TAG -\> TNN

``` bash
sed -i '/>Memecylon_tr_combo__TC04G019470__hc-LOC07251/,+1s/TCTATGTAGAGATGGGTTGCAGTTGTTTTTGTATTC/TCTATGTNNAGATGGGTTGCAGTTGTTTTTGTATTC/g' final_sequences_set.FNA
```

#### Tibouchina\_\_6563\_\_hc-LOC07413

TAA -\> NAN

``` bash
sed -i '/>Tibouchina__6563__hc-LOC07413/,+1s/TCCAAGCTTCCAAAGTAAAGGACTTCTGTGGTGACAATATTG/TCCAAGCTTCCAAAGNANAGGACTTCTGTGGTGACAATATTG/g' final_sequences_set.FNA
```

#### Torricelli_tr_combo\_\_TC01G021540\_\_hc-LOC07589

TGA -\> NGN

``` bash
sed -i '/>Torricelli_tr_combo__TC01G021540__hc-LOC07589/,+1s/TCATGACTTTTCTGTGCATCCATCAACTTCCTA/TCANGNCTTTTCTGTGCATCCATCAACTTCCTA/g' final_sequences_set.FNA
```

#### Tibouchina\_\_NM_001198514.1\_\_hc-LOC08266

Remove the stop codon and the region with a lot of ambiguous characters
just before.

``` bash
sed -i '/>Tibouchina__NM_001198514.1__hc-LOC08266/,+1s/TTYWYSSYSTYYYSTKYKTTCTGA//g' final_sequences_set.FNA
```

#### Memecylon_tr_combo\_\_TC09G010970\_\_hc-LOC08580

The whole tail of the sequence, containing 2 stop codons badly aligns to
the rest of the TSs. Remove the whole tail (51 nucleotides).

``` bash
sed -i '/>Memecylon_tr_combo__TC09G010970__hc-LOC08580/,+1s/GTAAGTGTTCATTATTTCACACGCCTTCTGAGTTGATTTTGTTGATTCTGT//g' final_sequences_set.FNA
```

#### Tibouchina\_\_TC05G027690\_\_hc-LOC09163

Stop codon in the tail of the sequence that badly aligns to the rest of
the TS. Remove the tail (42 nucleotides).

``` bash
sed -i '/>Tibouchina__TC05G027690__hc-LOC09163/,+1s/ATGATCGGTAGAGCACACGTCTGAACTCCAGTCACGAGATTC//g' final_sequences_set.FNA
```

#### Tibouchina\_\_7572\_\_hc-LOC09834

The stop codon is at the end of a hit in a badly aligning region. Remove
the region (24 nucleotides).

``` bash
sed -i '/>Tibouchina__7572__hc-LOC09834/,+1s/CAAGGAGTTATAACGTGACTCAAG//g' final_sequences_set.FNA
```

#### Tibouchina\_\_5942\_\_hc-LOC14673

The other stop codon is in a badly aligning region with ambiguous
characters. Remove the region (54 nucleotides).

``` bash
sed -i '/>Tibouchina__5942__hc-LOC14673/,+1s/RGCMRRRRKCMSYRSSSSYRGMARRARRCMSKWMGMCCAGAAGAAGACCGGTAG//g' final_sequences_set.FNA
```

#### Tibouchina\_\_AT5G13630\_\_hc-LOC15950

Stop codon in a region of the sequence that badly aligns to the rest of
the TS. Remove the region (33 nucleotides).

``` bash
sed -i '/>Tibouchina__AT5G13630__hc-LOC15950/,+1s/GTATAGCCGGCTGMSMKRGYSRTSAARRTKGTK//g' final_sequences_set.FNA
```

#### Brachyotum\_\_AT3G02760\_\_hc-LOC24270

Stop codon in a region of the sequence that badly aligns to the rest of
the TS. Remove the region (33 nucleotides).

``` bash
sed -i '/>Brachyotum__AT3G02760__hc-LOC24270/,+1s/AATTCCTAGAATTCATKSKKKCKYTAWKRWRRCSWKMWA//g' final_sequences_set.FNA
```

## 2. Recycle the non matching TS

The TS that do not align to any sequence in the reference set can still
be included in the final set. We need to make sure they are in the
correct frame and that they do not include any stop codon. To do so, we
will use a [script developed by Chris
Jackson](https://github.com/mossmatters/HybPiper/issues/82#issuecomment-1158575310)
that identifies the first forward frame that has no stop codon and
recovers the in-frame sequence while removing trailing 3??? nucleotides.

``` bash
path_to_sequence_set=$path_to_wd/melasto689_on_full_ref
cd $path_to_sequence_set/nomatch_TS
cat *nomatch.FNA > nomatches_all.FNA
python /blue/soltis/dagallierl/PROGRAMS/get_inframe_targetfile.py nomatches_all.FNA
```

The retrieved in-frame sequences are renamed so that the old locus name
is included in the TS name and so that their new locus name is in the
form `-outlier###` (with `#` being a digit).

``` bash
cat nomatches_all_inframe.FNA | sed 's/-/_/g' | sed 's/\./_/g' | seqkit replace -w0 -p $ -r "-outlier{nr}" --nr-width 3 > nomatches_recovered.FNA
```

We still want to check whether these recovered TS could actually align
to a sequence in the final set. To do so, we extract their old locus
name and for each locus we retrieve all the sequences from the final set
that belong to this locus into a separate file.

``` bash
seqkit seq -n nomatches_recovered.FNA | sed 's/-.*//g' | sed 's/.*_//g'
> loci_to_check.txt
for locus in $(cat loci_to_check.txt)
do
  seqkit grep -w0 -nrp $locus final_sequences_set.FNA > $locus"_test.FNA"
  seqkit grep -w0 -nrp $locus nomatches_recovered.FNA >> $locus"_test.FNA"
done
```

We then check manually each file by aligning the sequences for every
locus and to look visually if the ???outlier??? sequence align to the rest
of the sequences. To do so we used [Muscle
v3.8.425](https://drive5.com/muscle/) as wrapped in the program
[AliView](https://ormbunkar.se/aliview/) to visually assess the
alignement.

It turns out that the sequence from the locus ???outliers015??? aligns to
the locus ???5463???. We thus change the locus name of this sequence to
???5463???.

``` bash
seqkit replace -w0 -p '-outlier015' -r '-5463' nomatches_recovered.FNA | sed 's/_5463-5463/-5463/g' > nomatches_recovered_renammed.FNA
rm -r *test.FNA
```

## 3. Rename the loci names

To simplify the process of building the set of reference sequences (Step
1), we used custom loci names in the form `LOC######` (where `#`
represents a digit).

In order to ease the readibility of the final probe set, we revert the
loci names to their old name (the one used in Melasto689 and Angio353).

### 3.1 Names matching table

Prepare a matching table of the loci names between the original names
and the custom loci names.

``` bash
cd $path_to_wd/CLEAN_PROBE_SET
path_to_sequence_set=$path_to_wd/melasto689_on_full_ref
cp $path_to_sequence_set/cleaned_sequences_set/TS_full_list.txt ./
```

Retrieve the list of custom loci names.

``` bash
cat TS_full_list.txt | sed 's/.*-//g' | sort | uniq > custom_loci_list.txt
```

Initialize the matching table and for every locus in the list of custom
loci names:  
- retrieve the corresponding locus (or loci) name(s)  
- populate the matching table with the custom loci name and its
corresponding locus (loci) name(s)

``` bash
:> locus_matching_table.csv
for customLocus in $(cat custom_loci_list.txt)
do
  correspondingLocus=$(grep -e -$customLocus TS_full_list.txt | grep -Po '(?<=__).*(?=__)' | sort | uniq)
  echo $customLocus, $correspondingLocus >> locus_matching_table.csv
done
```

Several custom loci have more than one corresponding locus. This happens
because in the Melasto689 probe set, TS from Angio353 loci actually have
high similarity to sequences from MarkerMiner loci. In these cases, we
manually check the sequences: for each locus we align the TS of the
locus and check whether they all align correctly together. If they do we
keep them under the same name (we arbitrarily choose to keep the
Angio353 name). If they do not align correctly together, we separate
them into different names. Here is the detail for these loci:

**customLocusID: originalLocusID = finalLocusID**  
LOC00757: KT377110.1 = 5870  
LOC07406: AT5G60590 = 6570  
LOC13529: AT1G31780 = 6098  
LOC14778: TC04G006100 = 4954  
LOC15996: AT5G06830 = 6175  
LOC16890: AT2G44660 = 6284  
LOC18108: AT2G43030 = 5857  
LOC19586: TC02G003400 = 6875  
LOC22857: AT5G67530 = 5264  
In the list above, customLocusID is the custom locus name in the
reference set, originalLocusID is the MarkerMiner locus name as in
Melasto689, and finalLocusID is the final locus name.

For the locus LOC00263, the template sequences should be re-separated
into different loci:  
**old TS name = new TS name**  
Tetrazygia\_\_NM_113708.2\_\_hc-LOC00263 = Tetrazygia_hc-NM_113708.2  
Tetrazygia\_\_NM_124699.3\_\_hit-LOC00263 = Tetrazygia_hit-NM_124699.3  
Tibouchina\_\_NM_113708.2\_\_hit-LOC00263 = Tibouchina_hit-NM_113708.2  
Tibouchina\_\_NM_124699.3\_\_hit-LOC00263 = Tibouchina_hit-NM_124699.3  
scaffold_WWQZ_2004803_Medinilla_magnifica-LOC00263 =
scaffold_WWQZ_2004803_Medinilla_magnifica-NM_113708.2  
scaffold_SWGX_2016500_Tetrazygia_bicolor-LOC00263 =
scaffold_SWGX_2016500_Tetrazygia_bicolor-NM_113708.2

We end up with the manually modified .csv file
[`locus_matching_table_CUSTOM.csv`](CLEAN_PROBE_SET/locus_matching_table_CUSTOM.csv)
containing the string to replace on the first column and the replacement
string on the second column.

### 3.2 Renaming

Build a replacement table that can be read by `seqkit replace`.

``` bash
cat locus_matching_table_CUSTOM.csv | sed 's/ //g' | awk -F ',' '{print $1"\t"$2}' > replacement.txt
```

Rename the sequences using `seqkit replace` and draw the list of all the
TS. This is done in 2 steps because for some TS only the locus name is
changed but for other TS the full TS name is changed.

``` bash
seqkit replace -w0 -p '-(.+)$' -r '-{kv}' -k replacement.txt --keep-key final_sequences_set.FNA > final_sequences_set_renamed.FNA

seqkit seq -n final_sequences_set_renamed.FNA > TS_full_list_renamed.txt

cat locus_matching_table_CUSTOM.csv | grep -e 'NM_113708.2\|NM_124699.3' | sed 's/ //g' | awk -F ',' '{print $1"\t"$2}' > extra_replacement.txt
seqkit replace -w0 -p '(.+)' -r '{kv}' -k extra_replacement.txt --keep-key final_sequences_set_renamed.FNA > final_sequences_set_renamed2.FNA
seqkit seq -n final_sequences_set_renamed2.FNA > TS_full_list_renamed2.txt
```

Most of the TS still contain their old locus name (which is actually
their final locus name as the the custom locus name were removed)
bordered with the ???\_\_??? (double underscore). We thus remove all the
???\_\_locus\_\_??? from the TS names.

``` bash
loci_list=$(seqkit seq -n final_sequences_set_renamed2.FNA | sed 's/.*-//g' | sort | uniq)
:> final_sequences_set_renamed3.FNA
for locus in $loci_list
do
  seqkit grep -w0 -nrp -$locus final_sequences_set_renamed2.FNA | sed '/^>/s/__'$locus'__/_/g' >> final_sequences_set_renamed3.FNA
done
sed -i 's/_-/-/g' final_sequences_set_renamed3.FNA
sed -i 's/__/_/g' final_sequences_set_renamed3.FNA

seqkit seq -n final_sequences_set_renamed3.FNA > TS_full_list_renamed3.txt
```

## 4. Append the recycled non matching TS

``` bash
cat final_sequences_set_renamed3.FNA nomatches_recovered_renammed.FNA > final_sequences_set_renamed4.FNA
```

## 5. Finalization

Sort the sequences according to locus name.

``` bash
cat final_sequences_set_renamed4.FNA | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' | seqkit sort -w0 -n | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > final_sequences_set_TO_CHECK.FNA
```

Translate the nucleotide sequence set into amino-acids.

``` bash
seqkit translate -w0 final_sequences_set_TO_CHECK.FNA > final_sequences_set_TO_CHECK_prot.FAA
```

The sequence set is ready to be checked by HybPiper.

## 6. Check the sequence set with HybPiper

HybPiper v2 has a function that allow to check the target file for low
complexity sequences and for the presence of stop codons and sequences
not multiple of 3 (see [HybPiper
wiki](https://github.com/mossmatters/HybPiper/wiki/Full-pipeline-parameters#70-hybpiper-check_targetfile)
for more details).

``` bash
cd $path_to_wd/CLEAN_PROBE_SET
hybpiper check_targetfile --targetfile_dna final_sequences_set_TO_CHECK.FNA > check_targetfile.log
hybpiper check_targetfile --targetfile_aa final_sequences_set_TO_CHECK_prot.FAA > check_targetfile_aa.log
```

Some sequences still have stop codons and are not multiples of 3. These
sequences are reference sequences that come from the mega353 dataset.

### 6.1 Correct the sequences that need to

Manually get the list of TS to correct from the log file, and store them
in the file `TS_to_correct.txt`.  
Separate the TS that need to be corrected (`TS_to_correct.FNA`) from
those that are ok `TS_checked_OK.FNA`.

``` bash
cd $path_to_wd/CLEAN_PROBE_SET
cat TS_to_correct.txt | sed 's/$/\n/' | sed 's/, /\n/g' > list_TS_to_correct.txt
seqkit grep -w0 -nrf list_TS_to_correct.txt final_sequences_set_TO_CHECK.FNA > TS_to_correct.FNA
seqkit grep -w0 -nrvf list_TS_to_correct.txt final_sequences_set_TO_CHECK.FNA > TS_checked_OK.FNA
```

Manually correct `TS_to_correct.FNA` into
`TS_to_correct_CORRECTED.FNA`: - remove the stop codon and the trailing
3??? nucleotides - for sequences not multiple of 3, remove the 3??? hanging
nucleotide(s)

All the correction commands are stored in the script
[`custom_TS_correction.sh`](scripts/custom_TS_correction.sh).

Then merge back `TS_checked_OK.FNA` with `TS_to_correct_CORRECTED.FNA`.

``` bash
cat TS_checked_OK.FNA TS_to_correct_CORRECTED.FNA > final_sequences_set_corrected.FNA
```

### 6.2 Remove low complexity sequences

Manually get the list of TS with low complexity from
`check_targetfile.log` and `check_targetfile_aa.log` store it in the
file `list_of_TS_with_low_complexity.txt`. For This step removes the TS
flagged to present low complexity for a given locus, but only if the
locus is also represented by other TS (not flagged with low complexity).

-   remove duplicates from `list_of_TS_with_low_complexity.txt`
-   retrieve the list of loci for the list of TS
-   initialize the file `TS_with_low_complexity_to_remove.txt`
-   for every locus in `list_of_TS_with_low_complexity.txt`:
    -   retrieve the full list of TS for this locus
    -   retrieve the list of TS to remove for this locus (i.e.??TS with
        low complexity)
    -   retrieve the list of remaining TS (i.e.??the full list minus the
        TS to remove)
    -   if there are any remaining TS in the list, then:
        -   append the list of TS to remove to the file
            `TS_with_low_complexity_to_remove.txt`
-   remove the low complexity TS from the final sequence set using the
    file `TS_with_low_complexity_to_remove.txt`

``` bash
dos2unix list_of_TS_with_low_complexity.txt
cat list_of_TS_with_low_complexity.txt | sort | uniq > list_of_TS_with_low_complexity.txt.temp
rm list_of_TS_with_low_complexity.txt
mv list_of_TS_with_low_complexity.txt.temp list_of_TS_with_low_complexity.txt

loci_list=$(cat list_of_TS_with_low_complexity.txt | sed 's/.*-//g' | sort | uniq)
:> TS_with_low_complexity_to_remove.txt
for locus in $loci_list;
do
  echo $locus
  seqkit grep -w0 -nrp -$locus final_sequences_set_corrected.FNA | seqkit seq -n > TS_full_list
  grep -e -$locus list_of_TS_with_low_complexity.txt > TS_to_remove
  cat TS_full_list TS_to_remove | sort | uniq -u > TS_remaining
    if [ -s TS_remaining ]; then
      echo Remaining TS exist for $locus, removing the low complexity TS...
      cat TS_to_remove >> TS_with_low_complexity_to_remove.txt
    else
      echo No remaining TS for $locus, keeping it as is...
    fi
  rm TS_full_list TS_to_remove TS_remaining
done

seqkit grep -w0 -nvf TS_with_low_complexity_to_remove.txt final_sequences_set_corrected.FNA > final_sequences_set_corrected_low_complexity_removed.FNA
```

All the TS flagged for low complexity were removed, except for locus
6128, for which all the TSs were flagged as presenting low complexity.

### 6.3 Finalize the file and translate

Sort the sequences according to locus name.

``` bash
cat final_sequences_set_corrected_low_complexity_removed.FNA | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' | seqkit sort -w0 -n | sed -E '/^>/s/(\w*\.*\w*)-(\w*\.*\w*)/\2-\1/' > final_sequences_set_corrected_low_complexity_removed_sorted.FNA
```

Translate to amino-acids.

``` bash
cat final_sequences_set_corrected_low_complexity_removed_sorted.FNA > PROBE_SET_CLEAN.FNA
seqkit translate -w0 PROBE_SET_CLEAN.FNA > PROBE_SET_CLEAN_prot.FAA
```

### 6.4 Extra check with HybPiper

This is an extra final check that the final sequence set is ok for
HybPiper.

``` bash
cd $path_to_wd/CLEAN_PROBE_SET
hybpiper check_targetfile --targetfile_dna PROBE_SET_CLEAN.FNA > check_targetfile_final.log
hybpiper check_targetfile --targetfile_aa PROBE_SET_CLEAN_prot.FAA > check_targetfile_aa_final.log
```

The sequence set is now cleaned and can be used with HybPiper or any
other program. **The final files are
[`PROBE_SET_CLEAN.FNA`](CLEAN_PROBE_SET/PROBE_SET_CLEAN.FNA)
(nucleotides) and
[`PROBE_SET_CLEAN_prot.FAA`](CLEAN_PROBE_SET/PROBE_SET_CLEAN_prot.FAA)
(amino-acids).**

# References

-   Carpenter EJ, Matasci N, Ayyampalayam S, Wu S, Sun J, Yu J,
    Jimenez??Vieira FR, Bowler C, Dorrell RG, Gitzendanner MA, Li L, Du
    W, K.??Ullrich K, Wickett NJ, Barkmann TJ, Barker MS, Leebens-Mack
    JH, Wong GK-S (2019) Access to RNA-sequencing data from 1,173 plant
    species: The 1000 Plant transcriptomes initiative (1KP). GigaScience
    8: giz126. <https://doi.org/10.1093/gigascience/giz126>
-   Jantzen JR, Amarasinghe P, Folk RA, Reginato M, Michelangeli FA,
    Soltis DE, Cellinese N, Soltis PS (2020) A two-tier bioinformatic
    pipeline to develop probes for target capture of nuclear loci with
    applications in Melastomataceae. Applications in Plant Sciences 8:
    e11345. <https://doi.org/10.1002/aps3.11345>
-   Johnson MG, Pokorny L, Dodsworth S, Botigu?? LR, Cowan RS, Devault A,
    Eiserhardt WL, Epitawalage N, Forest F, Kim JT, Leebens-Mack JH,
    Leitch IJ, Maurin O, Soltis DE, Soltis PS, Wong GK, Baker WJ,
    Wickett NJ (2019) A Universal Probe Set for Targeted Sequencing of
    353 Nuclear Genes from Any Flowering Plant Designed Using k-Medoids
    Clustering. Systematic Biology 68: 594???606.
    <https://doi.org/10.1093/sysbio/syy086>
-   McLay TGB, Birch JL, Gunn BF, Ning W, Tate JA, Nauheimer L, Joyce
    EM, Simpson L, Schmidt-Lebuhn AN, Baker WJ, Forest F, Jackson
    CJ (2021) New targets acquired: Improving locus recovery from the
    Angiosperms353 probe set. Applications in Plant Sciences 9.
    <https://doi.org/10.1002/aps3.11420>
-   One Thousand Plant Transcriptomes Initiative (2019) One thousand
    plant transcriptomes and the phylogenomics of green plants. Nature
    574: 679???685. <https://doi.org/10.1038/s41586-019-1693-2>
-   Reginato M, Michelangeli FA (2016) Primers for low-copy nuclear
    genes in the Melastomataceae. Applications in Plant Sciences
    4: 1500092. <https://doi.org/10.3732/apps.1500092>
