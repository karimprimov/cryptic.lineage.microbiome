
#260522:
######################################################################################################################################################
######################################################################################################################################################

#1 Acropora tenuis
export BioProject=PRJNA849642
$HOME/edirect/esearch -db sra -query $BioProject | efetch -format runinfo | cut -d "," -f 1 | grep SRR > PRJNA849642.SRR

export accs=$(cat PRJNA849642.SRR)
>fetch
for i in $accs; do
echo "prefetch $i --max-size u --output-directory ./sra" >> fetch; done
ls6_launcher_creator.py -j fetch -n fetch -a IBN21018 -e kdprimov@utexas.edu -t 48:00:00 -w 48 -q normal
sbatch fetch.slurm


#260525:
#got an error message about this prefetch job, fixing it now:
awk '{print "prefetch " $1 " --max-size u --output-directory ./tenuis.new.sra"}' PRJNA849642.SRR > tenuis.new.fetch
python3 $HOME/bin/ls6_launcher_creator.py -j tenuis.new.fetch -n tenuis.new.fetch -a IBN21018 -e kdprimov@utexas.edu -t 48:00:00 -w 4 -N 5 -q normal
sbatch tenuis.new.fetch.slurm

sed -i 's/\r//g' PRJNA849642.SRR
export accs=$(cat PRJNA849642.SRR)
> atenuis.dump
for i in $(cat PRJNA849642.SRR); do
  echo "fasterq-dump --split-files --threads 8 --temp /scratch/08554/kprimov/tmp_$i ./tenuis.new.sra/$i/$i.sra -O ./ && pigz -p 8 ${i}_*.fastq" >> atenuis.dump
done
ls6_launcher_creator.py -j atenuis.dump -n atenuis.dump -a IBN21018 -e kdprimov@utexas.edu -t 24:00:00 -w 8 -N 3 -q normal
sbatch atenuis.dump.slurm


#260529: trimgalore
export allo=IBN21018 #change to your allocation
export email=kdprimov@utexas.edu #change to your email
>260529.acro.tenuis.trimgalore
for file in *_1.fastq.gz; do
  base=${file/_1.fastq.gz/}
  echo "trim_galore --paired --nextera --nextseq 20 -o 260529.tenuis.trimmed.output ${base}_1.fastq.gz ${base}_2.fastq.gz" >> 260529.acro.tenuis.trimgalore; done
python3 $HOME/bin/ls6_launcher_creator.py -j 260529.acro.tenuis.trimgalore -n 260529.acro.tenuis.trimgalore -a $allo -e $email -t 02:00:00 -N 5 -w 8 -q development
sbatch 260529.acro.tenuis.trimgalore.slurm

#260530:
#fastqc code:
#pwd:/scratch/08554/kprimov/260522.atenuis.data/260529.tenuis.trimmed.output
export allo=IBN21018 #change to your allocation
export email=kdprimov@utexas.edu #change to your email
>260530.fastqc
for i in *.fq.gz; do
  echo "fastqc $i" >> 260530.fastqc; done
python3 $HOME/bin/ls6_launcher_creator.py -j 260530.fastqc -n 260530.fastqc -a $allo -e $email -t 02:00:00 -N 1 -w 16 -q development
sbatch 260530.fastqc.slurm

#copying the fastqc files over to local CPU:
TACC_SYS=ls6
# TACC_SYS=ls5
TACC_USER=kprimov
TACC_URL=tacc.utexas.edu
TACC_SCRATCH=/scratch/08554/kprimov
TACC_LOGON=$TACC_USER@$TACC_SYS.$TACC_URL
SCP_SCRATCH=$TACC_LOGON:$TACC_SCRATCH
# Set your TACC and Local file directories
TACC_DIR=260522.atenuis.data/260529.tenuis.trimmed.output
LOCAL_DIR=/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/260530.fastqc/atenuis
# The file/files you want to copy
FILE=*.html
#Copy from TACC to local
scp $SCP_SCRATCH/$TACC_DIR/$FILE $LOCAL_DIR/


#260601:
#alignment:
#pwd:/scratch/08554/kprimov/260522.atenuis.data/260529.tenuis.trimmed.output
export GENOME_FASTA=/work/08554/kprimov/ls6/concatenated.acropora.tenuis.ref.genome/acropora.tenuis.concat.ref
export allo=IBN21018 #change to your allocation
export email=kdprimov@utexas.edu #change to your email
>260601.tenuis.align
for i in *_1_val_1.fq.gz; do
    base=$(basename $i "_1_val_1.fq.gz")
    REF=$GENOME_FASTA
    # Align reads and extract unaligned reads using Bowtie2
    echo "bowtie2 -x $REF -1 ${base}_1_val_1.fq.gz -2 ${base}_2_val_2.fq.gz -S ${base}.coral.aligned.sam --un-conc ${base}_unaligned.fastq -p 16" >> 260601.tenuis.align
done
python3 $HOME/bin/ls6_launcher_creator.py -j 260601.tenuis.align -n 260601.tenuis.align -a $allo -e $email -t 36:00:00 -N 6 -w 4 -q normal
sbatch 260601.tenuis.align.slurm



######################################################################################################################################################
######################################################################################################################################################



#2 Acropora hyacinthus
export BioProject=PRJNA657822
$HOME/edirect/esearch -db sra -query $BioProject | efetch -format runinfo | cut -d "," -f 1 | grep SRR > PRJNA657822.SRR

export accs=$(cat PRJNA657822.SRR)
>fetch
for i in $accs; do
echo "prefetch $i --max-size u --output-directory ./sra" >> fetch; done
ls6_launcher_creator.py -j fetch -n fetch -a IBN21018 -e kdprimov@utexas.edu -t 24:00:00 -w 48 -q normal
sbatch fetch.slurm

#260525:
#got an error message about this prefetch job, fixing it now:
awk '{print "prefetch " $1 " --max-size u --output-directory ./sra"}' PRJNA657822.SRR > fetch
python3 $HOME/bin/ls6_launcher_creator.py -j fetch -n fetch -a IBN21018 -e kdprimov@utexas.edu -t 36:00:00 -w 4 -N 3 -q normal
sbatch fetch.slurm


sed -i 's/\r//g' PRJNA657822.SRR
export accs=$(cat PRJNA657822.SRR)
> ahyacinthus.dump
for i in $(cat PRJNA657822.SRR); do
  echo "fasterq-dump --split-files --threads 8 --temp /scratch/08554/kprimov/tmp_$i ./sra/$i/$i.sra -O ./ && pigz -p 8 ${i}_*.fastq" >> ahyacinthus.dump
done
ls6_launcher_creator.py -j ahyacinthus.dump -n ahyacinthus.dump -a IBN21018 -e kdprimov@utexas.edu -t 12:00:00 -w 8 -N 2 -q normal
sbatch ahyacinthus.dump.slurm


#260529 trimgalore:
export allo=IBN21018 #change to your allocation
export email=kdprimov@utexas.edu #change to your email
>260529.acro.hyacinthus.trimgalore
for file in *_1.fastq.gz; do
  base=${file/_1.fastq.gz/}
  echo "trim_galore --paired --nextera --nextseq 20 -o 260529.hyacinthus.trimmed.output ${base}_1.fastq.gz ${base}_2.fastq.gz" >> 260529.acro.hyacinthus.trimgalore; done
python3 $HOME/bin/ls6_launcher_creator.py -j 260529.acro.hyacinthus.trimgalore -n 260529.acro.hyacinthus.trimgalore -a $allo -e $email -t 02:00:00 -N 5 -w 16 -q development
sbatch 260529.acro.hyacinthus.trimgalore.slurm

#260530:
#fastqc code:
#pwd:/scratch/08554/kprimov/260522.atenuis.data/260529.tenuis.trimmed.output
export allo=IBN21018 #change to your allocation
export email=kdprimov@utexas.edu #change to your email
>260530.fastqc
for i in *.fq.gz; do
  echo "fastqc $i" >> 260530.fastqc; done
python3 $HOME/bin/ls6_launcher_creator.py -j 260530.fastqc -n 260530.fastqc -a $allo -e $email -t 02:00:00 -N 1 -w 16 -q development
sbatch 260530.fastqc.slurm

#copying the fastqc files over to local CPU:
TACC_SYS=ls6
# TACC_SYS=ls5
TACC_USER=kprimov
TACC_URL=tacc.utexas.edu
TACC_SCRATCH=/scratch/08554/kprimov
TACC_LOGON=$TACC_USER@$TACC_SYS.$TACC_URL
SCP_SCRATCH=$TACC_LOGON:$TACC_SCRATCH
# Set your TACC and Local file directories
TACC_DIR=260522.ahyacinthus.data/260529.hyacinthus.trimmed.output
LOCAL_DIR=/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/260530.fastqc/ahyacinthus
# The file/files you want to copy
FILE=*.html
#Copy from TACC to local
scp $SCP_SCRATCH/$TACC_DIR/$FILE $LOCAL_DIR/

#260601:
#alignment:
#pwd:/scratch/08554/kprimov/260522.ahyacinthus.data/260529.hyacinthus.trimmed.output
export GENOME_FASTA=/work/08554/kprimov/ls6/260219.amil.wgs.symb.concat.genome/a.hyacinthus.a.millepora.ref
export allo=IBN21018 #change to your allocation
export email=kdprimov@utexas.edu #change to your email
>260601.hyacinthus.align
for i in *_1_val_1.fq.gz; do
    base=$(basename $i "_1_val_1.fq.gz")
    REF=$GENOME_FASTA
    # Align reads and extract unaligned reads using Bowtie2
    echo "bowtie2 -x $REF -1 ${base}_1_val_1.fq.gz -2 ${base}_2_val_2.fq.gz -S ${base}.coral.aligned.sam --un-conc ${base}_unaligned.fastq -p 16" >> 260601.hyacinthus.align
done
python3 $HOME/bin/ls6_launcher_creator.py -j 260601.hyacinthus.align -n 260601.hyacinthus.align -a $allo -e $email -t 24:00:00 -N 4 -w 4 -q normal
sbatch 260601.hyacinthus.align.slurm

######################################################################################################################################################
######################################################################################################################################################
#tara oceans data:
#pwd: /scratch/08554/kprimov/2604.tara.oceans.data
#pwd for truncated sam files: /scratch/08554/kprimov/2604.tara.oceans.data/tara.coral.samples/260408.tara.trim.output/trimmed.labeled.files/aligned.sam.files
#need to redownload these for proper pairwise fst calculation

wget "https://zenodo.org/record/6299409/files/TARA-PACIFIC_samples-provenance_20200731d.txt?download=1"
mv TARA-PACIFIC_samples-provenance_20200731d.txt* TARA-PACIFIC_samples-provenance_20200731d.txt
# 1. Filter rows using column 57 specifically
awk -F'\t' 'NR==1 || $57 ~ /Porites/ || $57 ~ /Pocillopora/' TARA-PACIFIC_samples-provenance_20200731d.txt > filtered_coral_rows.tsv
# 2. Extract ENA Sample Accessions (Column 2) from your newly filtered file
awk -F'\t' 'NR>1 {print $2}' filtered_coral_rows.tsv > coral_ena_ids.txt
wc -l filtered_coral_rows.tsv
#unique rows 
awk -F'\t' '$59 ~ /SEQ/ {print $2}' filtered_coral_rows.tsv | sort -u > coral_ena_ids.txt


#extracting the sra run data to get the ERR numbers for proper downloading of the Porites and Pocillopora samples:
wget "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=ERP137075" -O ERP137075_runinfo.csv





#trying to get all the metadata to cleanly extract the 200+ samples of Porites and Pocillopora
grep -E "Porites|Pocillopora" TARA-PACIFIC_samples-provenance_20200731d.txt > coral.samples.txt
cut -f3 coral.samples.txt | sort | uniq > coral_biosamples.txt
grep -F -f coral_biosamples.txt ERP137075_runinfo.csv > coral_runs.csv
cut -d',' -f1 coral_runs.csv > ERR_ids.txt
# 1. Convert your newline-separated ERR IDs into a comma-separated string for the API
ERR_LIST=$(paste -sd, ERR_ids.txt)
curl -L "https://ebi.ac.uk{ERR_LIST}&fields=fastq_ftp&format=tsv" | tail -n +2 | tr ';' '\n' | sed 's/^/ftp:\/\//' > coral_download_urls.txt
vim download.slurm

#slurm file:

#!/bin/bash
#SBATCH -J tara_download          # Job name
#SBATCH -n 24                     # Core/task count matching your standard template
#SBATCH -N 1                      # Set to 1 node (1 node has full access to the proxy pipe)
#SBATCH -p normal                 # Back to your standard queue
#SBATCH -o tara_download.o%j       # Standard output log
#SBATCH -e tara_download.e%j       # Error log
#SBATCH -t 48:00:00               # Keeping your 48-hour window
#SBATCH -A IBN21018               # Your active allocation ID verified by the system

# Change to the directory where you submitted the job (Ensure this is in $WORK or $SCRATCH)
cd $SLURM_SUBMIT_DIR

# 1. Convert your newline-separated ERR IDs into a comma-separated string for the ENA API
ERR_LIST=$(paste -sd, ERR_ids.txt)

# 2. Query ENA backend to fetch the exact URL paths for your files
echo "Fetching FTP paths from ENA..."
curl -L "https://ebi.ac.uk{ERR_LIST}&fields=fastq_ftp&format=tsv" | tail -n +2 | tr ';' '\n' | sed 's/^/ftp:\/\//' > coral_download_urls.txt

# 3. Use wget to download all the files sequentially
echo "Starting raw FASTQ downloads via wget..."
wget -c -i coral_download_urls.txt

echo "Download job completed successfully at $(date)"
##################################################################
##################################################################
sbatch download.slurm













#
















#pairwise fst pipeline:

#first: need to calculate the SAF for each population individually (angsd)
  #eg: ./angsd -b list10  -anc hg19ancNoChr.fa -out pop1 -dosaf 1 -gl 1
#second: need to calculate the 2D SFS including all populations (angsd realsfs)
  #pairwise for each set of population pairs
  #eg: ./misc/realSFS pop1.saf.idx pop2.saf.idx -P 24 >pop1.pop2.ml
      #./misc/realSFS pop1.saf.idx pop3.saf.idx -P 24 >pop1.pop3.ml
#third: get fst index (angsd realsfs)
  #eg: /misc/realSFS fst index pop1.saf.idx pop2.saf.idx pop3.saf.idx -sfs pop1.pop2.ml -sfs pop1.pop3.ml -sfs pop2.pop3.ml -fstout here
#fourth: calculate fst (angsd realsfs)
  #.eg: ./misc/realSFS fst stats2 here.fst.idx -win 50000 -step 10000 >slidingwindow


































































































































































