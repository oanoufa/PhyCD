#!/bin/bash
#SBATCH --job-name=map_read
#SBATCH --output=/nfs/research/goldman/anoufa/data/nicola_MNMs_paper/_logs/%x.out
#SBATCH --error=/nfs/research/goldman/anoufa/data/nicola_MNMs_paper/_logs/%x.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

read_name="SRR22130132"

# Define paths to tools and files
minimap2="/hps/software/users/goldman/anoufa/minimap2-2.29_x64-linux/minimap2"
samtools="/hps/software/users/goldman/samtools/samtools-1.13/samtools"

ref_path="/nfs/research/goldman/anoufa/data/NC_045512.2.fasta"
read_path_1="/nfs/research/goldman/anoufa/data/nicola_MNMs_paper/${read_name}/${read_name}_1.fastq.gz"
read_path_2="/nfs/research/goldman/anoufa/data/nicola_MNMs_paper/${read_name}/${read_name}_2.fastq.gz"
bam_path="/nfs/research/goldman/anoufa/data/nicola_MNMs_paper/${read_name}/${read_name}.bam"


# Count reads in FASTQ files
echo "Counting reads in input files..."
reads_1=$(zcat "$read_path_1" | wc -l)
reads_2=$(zcat "$read_path_2" | wc -l)
echo "Reads in $read_path_1: $((reads_1 / 4))"
echo "Reads in $read_path_2: $((reads_2 / 4))"

# Mapping + processing
"$minimap2" -t 8 -a -x sr "$ref_path" "$read_path_1" "$read_path_2" | \
"$samtools" fixmate -u -m - - | \
"$samtools" sort -T /tmp/sorting -@4 - | \
"$samtools" markdup -@8 --reference "$ref_path" - "$bam_path"

# Index BAM
"$samtools" index "$bam_path"

# Count mapped reads
echo "Mapped reads:"
"$samtools" view -c -F 4 -F 256 -F 2048 "$bam_path"

# Flagstat
echo "Flagstat output:"
"$samtools" flagstat "$bam_path"

# Job done
echo "Job done"