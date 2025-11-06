#!/bin/bash
#SBATCH --job-name=map_reads
#SBATCH --output=/nfs/research/goldman/anoufa/data/samples/z_logs/%x_%j.out
#SBATCH --error=/nfs/research/goldman/anoufa/data/samples/z_logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# --- User input: list of sample names ---
sample_list="/nfs/research/goldman/anoufa/data/samples/z_logs/sample_names.txt"

# --- Tools ---
minimap2="/hps/software/users/goldman/anoufa/minimap2-2.29_x64-linux/minimap2"
samtools="/hps/software/users/goldman/samtools/samtools-1.13/samtools"
ref_path="/nfs/research/goldman/anoufa/data/MAPLE_input/MAPLE_REF.fasta"

# --- Loop over sample names ---
while read -r sample_name; do
    echo "Processing sample: $sample_name"

    data_dir="/nfs/research/goldman/anoufa/data/samples/${sample_name}"
    read_path_1="${data_dir}/${sample_name}_1.fastq.gz"
    read_path_2="${data_dir}/${sample_name}_2.fastq.gz"
    single_path="${data_dir}/${sample_name}.fastq.gz"
    bam_path="${data_dir}/${sample_name}.bam"
    bai_path="${data_dir}/${sample_name}.bam.bai"

    # --- Skip if BAM and BAI already exist ---
    if [[ -f "$bam_path" && -f "$bai_path" ]]; then
        echo "BAM and BAI already exist for $sample_name. Skipping..."
        echo "---------------------------------------"
        continue
    fi

    echo "Checking read files..."
    if [[ -f "$read_path_2" ]]; then
        echo "Paired-end data detected."
        
        # Count reads
        reads_1=$(zcat "$read_path_1" | wc -l)
        reads_2=$(zcat "$read_path_2" | wc -l)
        echo "Reads in $read_path_1: $((reads_1 / 4))"
        echo "Reads in $read_path_2: $((reads_2 / 4))"

        # Run PE pipeline
        "$minimap2" -t 8 -a -x sr "$ref_path" "$read_path_1" "$read_path_2" | \
        "$samtools" fixmate -u -m - - | \
        "$samtools" sort -T /tmp/sorting -@4 - | \
        "$samtools" markdup -@8 --reference "$ref_path" - "$bam_path"

    elif [[ -f "$single_path" ]]; then
        echo "Single-end data detected."

        # Count reads
        reads=$(zcat "$single_path" | wc -l)
        echo "Reads in $single_path: $((reads / 4))"

        # Run SE pipeline
        "$minimap2" -t 8 -a -x map-ont "$ref_path" "$single_path" | \
        "$samtools" sort -T /tmp/sorting -@4 - | \
        "$samtools" markdup -@8 --reference "$ref_path" - "$bam_path"

    else
        echo "No input FASTQ files found for $sample_name"
        continue
    fi

    # Index BAM
    "$samtools" index "$bam_path"

    # Count mapped reads
    echo "Mapped reads:"
    "$samtools" view -c -F 4 -F 256 -F 2048 "$bam_path"

    # Flagstat
    echo "Flagstat output:"
    "$samtools" flagstat "$bam_path"

    echo "Sample $sample_name done."
    echo "---------------------------------------"

done < "$sample_list"

echo "All jobs done."
