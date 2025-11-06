
minimap2="/hps/software/users/goldman/anoufa/minimap2-2.29_x64-linux/minimap2"
samtools="/hps/software/users/goldman/samtools/samtools-1.13/samtools"

read_name="ERR4462326"
ref_path="/nfs/research/goldman/anoufa/data/NC_045512.2.fasta"
read_path_1="/nfs/research/goldman/anoufa/data/nicola_MNMs_paper/${read_name}_1.fastq.gz"
read_path_2="/nfs/research/goldman/anoufa/data/nicola_MNMs_paper/${read_name}_2.fastq.gz"
bam_path="/nfs/research/goldman/anoufa/data/nicola_MNMs_paper/$read_name.bam"

# First option

$minimap2 -t 4 -a -x sr "$ref_path" "$read_path_1" "$read_path_2" | \
$samtools fixmate -u -m - - | \
$samtools sort -T /tmp/sorting - | \
$samtools markdup -@8 --reference "$ref_path" - "$bam_path"

$samtools index $bam_path

$samtools view -c -F 4 $bam_path


# Second option

# vdn_cmd_mhunt="singularity exec /nfs/research/zi/mhunt/Containers/viridian_v1.5.1.img viridian run_one_sample --run_accession ERR8959214 --keep_bam --outdir /nfs/research/goldman/anoufa/data/nicola_MNMs_paper/vdn_1"
# vdn_cmd="singularity exec /nfs/research/zi/mhunt/Containers/viridian_v1.5.1.img viridian run_one_sample --run_accession ERR4462326 --keep_bam --outdir /nfs/research/goldman/anoufa/data/nicola_MNMs_paper/vdn"

# sbatch -J vdn_p_$read_name -t 00:30:00 --mem=16G \
#     -o /nfs/research/goldman/anoufa/data/nicola_MNMs_paper/vdn.out \
#     -e /nfs/research/goldman/anoufa/data/nicola_MNMs_paper/vdn.err \
#     --wrap="$vdn_cmd_mhunt"