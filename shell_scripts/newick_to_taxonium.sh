# Convert the Newick tree to Taxonium format with metadata info
inputTree="/nfs/research/goldman/anoufa/data/MAPLE_output/masked/mas_tree.tree"
outputTree="/nfs/research/goldman/anoufa/data/MAPLE_output/masked/mas_tree.jsonl"
metadataFile="/nfs/research/goldman/anoufa/data/MAPLE_input/others/viridian_samples.metadata.tsv"

sbatch -J ntt -t 04:00:00 --mem=32G \
    -o /nfs/research/goldman/anoufa/data/out_err/ntt.out \
    -e /nfs/research/goldman/anoufa/data/out_err/ntt.err \
    --wrap="/homes/anoufa/.pyenv/versions/3.11.6/bin/python3 -m taxoniumtools.newick_to_taxonium \
        -i $inputTree \
        -m $metadataFile \
        -o $outputTree \
        -c Country,Viridian_pangolin_1.29 \
        --key_column Run"