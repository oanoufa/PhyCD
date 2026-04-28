# PhyCD: Phylogeny-aided Contamination Detection

PhyCD is a computational approach for detecting consensus sequence errors
caused by the combination of sample contamination and amplicon dropout in
SARS-CoV-2 (and potentially other amplicon-sequenced pathogen) genomes.

PhyCD masks consensus genome positions associated with suspicious
sequencing read coverage drops, then leverages pandemic-scale
phylogenetic placement to identify putative contamination events.

## Overview

Contamination occurs when multiple pathogen genomes from different hosts
are accidentally mixed within one sequenced sample. When combined with
amplicon dropout due to mutations at primer binding sites that prevent
amplification of the majority genome, the minority contaminant's reads
can dominate at dropout regions, producing a hybrid consensus sequence
that resembles recombination and introduces artefactual mutations.

PhyCD addresses this by:

- Masking genome regions exhibiting suspicious relative drops in
  sequencing coverage
- Phylogenetically placing masked and unmasked sequences onto a
  pandemic-scale reference tree
- Identifying samples where masking substantially reduces phylogenetic
  placement distance, indicating that the masked positions likely
  contained artefactual mutations
- Searching for plausible contaminant genomes and statistically
  assessing contamination hypotheses

Applied to nearly 5 million SARS-CoV-2 genomes, PhyCD identified 10,942
putative contamination events under conservative parameters, twice as
many as a randomly masked control baseline. This corresponds to the removal of 64,753 potentially
artefactual substitutions from phylogenetic placement branch lengths.

## Pipeline

The pipeline consists of the following steps:

1. **Filtering**: Samples with few heterozygous sites and no large
   coverage drops are selected to build a reference phylogenetic tree
   using [MAPLE](https://github.com/NicolaDM/MAPLE).

2. **Dropout masking**: For each unfiltered sample, genome positions with
   coverage below a threshold fraction (ρ) of the sample median
   coverage (or below 50X) are masked. Three consensus sequences are
   generated per sample:
   - *Dropout unmasked*: no dropout masking applied
   - *Dropout masked*: positions below ρ * median coverage are masked
   - *Randomly masked*: same number of positions masked, but shifted
     randomly along the genome (control baseline)

3. **Phylogenetic placement**: All three sequences are placed onto the
   reference tree using MAPLE. Samples where the *dropout masked* sequence
   has a placement branch length below 5 and 2 or more substitutions
   shorter than the *dropout unmasked* sequence are flagged as putatively
   contaminated.

4. **Contaminant search**: For each flagged sample, plausible contaminant
   genomes are identified from the filtered dataset based on mutation
   matching at *dropout masked* regions and minor allele matching at
   heterozygous sites.

5. **Likelihood-based assessment**: A probabilistic model (adapted from
   [Eyre et al., 2013](https://doi.org/10.1371/journal.pcbi.1003059))
   statistically evaluates each sample-contaminant pair to estimate a
   posterior probability of contamination.

## Installation

    git clone https://github.com/oanoufa/PhyCD.git
    cd PhyCD
    python -m venv .phycd_env
    source .phycd_env/bin/activate
    pip install -r pyenv_requirements.txt

### Dependencies

- Python >= 3.8
- [MAPLE](https://github.com/NicolaDM/MAPLE) v0.7.5
- [Snakemake](https://snakemake.github.io/) v8.5.2

## Usage

### Running the pipeline

PhyCD can be run end-to-end using the provided `run_phycd.sh` script, which processes all samples in a single run. The provided data is a small example dataset made for testing purposes only. It is composed of 300 samples randomly drawn within lineage B.1.1.7.

PhyCD needs pandemic-scale data and therefore does not produce meaningful results on this example dataset. For larger datasets, the pipeline should be ran on clusters.

Before running the pipeline, the `root_dir` variable in the `Snakefile` should be set to the path to your local copy of the PhyCD repository. Parameters are set for the example pipeline run, but can be modified as needed (see below).

    bash run_phycd.sh

### Parameters

| Parameter | Symbol | Default | Description |
|-----------|--------|---------|-------------|
| Dropout masking threshold | ρ | 0.05 | Fraction of median coverage below which positions are masked |
| Coverage filter threshold | κ | 0 | Max positions below ρ*M allowed in filtered samples |
| Heterozygosity threshold | η | 0.10 | Minor allele proportion above which a site is heterozygous |
| Heterozygous sites filter | θ | 3 | Max heterozygous sites allowed in filtered samples |
| Branch length threshold | - | 5 | Max placement branch length for flagged samples |
| Distance difference threshold | - | 2 | Min reduction in placement branch length to flag a sample |

Lower values of ρ (e.g., 5%) offer the highest signal-to-noise ratio
and masking efficiency. Higher values (e.g., 20%) flag more samples but
introduce proportionally more false positives. See the manuscript for a
detailed parameter exploration.

## Input data

PhyCD expects as input:

- **Quality control files** assembled by
  [Viridian](https://github.com/iqbal-lab-org/viridian)
- **Reference sequence** of the pathogen (e.g., SARS-CoV-2 reference genome)

## Output

PhyCD produces:

- `PhyCD/pipeline/{parameters}/data/8/masked_samples_{parameters}.tsv`: The tsv file containing the flagged putatively contaminated samples with various metrics, candidate contaminant genomes, posterior probability of contamination, initial and final placements;
- `PhyCD/pipeline/{parameters}/data/2/alignment_files`: [MAPLE](https://github.com/NicolaDM/MAPLE) alignment files for the filtered samples used to build the reference tree, the masked, unmasked, and randomly masked sequences for all samples except the filtered ones;
- `PhyCD/pipeline/{parameters}/data/clean_tree/clean_{parameters}_tree.jsonl`: The reference phylogenetic tree built from the filtered samples with country and lineage annotations. The tree can be visualized on [Taxonium](https://www.taxonium.org/);
- `PhyCD/pipeline/{parameters}/data/figs`: various figures summarizing the results.

## Reproducing the manuscript results

The results in the manuscript were generated using the
[Viridian](https://github.com/iqbal-lab-org/viridian) dataset of
4,952,451 publicly shared SARS-CoV-2 sequenced samples
([Hunt et al., 2026](https://doi.org/10.1038/s41592-025-02947-1)).

A Zenodo repository is publicly available at [Zenodo repository](https://zenodo.org/records/19815950). The repository contains:
- the reference phylogenetic tree of 785,011 samples;
- the *dropout masked* alignment augmented with the filtered samples (4,952,451 samples + reference);
- the list of the 10,942 flagged putatively contaminated samples.
All those files were generated using PhyCD with default parameters (0.05, 0, 0.10, 3).

## License

This project is licensed under the GNU General Public License v3.0.
See the [LICENSE](LICENSE) file for details.

## Funding

This work was supported by the European Molecular Biology Laboratory - European Bioinformatics Institute (EMBL-EBI) and the French Embassy in London through their internship program with EMBL-EBI.
