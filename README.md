# Challenges in capturing the mycobiome from shotgun metagenomic data

## This is a collection of all scripts and resulting files related to a publication [Avershina et al., 2024](INSERT_LINK_HERE)


## Contents

- [Folder structure](#folder-structure)
- [Data analysis](#data-analysis)

### Folder structure
1. `./mock_conda`: 
    
    a) snakemake_conda.yaml: configuration file for conda environment with snakemake and all dependencies (ART/Kraken2/MetaPhlAn4) for mock data generation and Kraken2/MetaPhlAn4 taxonomy analysis

    b) eukdetect_equal_{reads/coverage}.yaml: configuration files used for taxonomy classification with EukDetect
2. `./config`: configuration file for the Snakemake
3. `./mock_profiles`: 

    a) unicellular_classes.tsv: List of taxIDs for unicellular classes within Acomycota and Basidiomycota that were initially searched for

    b) manually_added_taxa.csv: List of fungal species manually added to the list of taxIDs

    c) genomes_to_download.csv: List of species IDs that were searched throuh the NCBI RefSeq database - generated using ./scripts/data_analysis/prepare_fastas.py

    d) final_genomes_summary.csv: List of species IDs that had genome assembly deposited in the NCBI RefSeq - generated using ./scripts/data_analysis/prepare_fastas.py

    e) metadata_ufcg.csv: Genome metadata file prepared for the UFCG phylogenetic tree generation

    f) `./profiles`: Equal reads mock community profiles - generated using ./scripts/data_analysis/prepare_fastas.py

    g) `./profiles_equal_cov`: Equal coverage mock community profiles - generated using ./scripts/data_analysis/prepare_fastas.py

4. `./workflow`: Snakefile for mock communities sequencing data generation, taxonomy classification and scripts related to data analysis. The description of the scripts is provided in [Data analysis](#data-analysis).

5. `./data`: All figures and tables generated during the analysis of mock communities data

    a)`./data/kraken` - Kraken2 reports and their summary

    b)`./data/metaphlan` - MetaPhlAn4 reports and their summary

    c) `./data/eukdetect` - EukDetect reports and their summary

    d)`./data/HMS` - MycobiomeScan2.0 reports and their summary

    e)`./data/all_methods-summary` - Statistics, summary reports and figures 


### Data analysis



