# Mock mycobiome communities - metagenomic data generation and analysis

### This is a collection of all scripts and resulting files related to "Challenges in capturing the mycobiome from shotgun metagenomic data" publication by [Avershina et al., 2024](INSERT_LINK_HERE)

Fastq files generated here are deposited in European Nucleotide Archive under [PRJEB79245](https://www.ebi.ac.uk/ena/browser/home)

### Contents

- [Folder structure](#folder-structure)
- [Data analysis](#data-analysis)

### Folder structure
1. `./mock_conda`: 
    
    a) `snakemake_conda.yaml`: configuration file for conda environment with snakemake and all dependencies (ART/Kraken2/MetaPhlAn4) for mock data generation and Kraken2/MetaPhlAn4 taxonomy analysis

    b) `eukdetect_equal_{reads/coverage}.yaml`: configuration files used for taxonomy classification with EukDetect
2. `./config`: configuration file for the Snakemake
3. `./mock_profiles`: 

    a) `unicellular_classes.tsv`: List of taxIDs for unicellular classes within Ascomycota and Basidiomycota that were initially searched for

    b) `manually_added_taxa.csv`: List of fungal species manually added to the list of taxIDs

    c) `genomes_to_download.csv`: List of species IDs that were searched throuh the NCBI RefSeq database - generated using ./scripts/data_analysis/prepare_fastas.py

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

All other analyses were performed in Python v3.12 using Spyder IDE v5.5.4 unless stated otherwise.
All scripts are located in `./workflow/scripts/data_analysis`. 

1. `prepare_fastas.py` - Download allspecies taxids, download NCBI RefSeq genomes, concatenate genome contigs in each file, and create profiles for the mock communities. Note that the script is dependent on the [ncbi-datasets](https://github.com/ncbi/datasets) conda environment

2. Reads simulation by ART, mock community metagenome data generation and its classification by Kraken2 and MetaPhlAn4 are wrapped into the Snakemake v8, conda is provided in `./mock_conda/snakemake_conda.yaml`. After all fastas, profiles and Kraken2/MetaPhlAn4 databases are available, provide paths to data in `./config/config.yaml` and run 
```bash
    cd ./Mock_fungal
    snakemake
```
The EukDetect and the HumanMycobiomeScan (MycobiomeScan v2.0) classifications are performed independently using the original pipeline/script. Configuration scripts used for EukDetect are located in `./mock_conda/eukdetect_equal_{reads/coverage}.yaml`

3. `{kraken/metaphlan/hms/eukdetect}_summary.py` - Summarize taxonomy predictions on species/genus and family levels; find true and false positives. These summaries will be used further for all the analyses

4. `abundance_estimation.py` - Calculate relative abundance on different taxonomy levels; calculate RMSE; compare equal reads vs equal coverage predictions; generate boxplot (Figure 2b)

5. `summarize_prec_recall.py` - Calculate precision, recall and F1 score for all tools; compare equal reads vs equal coverage predictions; generate boxplot (Figure 2a); Calculate Pearson correlation to community richness (Figure 2c). Note that for this part, RMSE results from `abundance_estimation.py` should be available

6. `check_presence_database.py` - Find which mock community species are deposited in the databases; make a summary of genomes characteristics (Figure 1)

7. `ER_EC_detection.py` - Find which species/genera/families were differentially detected between equal reads and equal coverage mock communities by tested tools

8. `create_datasets_itol.py` - Generate files for visualization of phylogenetic tree with iTOL

