# Mock mycobiome communities - metagenomic data generation and analysis


A collection of all scripts and resulting files related to ***"Challenges in capturing the mycobiome from shotgun metagenomic data: lack of software and databases"***.  

### Contents
 
- [Folder structure](#folder-structure)
- [Data analysis](#data-analysis)

### Folder structure
1. `./Mock_fungal/mock_conda`: 
    
    a) `snakemake_conda.yaml`: configuration file for conda environment with snakemake and all dependencies (ART/Kraken2/MetaPhlAn4) for mock data generation and Kraken2/MetaPhlAn4 taxonomy analysis

    b) `eukdetect_equal_{reads/coverage}.yaml`: configuration files used for taxonomy classification with EukDetect
2. `./Mock_fungal/config`: configuration file for the Snakemake
3. `./Mock_fungal/mock_profiles`: 

    a) `unicellular_classes.tsv`: List of taxIDs for unicellular classes within Ascomycota and Basidiomycota that were initially searched for

    b) `manually_added_taxa.csv`: List of fungal species manually added to the list of taxIDs

    c) `genomes_to_download.csv`: List of species IDs that were searched through the NCBI RefSeq database - generated using `./Mock_fungal/scripts/data_analysis/prepare_fastas.py`

    d) `final_genomes_summary.csv`: ***List of NBCI RefSeq genomes that were used for the mock communities**** - generated using `./Mock_fungal/scripts/data_analysis/prepare_fastas.py`

    e) `metadata_ufcg.csv`: Genome metadata file prepared for the UFCG phylogenetic tree generation

    f) `/profiles`: ***Equal reads mock community profiles*** - generated using `./Mock_fungal/scripts/data_analysis/prepare_fastas.py`

    g) `/profiles_equal_cov`: ***Equal coverage mock community profiles*** - generated using `./Mock_fungal/scripts/data_analysis/prepare_fastas.py`

4. `./Mock_fungal/workflow`: Scripts for mock communities sequencing data generation, taxonomy classification, and data analysis. Detailed info is provided in [Data analysis](#data-analysis).

5. `./Mock_fungal/data`: All figures and tables generated during the analysis of mock communities data

    a)`/kraken` - Kraken2 reports and their summary

    b)`/metaphlan` - MetaPhlAn4 reports and their summary

    c) `/eukdetect` - EukDetect reports and their summary

    d)`/HMS` - MycobiomeScan2.0 reports and their summary

    e)`/funomic` - FunOMIC reports and their summary

    f)`/micop` - MiCoP reports and their summary

    g)`/all_methods_summary` - Statistics, summary reports and figures 


### Data analysis

All python scripts were executed using Spyder IDE v5.5.4 (Python v3.12). Scripts are located in `./Mock_fungal/workflow/scripts` and `./Mock_fungal/workflow/scripts/data_analysis`. 

1. `prepare_fastas.py` - Download species taxids based on class taxids; download existing NCBI RefSeq genomes; concatenate contigs in each fasta file; and create profiles for the mock communities. Note that the script is dependent on the [ncbi-datasets](https://github.com/ncbi/datasets) conda environment

2. Reads simulation by ART, mock community metagenome data generation and its classification by Kraken2 and MetaPhlAn4 are wrapped into the Snakemake v8, conda is provided in `./Mock_fungal/mock_conda/snakemake_conda.yaml`. When all fastas, profiles and Kraken2/MetaPhlAn4 databases are available, provide paths to data in `./Mock_fungal/config/config.yaml` and run 
```bash
    cd ./Mock_fungal
    snakemake
```

3. The EukDetect and the HumanMycobiomeScan (MycobiomeScan v2.0) classifications are performed independently using the original pipeline/commands. Configuration scripts for EukDetect are located in `./Mock_fungal/mock_conda/eukdetect_equal_{reads/coverage}.yaml`

4. `{kraken/metaphlan/hms/eukdetect/funomic/micop}_summary.py` - Summarize taxonomy predictions on species/genus and family levels; find true and false positives. These summaries will be used further for all the analyses

5. `check_presence_database.py` - Find which mock community species are deposited in the databases; make a summary of genomes characteristics (*Figure 1*)

6. `abundance_estimation.py` - Calculate relative abundance on different taxonomy levels; calculate RMSE; compare equal reads vs equal coverage predictions; generate boxplot (*Figure 2b*)

7. `summarize_prec_recall.py` - Calculate precision, recall and F1 score for all tools; compare equal reads vs equal coverage predictions; generate boxplot (*Figure 2a*); Calculate Pearson correlation to community richness (*Figure 2c*). Note that for this part, RMSE results from `abundance_estimation.py` should be available

8. `create_datasets_itol.py` - Generate files for visualization of phylogenetic tree with [iTOL](https://itol.embl.de) (used for *Figure 3a*)

9. `taxa_detection.py` - Find which species/genera/families were differentially detected between equal reads and equal coverage mock communities by tested tools; which genera are not detected in case there are more than 1 species per genus (*Figure 3b*); presense of unidentified genera in the databases




