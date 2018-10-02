# Metformin CRP Project

## *E. coli* growth analysis
- BacGrowth_Filipe_D1_D5_D8.R - Bacterial growth analysis for figure 1 (growth of bacteria from gut)
- BacGrowth_GutSize_Figure1.R - Bacterial growth analysis for figure 1 (Gut size changes, data extracted from PPT slides)
- BacGrowth_Gut_Figure1.R - Bacterial growth analysis for figure 1
- BacGrowth_Media.R - Bacterial growth in different media
- BacGrowth_Rosie_Acetoacetate.R - Bacterial growth analysis for figure 4 (Acetoacetate supplementation)
- BacGrowth_Rosie_Arg_mutants.R - Bacterial growth analysis for figure 7 (Arginine degradation pathway knockouts)
- BacGrowth_Rosie_IPTG_oeCRP_cyaA.R - Bacterial growth for figure 6 (CRP overexpression)
- BacGrowth_Rosie_PTS.R - Bacterial growth analysis for figure 6 (PTS knockouts)
- BacGrowth_Rosie_Resistance.R - Bacterial growth analysis for figure 1 (OP50-MR, OP50-R1, OP50-R2)
- BacGrowth_Rosie_TF.R - Bacterial growth analysis for figure 5 (TF knockouts)

## 4-way screen analysis
- Biolog_analysis_NGM_Ecoli.R - 4-way screen analysis for E. coli
- Biolog_analysis_NGM_Worm_Imaging.R - 4-way screen analysis for C. elegans
- Biolog_analysis_NGM_combined.R - 4-way screen analysis joint and final figures
- Metformin_Worm_Imaging_training.R - Analysis of nematode fluorescence distributions for adaptive thresholding

## *C. elegans* fluorescence analysis
- Fluorescence_Glycerol.R - Nematode fluorescence analysis for glycerol supplementation
- Fluorescence_TF_acs2.R - Nematode fluorescence analysis with bacterial TF knockouts
- Fluorescence_collect.R - Nematode fluorescence data collect. Performs global normalisation of data. Multiple other scripts relie on its output.
- Fluorescence_collect_Acetoacetate_transgenes_Figure4.R - Nematode fluorescence collect and analysis for figure 4 (Acetoacetate supplementation)
- Fluorescence_collect_ArgAgm_Figure7.R - Nematode fluorescence analysis for figure 7
- Fluorescence_figures_redo.R - Nematode fluorescence re-analysis for various one-off figures
- Fluorescence_transgenes.R - Nematode fluorescence analysis for all transgenes (depends on Fluorescence_collect.R)
- Fluorescence_vha6_Fig5-6.R - Nematode fluorescence analysis for vha-6 reporter

## *E. coli* transcriptional reporter library screen
- GFP_reporter_NGM.R - Old version
- GFP_reporter_NGM_new.R - Tidyverse version (not complete due to clashes between packages)
- UAL_reporters_Volcano_plots.R - Generate volcano plots from results

## *C. elegans* lifespan statistical analysis in R
- Lifespans.R - C. elegans lifespan statistical analysis in R

## *E. coli* and *C. elegans* metabolomics analysis
- Metabolomics_All_comparison.R - Comparison of various metabolomics datasets
- Metabolomics_Celegans_AA.R - *C. elegans* amino acids
- Metabolomics_Celegans_AA_old.R - *C. elegans* amino acids
- Metabolomics_Celegans_FA.R - *C. elegans* fatty acids
- Metabolomics_Ecoli_100met.R - *E. coli* 100 metabolites dataset
- Metabolomics_Ecoli_100met_new.R - *E. coli* 100 metabolites dataset
- Metabolomics_Ecoli_AA.R - *E. coli* amino acids
- Metabolomics_Ecoli_AA_new.R - *E. coli* amino acids
- Metabolomics_Ecoli_HMT.R - *E. coli* Human Metabolome technologies
- Metabolomics_Ecoli_nucleotide.R - *E. coli* nucleotide metabolomics
- Metabolomics_Ecoli_nucleotide_New_extended.R - *E. coli* nucleotide metabolomics
- Metabolomics_Ecoli_nucleotide_metabolomics.R - *E. coli* nucleotide metabolomics

## *E. coli* proteomics analysis
- Proteomics_Ecoli.R - Data analysis
- Proteomics_Ecoli_update.R - Migration to tidyverse (not complete)
- Proteomics_TF.R - TF enrichment analysis

## *C. elegans* RNA-Seq analysis
- RNAseq-Metabolomics_Joint_PCA.R - Joint PCA, MDS plots for *C. elegans* metabolomics and RNA-Seq.
- RNAseq_Celegans.R - *C. elegans* Metformin RNA-Seq analysis
- RNAseq_Celegans2.R - *C. elegans* Metformin RNA-Seq analysis (experimentation with parameters, etc.)
- RNAseq_Celegans_Ballgown.R - RNA-Seq data analysis using Ballgown (not recommended, as EdgeR has more functionality)
- RNAseq_Celegans_DR.R - *Heintz et al. 2017* *C. elegans* DR dataset analysis
- RNAseq_Celegans_HTseq.R - RNA-Seq data analysis using HTSeq read counts
- RNAseq_Celegans_MTA.R - Metformin *C. elegans* RNA-Seq data MTA summary
- RNAseq_Celegans_Metf_and_DR_joint.R - *C. elegans* Metformin RNA-Seq and *Heintz et al. 2017* analysis results join
- RNAseq_Celegans_Rotenone.R - Analysis of *Schmeisser et al. 2013* *C. elegans* RNA-Seq data
- RNAseq_Celegans_Stringtie_count_prep.R - Preparation of read-counts from Stringtie
- RNAseq_Celegans_Stringtie_count_annot_prep.R - Preparation of read-counts from Stringtie and their quality control
- RNAseq_Celegans_permutations.R - Metformin *C. elegans* RNA-Seq re-analysis for MTA
- RNAseq_Tuxedo.sh - Commands for running Tuxedo suite (HiSat + StringTie)
- RNAseq_Tuxedo_Paired.sh - Commands for running Tuxedo suite with paired-end reads (*Heintz 2017*) (HiSat + StringTie)
- RNAseq_Tuxedo_Paired_part2.sh - Commands for running Tuxedo suite (*Heintz 2017*) (HiSat + StringTie)
