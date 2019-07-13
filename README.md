# The 4-way paper repository
Scripts used in the paper "Host-Microbe-Drug-Nutrient Screen Identifies Bacterial Effectors of Metformin Therapy".

## 4-way screen analysis (4-way_screen/)
- **Biolog_analysis_NGM_Ecoli.R** - 4-way screen analysis for *E. coli* (Fig. 1, S1)
- **Biolog_analysis_NGM_Worm_Imaging.R** - 4-way screen analysis for *C. elegans* (Fig. 1, S1)
- **Biolog_analysis_NGM_combined.R** - Join *C. elegans* and *E. coli* 4-way screen analysis and final figures (Fig. 1, S1)

## *E. coli* growth analysis (bacterial_growth/)
- **BacGrowth_media.R** - Bacterial growth in different media (Fig. 1, S1)
- **BacGrowth_TF.R** - Bacterial growth in TF knockouts (Fig. S3)
- **BacGrowth_PTS.R** - Bacterial growth in PTS knockouts (Fig. S4)
- **BacGrowth_Arg_mutants.R** - Bacterial growth in arginine degradation pathway knockouts (Fig. S5)
- **BacGrowth_acetoacetate.R** - Bacterial growth acetoacetate supplementation (Fig. S7)

## *C. elegans* fluorescence analysis (fluorescence/)
- **Fluorescence_collect.R** - Collect nematode fluorescence data and perform global normalisation. Other scripts rely on its output.
- **Fluorescence_media.R** - Nematode fluorescence in different media, ribose and glycerol supplementation, Pcrp (oe) (Fig. 1, S1, S2, 3, S4)
- **Fluorescence_TF_acs2.R** - Nematode fluorescence with bacterial TF knockouts (Fig. 2, S2, S3)
- **Fluorescence_glycerol.R** - Nematode fluorescence in glycerol supplementation (Fig. S2)
- **Fluorescence_arg_agm.R** - Nematode fluorescence with arginine and agmatine pathway bacterial knockouts (Fig. 4, S5)
- **Fluorescence_transgenes.R** - Nematode fluorescence for various transgenes (depends on **Fluorescence_collect.R**) (Fig. 6, S6)
- **Fluorescence_vha6.R** - Nematode fluorescence analysis for vha-6 reporter (Fig. 6, S6)
- **Fluorescence_acetoacetate.R** - Nematode fluorescence in acetoacetate supplementation (Fig. S7)

## *E. coli* and *C. elegans* metabolomics analysis (metabolomics/)
- **Metabolomics_Ecoli_HMT.R** - *E. coli* metabolomics (Fig. 2, 4, S5)
- **Metabolomics_Celegans_FA.R** - *C. elegans* fatty acid metabolomics (Fig. S7)

## *E. coli* proteomics analysis (proteomics/)
- **Proteomics_Ecoli.R** - Proteomics dataset analysis (Fig. 2)
- **Proteomics_TF.R** - TF enrichment (Fig. 2)

##  *C. elegans* RNA-Seq analysis (rna-seq/)
- **RNAseq_Celegans.R** - *C. elegans* N2 nematode RNA-Seq analysis (Fig. 6)
