# "Host-Microbe-Drug-Nutrient Screen Identifies Bacterial Effectors of Metformin Therapy"

## bacterial_growth - *E. coli* growth analysis
- **BacGrowth_media.R** - Bacterial growth in different media (Fig. 1, S1)
- **BacGrowth_TF.R** - Bacterial growth in TF knockouts (Fig. S3)
- **BacGrowth_PTS.R** - Bacterial growth in PTS knockouts (Fig. S4)
- **BacGrowth_Arg_mutants.R** - Bacterial growth in arginine degradation pathway knockouts (Fig. S5)
- **BacGrowth_acetoacetate.R** - Bacterial growth acetoacetate supplementation (Fig. S7)

## 4-way_screen - 4-way screen analysis (Fig. 1, S1)
- **Biolog_analysis_NGM_Ecoli.R** - 4-way screen analysis for *E. coli*
- **Biolog_analysis_NGM_Worm_Imaging.R** - 4-way screen analysis for *C. elegans*
- **Biolog_analysis_NGM_combined.R** - Join *C. elegans* and *E. coli* 4-way screen analysis and final figures

## fluorescence - *C. elegans* fluorescence analysis
- **Fluorescence_collect.R** - Collect nematode fluorescence data and perform global normalisation. Other scripts rely on its output.
- **Fluorescence_media.R** - Nematode fluorescence in different media, ribose and glycerol supplementation, Pcrp (oe) (Fig. 1, S1, S2, 3, S4)
- **Fluorescence_TF_acs2.R** - Nematode fluorescence with bacterial TF knockouts (Fig. 2, S2, S3)
- **Fluorescence_glycerol.R** - Nematode fluorescence in glycerol supplementation (Fig. S2)
- **Fluorescence_arg_agm.R** - Nematode fluorescence with arginine and agmatine pathway bacterial knockouts (Fig. 4, S5)
- **Fluorescence_transgenes.R** - Nematode fluorescence for various transgenes (depends on **Fluorescence_collect.R**) (Fig. 6, S6)
- **Fluorescence_vha6.R** - Nematode fluorescence analysis for vha-6 reporter (Fig. 6, S6)
- **Fluorescence_acetoacetate.R** - Nematode fluorescence in acetoacetate supplementation (Fig. S7)

## metabolomics - *E. coli* and *C. elegans* metabolomics analysis
- **Metabolomics_Ecoli_HMT.R** - *E. coli* metabolomics (Fig. 2, 4, S5)
- **Metabolomics_Celegans_FA.R** - *C. elegans* fatty acid metabolomics (Fig. S7)

## proteomics - *E. coli* proteomics analysis
- **Proteomics_Ecoli.R** - Proteomics dataset analysis (Fig. 2)
- **Proteomics_TF.R** - TF enrichment (Fig. 2)

## rna-seq - *C. elegans* RNA-Seq analysis
- **RNAseq_Celegans.R** - *C. elegans* N2 nematode RNA-Seq analysis (Fig. 6)