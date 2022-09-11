# Preliminary processing and analysis for HIV/AS scRNA-seq dataset GSE198339

__Link to dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198339__

## Outline of this README
- Experiment information from GEO
- Analysis summary
- Software/R package information

## Experiment Information

| **Field** | **Information (Reproduced from GSE198339)** |
| ------------- | ------------- |
| **Title**	| Single-cell RNA sequencing of peripheral blood mononuclear cells reveals immune signaling dysregulations in people living with HIV and atherosclerosis |
| **Organism** | *Homo sapiens* |
| **Experiment type**	| Expression profiling by high throughput sequencing |
| **Summary** | **Background:** Atherosclerosis (AS)-associated cardiovascular disease is an important cause of mortality in an aging population of people living with HIV (PLWH). This elevated risk of atherosclerosis has been attributed to viral infection, prolonged usage of anti-retroviral therapy, and subsequent chronic inflammation.<br />**Methods:** To investigate dysregulated immune signaling in PLWH with and without AS, we sequenced 9368 peripheral blood mononuclear cells (PBMCs) from 8 PLWH, 4 of whom also had atherosclerosis (AS+).  To develop executable models of signaling pathways that drive cellular states in HIV-associated atherosclerosis, we developed the single-cell Boolean Omics Network Invariant Time Analysis (scBONITA) algorithm. ScBONITA (a) uses single-cell RNA sequencing data to infer Boolean rules for topologically characterized networks, (b) prioritizes genes based on their impact on signaling, (c) performs pathway analysis, and (d) maps sequenced cells to characteristic signaling states. We used scBONITA to identify dysregulated pathways in different cell-types from AS+ PLWH and AS- PLWH. To compare our findings with pathways associated with HIV infection, we used scBONITA to analyze a publicly available dataset of PBMCs from subjects before and after HIV infection. Additionally, the executable Boolean models characterized by scBONITA were used to analyze observed cellular states corresponding to the steady states of signaling pathways.<br /> **Results:** We identified an increased subpopulation of CD8+ T cells and a decreased subpopulation of monocytes in AS+ PLWH. Dynamic modeling of signaling pathways and pathway analysis with scBONITA provided a new perspective on the mechanisms of HIV-associated atherosclerosis. Lipid metabolism and cell migration pathways are induced by AS rather than by HIV infection. These pathways included AGE-RAGE and PI3K-AKT signaling in CD8+ T cells, and glucagon and cAMP signaling pathways in monocytes. Further analysis of other cell subpopulations suggests that the highly interconnected PI3K-AKT signaling pathway drives cell migratory state in response to dyslipidemia. scBONITA attractor analysis mapped cells to pathway-specific signaling states that correspond to distinct cellular states.<br />**Conclusions**: Dynamic network modeling and pathway analysis with scBONITA indicates that dysregulated lipid signaling regulates cell migration into the vascular endothelium in AS+ PLWH. Attractor analysis with scBONITA facilitated pathway-based characterization of cellular states that are not apparent in gene expression analyses.| 
| **Overall design**	| scRNA sequencing data for peripheral blood mononuclear cell samples from male HIV+ donors with and without atherosclerosis ***Raw data not available due to patient privacy concerns*** |
| **Contributor(s)**	| Palshikar MG, Palli R, Tyrell A, Maggirwar S, Schifitto G, Singh MV, Thakar J |

## Analysis summary:

### See the R notebook HIV_AS_dataset_analysis.rmd for code and all parameters

- Initialize Seurat object with raw counts and metadata (ie, clinical covariates). 
- Perform basic filtering and quality control, followed by dimensionality reduction and clustering. 
- Identify differentially populated cell clusters, ie, which clusters have disproportionate numbers of AS+ or AS- cells. 
- Identify cell cluster markers for each cluster, then characterize these markers using GSEA, over-representation analysis, and enricher with KEGG gene sets downloaded from MSigDB. 
- Similarly, perform GSEA and Enricher on genes DE between AS+ and AS- participants, for each subpopulation separately. 
- Generate figures similar to those in our publication. 

## Software/R package information:

Please note the version numbers of all the packages.

```

R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22000)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] plyr_1.8.6             knitr_1.29             circlize_0.4.10       
 [4] ComplexHeatmap_2.4.3   enrichplot_1.8.1       scales_1.1.1          
 [7] msigdbr_7.2.1          clusterProfiler_3.16.1 ggpubr_0.4.0          
[10] openxlsx_4.1.5         RColorBrewer_1.1-2     gplots_3.0.4          
[13] wesanderson_0.3.6      reshape_0.8.8          dplyr_1.0.9           
[16] Seurat_2.3.4           Matrix_1.2-18          cowplot_1.0.0         
[19] ggplot2_3.3.2         

loaded via a namespace (and not attached):
  [1] prabclus_2.3-2       R.methodsS3_1.8.1    tidyr_1.1.2         
  [4] bit64_4.0.2          irlba_2.3.3          multcomp_1.4-15     
  [7] R.utils_2.10.1       data.table_1.13.0    rpart_4.1-15        
 [10] generics_0.1.0       metap_1.8            snow_0.4-3          
 [13] BiocGenerics_0.34.0  TH.data_1.0-10       RSQLite_2.2.0       
 [16] RANN_2.6.1           europepmc_0.4        proxy_0.4-24        
 [19] bit_4.0.4            mutoss_0.1-12        xml2_1.3.2          
 [22] assertthat_0.2.1     viridis_0.5.1        xfun_0.22           
 [25] hms_0.5.3            evaluate_0.14        DEoptimR_1.0-8      
 [28] fansi_0.4.1          progress_1.2.2       caTools_1.18.0      
 [31] readxl_1.3.1         igraph_1.3.0         DBI_1.1.0           
 [34] tmvnsim_1.0-2        htmlwidgets_1.5.2    stats4_4.0.2        
 [37] purrr_0.3.4          ellipsis_0.3.2       backports_1.1.9     
 [40] gbRd_0.4-11          vctrs_0.4.1          Biobase_2.48.0      
 [43] ROCR_1.0-11          abind_1.4-5          withr_2.2.0         
 [46] ggforce_0.3.2        triebeard_0.3.0      robustbase_0.93-6   
 [49] checkmate_2.0.0      prettyunits_1.1.1    mclust_5.4.7        
 [52] mnormt_2.0.2         cluster_2.1.0        DOSE_3.14.0         
 [55] ape_5.4-1            segmented_1.3-0      crayon_1.3.4        
 [58] hdf5r_1.3.3          labeling_0.3         pkgconfig_2.0.3     
 [61] qqconf_1.0.0         tweenr_1.0.1         nlme_3.1-149        
 [64] nnet_7.3-14          rlang_1.0.4          diptest_0.75-7      
 [67] lifecycle_1.0.1      sandwich_3.0-0       downloader_0.4      
 [70] mathjaxr_1.0-1       doSNOW_1.0.19        cellranger_1.1.0    
 [73] polyclip_1.10-0      lmtest_0.9-38        urltools_1.7.3      
 [76] carData_3.0-4        zoo_1.8-8            base64enc_0.1-3     
 [79] ggridges_0.5.2       GlobalOptions_0.1.2  png_0.1-7           
 [82] viridisLite_0.3.0    rjson_0.2.20         bitops_1.0-6        
 [85] R.oo_1.24.0          KernSmooth_2.23-17   blob_1.2.1          
 [88] shape_1.4.4          lars_1.2             stringr_1.4.0       
 [91] qvalue_2.20.0        jpeg_0.1-8.1         rstatix_0.6.0       
 [94] gridGraphics_0.5-0   S4Vectors_0.26.1     ggsignif_0.6.0      
 [97] memoise_1.1.0        magrittr_1.5         ica_1.0-2           
[100] gdata_2.18.0         compiler_4.0.2       scatterpie_0.1.4    
[103] clue_0.3-57          plotrix_3.7-8        fitdistrplus_1.1-1  
[106] cli_3.3.0            dtw_1.22-3           pbapply_1.4-3       
[109] htmlTable_2.1.0      Formula_1.2-4        MASS_7.3-52         
[112] tidyselect_1.1.2     stringi_1.4.6        forcats_0.5.0       
[115] highr_0.8            yaml_2.2.1           GOSemSim_2.14.1     
[118] latticeExtra_0.6-29  ggrepel_0.8.2        fastmatch_1.1-0     
[121] tools_4.0.2          parallel_4.0.2       rio_0.5.16          
[124] rstudioapi_0.11      foreach_1.5.1        foreign_0.8-80      
[127] gridExtra_2.3        farver_2.0.3         Rtsne_0.15          
[130] ggraph_2.0.3         digest_0.6.25        rvcheck_0.1.8       
[133] BiocManager_1.30.10  fpc_2.2-8            Rcpp_1.0.8.3        
[136] car_3.0-9            broom_0.7.0          SDMTools_1.1-221.2  
[139] httr_1.4.2           AnnotationDbi_1.50.3 kernlab_0.9-29      
[142] Rdpack_2.1           colorspace_1.4-1     reticulate_1.16     
[145] IRanges_2.22.2       splines_4.0.2        yulab.utils_0.0.4   
[148] sn_1.6-2             graphlayouts_0.7.0   multtest_2.44.0     
[151] ggplotify_0.1.0      flexmix_2.3-17       jsonlite_1.7.1      
[154] tidygraph_1.2.0      modeltools_0.2-23    R6_2.4.1            
[157] TFisher_0.2.0        Hmisc_4.4-1          pillar_1.8.0        
[160] htmltools_0.5.3      glue_1.6.2           fastmap_1.1.0       
[163] BiocParallel_1.22.0  class_7.3-17         codetools_0.2-16    
[166] fgsea_1.14.0         tsne_0.1-3           mvtnorm_1.1-1       
[169] utf8_1.2.2           lattice_0.20-41      tibble_3.0.3        
[172] mixtools_1.2.0       numDeriv_2016.8-1.1  curl_4.3            
[175] gtools_3.8.2         zip_2.1.1            GO.db_3.11.4        
[178] survival_3.2-3       rmarkdown_2.14       munsell_0.5.0       
[181] DO.db_2.9            GetoptLong_1.0.2     iterators_1.0.13    
[184] haven_2.3.1          reshape2_1.4.4       gtable_0.3.0        
[187] rbibutils_2.0   

```


