# T-cell-exhaustion-lung-breast-cancer-bioinformatics
In the present study, we discovered conserved biomarkers associated with T cell exhaustion in lung and breast cancer. By regulating these biomarkers, we may be able to treat resistance and create successful immunotherapies.

## Objectives:

-To identify key genes associated with T cell exhaustion in lung and breast cancers by constructing and analyzing protein–protein interaction networks.

-To evaluate the clinical relevance of these key genes by examining their association with patient survival and treatment response, and by determining whether their effects are shared across cancers or specific to individual cancer types.

-To identify key genes associated with disease progression using machine learning approaches such as Random Forest and LASSO, and evaluate their diagnostic potential as statistically supported candidate biomarkers for future immunotherapy research.
## Data Collection:
We downloaded eight microarray datasets from Gene Expression Omnibus (GEO). 

- Lung cancer: GSE10072, GSE18842, GSE19804, GSE32863  
- Breast cancer: GSE10780, GSE42568, GSE65194, GSE70947  

## Methods:
We summarize the key work we carried out as follows:

    • Identification of differentially expressed markers in eight independent lung and breast cancer cohorts.
    
    • ssGSEA-based exhaustion scoring and WGCNA were integrated to identify exhaustion-associated candidate genes, followed by network, survival, and machine-learning analyses to validate and prioritize biomarkers.
    
    • Additional exhaustion-related genes were identified, and their diagnostic potential was evaluated using LASSO regression and a random forest.

## Results
    • Our results revealed a conserved regulatory pathway that includes CCL5 and HAVCR2. This pathway links the initial recruitment of T cells to terminal exhaustion. We consistently identified significant immune checkpoint markers, including HAVCR2, PDCD1, and CTLA4, which support the reliability of our analysis.
    
    • Analysis of immune cell infiltration indicated that upregulated levels of FCGR2A and FCGR2B in myeloid cells and macrophages are associated with an immunosuppressive environment.
    
    • We propose that effective immunotherapeutic strategies should target both T cell exhaustion pathways and myeloid-mediated suppressive networks to restore durable antitumor immunity, offering a more effective treatment approach for both malignancies.

## Authors
- Nilofer Shaikh  
- Alankrita Srivastava  

## Code Availability

All scripts used for data preprocessing, analysis, and modeling are available in this repository for reproducibility.
