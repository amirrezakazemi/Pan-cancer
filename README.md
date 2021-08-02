# Pan-cancer
This research project is about introducing a machine learning technique to cluster cancer samples based on their somatic mutation profiles. The problem is unsupervised clustering. We performed statistical biomarker analyses to evaluate our subtyping approach.
## Preprocessing
Gene annotation and obtaining mutation profile of samples. corresponding codes are in "preprocessing" folder

Codes: 
add_gene.R & add_genemotif.R: in thses two scripts we use Ensembl data to annotate DNAs and idenify start and end of each gene. Input: ICGC data, Output: A dataframe containing new columns.
count_genes.R & count_genemotifs.R: these two scripts are used to make a matrix that rows are samples and columns are genes/genemotifs. Input: ICGC data , Output: the matrix
health samples.R: to filter only patients
intgr_donor.R: preprocessing for clinical analysis.
summarize_data.R: to filter only somatic mutations in ICGC dataset
intEGRATE.R: preprocess ensemBL DATASET AND FILtER CODING GENES.

## Extract meaningful features.
The mutation profile of cancer samples is highly dispersing, so we need to identify genes that are significantly mutated and obtained genes are our candidate features. Since our project is on all cancer types we identified genes in each cancer separately. To do this, we performed a background statistical method (negative binomial distribution) to fit mutation data of each cancer. then based on a fixed P-value we obtained significantly mutated genes.corresponding codes are in "sgnf_genes" folder.

CODES:
DECIDE_DISTR.R: TO CHOOSE THE BEST FITTING DISTRIBUTION ON MUTATIONAL DATA FOR EACH CANCERTYPE
GENES_BG.R: FIT THE NEGATIVE BINOMIAL DISTRIBUTION AND USE THE THRESHOLD TO FILTER AND EXTRACT THE MOST SIGNIFICANT GENES IN EACH CANCERTYPE.
GENEMOTIFFS_BG.R: FIT THE BETA BINOMIAL distrIBUTion and use tHE ThresholD To filTEr and extrACT tHE MOST significant GENEMOTifs. after these two experiments we only use genes as our features. so the input for clustering process is a matrix containing number of mutations in each significant gene for each sample.


## Clustering
After obtaining features, we used a classical machine learning method to cluster samples. A hierarchical mixture model is performed in this step. We also examined other clustering methods such as k-means and db-scan. corresponding codes are in "clustering" folder.

codes:
clUSTRing1.R: tHIS CODE IS PERFORMING THE PROCESS OF CLUSTRING WITH KMEANS AND MCLUST
SUBCLUSTERS.R: THIS CODE IS FOR SETTING THE THRESHOLD OF 95 PERCENT. WE EXAMINED 2 THRESHOLD AND IDENTIFIED TWO SETS OF SUBTYPES AND THEN WE CHOSE THE BEST ONE.
CLusteR_valIDATion.r: performing stATistICAL analYSIS ON OUR SUBTypes tO VALidate The clUSTering performance by mclUST package.
clUSTer_analYSIS.R: CLusteR SAMPLEs using mclUST and plOTtING OUR SUBTypes and SAVE The resulT in bio_res folDER.


## Statistical biomarker analysis
In this step, we attempted to validate our clustering using biological interpretations. Different biomarker analyses such as enriched genes, enriched gene motifs, enriched gene ontologies, signature and survival analysis are performed and clusters had distinguished characteristics. corresponding codes are in "biomarker_analysis" folder.






