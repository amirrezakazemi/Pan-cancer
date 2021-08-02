# pan-cancer
this research project is about introducing a machine learning technique to cluster cancer samples based on their somatic mutation profiles. the problem is unsupervised clustering. we performed statistical biomarker analyses to evaluate our subtyping approach.
## preprocessing
gene annotation and obtaining mutation profile of samples. corresponding codes are in "preprocessing" folder

codes: 
add_gene.r & add_genemotif.r: in thses two scripts we use ensembl data to annotate dnas and idenify start and end of each gene. input: icgc data, output: a dataframe containing new columns.
count_genes.r & count_genemotifs.r: these two scripts are used to make a matrix that rows are samples and columns are genes/genemotifs. input: icgc data , output: the matrix
health samples.r: to filter only patients
intgr_donor.r: preprocessing for clinical analysis.
summarize_data.r: to filter only somatic mutations in icgc dataset
integrate.r: preprocess ensembl dataset and filter coding genes.

## extract meaningful features.
the mutation profile of cancer samples is highly dispersing, so we need to identify genes that are significantly mutated and obtained genes are our candidate features. since our project is on all cancer types we identified genes in each cancer separately. to do this, we performed a background statistical method (negative binomial distribution) to fit mutation data of each cancer. then based on a fixed p-value we obtained significantly mutated genes.corresponding codes are in "sgnf_genes" folder.

codes:
decide_distr.r: to choose the best fitting distribution on mutational data for each cancertype
genes_bg.r: fit the negative binomial distribution and use the threshold to filter and extract the most significant genes in each cancertype.
genemotiffs_bg.r: fit the beta binomial distribution and use the threshold to filter and extract the most significant genemotifs. after these two experiments we only use genes as our features. so the input for clustering process is a matrix containing number of mutations in each significant gene for each sample.


## clustering
after obtaining features, we used a classical machine learning method to cluster samples. a hierarchical mixture model is performed in this step. we also examined other clustering methods such as k-means and db-scan. corresponding codes are in "clustering" folder.

codes:
clustring1.r: this code is performing the process of clustring with kmeans and mclust
subclusters.r: this code is for setting the threshold of 95 percent. we examined 2 threshold and identified two sets of subtypes and then we chose the best one.
cluster_validation.r: performing statistical analysis on our subtypes to validate the clustering performance by mclust package.
cluster_analysis.r: cluster samples using mclust and plotting our subtypes and save the result in bio_res folder.


## statistical biomarker analysis
in this step, we attempted to validate our clustering using biological interpretations. different biomarker analyses such as enriched genes, enriched gene motifs, enriched gene ontologies, signature and survival analysis are performed and clusters had distinguished characteristics. corresponding codes are in "biomarker_analysis" folder.






